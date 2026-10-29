suppressMessages(library(DESeq2))
suppressMessages(library(pheatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(getopt))

spec <- matrix(c(
    "in_dir",       "i", 1, "character", "expression matrix file",
    "out_dir",      "o", 1, "character", "output directory",
    "case_control", "c", 1, "character", "case control csv",
    "group",        "g", 1, "character", "group csv",
    "orgdb",        "d", 1, "character", "OrgDb package name",
    "pvalue",       "p", 2, "double",    "padj threshold (default 0.05)",
    "log2fc",       "f", 2, "double",    "log2FC threshold (default 1)",
    "mincounts",    "m", 2, "integer",   "min count (default 2)",
    "topgene",      "t", 2, "integer",   "top gene number (default 50)"
), byrow=TRUE, ncol=5)
args <- getopt(spec)

in_dir <- args$in_dir
out_dir <- args$out_dir
case_control <- args$case_control
group_file <- args$group
OrgDb <- args$orgdb
pValue <- ifelse(is.null(args$pvalue), 0.05, args$pvalue)
log2FC <- ifelse(is.null(args$log2fc), 1, args$log2fc)
minCounts <- ifelse(is.null(args$mincounts), 2, args$mincounts)
topGene <- ifelse(is.null(args$topgene), 50, args$topgene)

dir.create(out_dir, showWarnings=FALSE)

draw_volcano <- function(data, pval, dir) {
    data$group <- ifelse(data$log2FoldChange > log2FC & data$pvalue < pval, 'up',
                        ifelse(data$log2FoldChange < -log2FC & data$pvalue < pval, 'down', 'none'))
    up <- subset(data, group=='up')
    down <- subset(data, group=='down')

    p <- ggplot(data, aes(log2FoldChange, -log10(pvalue), color=group)) +
        geom_point(alpha=0.8, size=1) +
        scale_color_manual(values=c('down'="#419A18", 'none'="gray", 'up'="#EE4E4A")) +
        geom_vline(xintercept=c(-log2FC, log2FC), lty=4) +
        geom_hline(yintercept=-log10(pval), lty=4) +
        labs(x="log2(FoldChange)", y="-log10(pvalue)",
             title=paste("up:", nrow(up), "down:", nrow(down))) +
        theme_bw()

    if(nrow(up) > 0) p <- p + geom_label_repel(data=head(up[order(up$pvalue),], 5), 
                                                aes(label=SYMBOL), size=2.5)
    if(nrow(down) > 0) p <- p + geom_label_repel(data=head(down[order(down$pvalue),], 5), 
                                                  aes(label=SYMBOL), size=2.5)

    ggsave(file.path(dir, "volcano.pdf"), p, width=7, height=6)
    ggsave(file.path(dir, "volcano.png"), p, width=7, height=6)
    return(data)
}

draw_heatmap <- function(vsd, dds, data, dir, topn) {
    select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:min(topn, nrow(dds))]
    mat <- assay(vsd)[select, ]
    df <- as.data.frame(colData(dds)[, "condition", drop=FALSE])
    pheatmap(mat, annotation_col=df, show_rownames=FALSE, 
             color=colorRampPalette(c("orange", "white", "purple"))(100),
             filename=file.path(dir, "heatmap_top.pdf"))

    diff_genes <- subset(data, group %in% c('up','down'))
    if(nrow(diff_genes) > 1) {
        diff_mat <- assay(vsd)[diff_genes$ENSEMBL[diff_genes$ENSEMBL %in% rownames(assay(vsd))], ]
        if(nrow(diff_mat) > 1) {
            pheatmap(diff_mat, scale='row', annotation_col=df, show_rownames=FALSE,
                     color=colorRampPalette(c("orange", "white", "purple"))(100),
                     filename=file.path(dir, "heatmap_diff.pdf"))
        }
    }
}

draw_pca <- function(vsd, dir) {
    pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
        geom_point(size=3) +
        geom_text_repel(aes(label=name), show.legend=FALSE) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_bw()
    ggsave(file.path(dir, "pca.pdf"), p, width=6, height=5)
    ggsave(file.path(dir, "pca.png"), p, width=6, height=5)
}

run_deseq2 <- function(counts, condition_df, case, control, dir) {
    dir.create(dir, showWarnings=FALSE)

    dds <- DESeqDataSetFromMatrix(countData=counts, colData=condition_df, design=~condition)
    keep <- rowSums(counts(dds)) >= minCounts
    dds <- dds[keep,]
    dds$condition <- relevel(dds$condition, ref=control)
    dds <- DESeq(dds)

    res <- results(dds, contrast=c("condition", case, control))
    res <- res[order(res$pvalue),]

    vsd <- vst(dds, blind=FALSE)
    write.csv(assay(vsd), file.path(dir, "count_vst.csv"), quote=FALSE)

    res_df <- as.data.frame(res)
    res_df$ENSEMBL <- rownames(res_df)
    anno <- bitr(res_df$ENSEMBL, fromType="ENSEMBL", toType=c("SYMBOL","ENTREZID"), OrgDb=OrgDb)
    res_df <- merge(res_df, anno, by="ENSEMBL", all.x=TRUE)

    res_df <- draw_volcano(res_df, pValue, dir)
    write.csv(res_df, file.path(dir, "diff_genes.csv"), row.names=FALSE, quote=FALSE)

    draw_heatmap(vsd, dds, res_df, dir, topGene)
    draw_pca(vsd, dir)

    sig <- subset(res_df, padj < pValue & abs(log2FoldChange) > log2FC)
    write.csv(sig, file.path(dir, paste0("sig_genes_padj", pValue, ".csv")), row.names=FALSE, quote=FALSE)

    return(list(up=sum(res_df$group=='up', na.rm=TRUE), down=sum(res_df$group=='down', na.rm=TRUE)))
}

main <- function() {
    counts <- read.csv(in_dir, row.names=1, check.names=FALSE)
    groups <- read.csv(group_file, colClasses=c("character","character"), check.names=FALSE)
    cc <- read.csv(case_control, colClasses=c("character","character"), check.names=FALSE)

    counts <- counts[rowSums(counts) > 0,]

    all_stats <- data.frame()

    for(i in 1:nrow(cc)) {
        case <- cc$case[i]
        control <- cc$control[i]
        cat("Processing:", case, "vs", control, "\n")

        samples <- subset(groups, group %in% c(case, control))
        condition_df <- data.frame(sample=samples$sample, condition=factor(samples$group, levels=c(control, case)))
        counts_sub <- counts[, condition_df$sample]

        dir <- file.path(out_dir, paste0(case, "_VS_", control))
        stats <- run_deseq2(counts_sub, condition_df, case, control, dir)

        all_stats <- rbind(all_stats, data.frame(comparison=paste0(case,"_VS_",control), 
                                                  up=stats$up, down=stats$down))
    }

    write.csv(all_stats, file.path(out_dir, "summary.csv"), row.names=FALSE, quote=FALSE)
}

main()
