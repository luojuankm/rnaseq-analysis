#!/bin/bash
set -e

in_dir=$1
star_index=$2
threads=$3
sample=$4

mkdir -p ${in_dir}/12.star/${sample}

R1=$(ls ${in_dir}/01.CleanData/${sample}/${sample}*_clean_R1.f*q.gz 2>/dev/null || ls ${in_dir}/01.CleanData/${sample}/${sample}*_clean_1*.f*q.gz 2>/dev/null | head -1)
R2=$(ls ${in_dir}/01.CleanData/${sample}/${sample}*_clean_R2.f*q.gz 2>/dev/null || ls ${in_dir}/01.CleanData/${sample}/${sample}*_clean_2*.f*q.gz 2>/dev/null | head -1)

if [ -n "$R1" ] && [ -n "$R2" ]; then
    read_type="--readFilesIn ${R1} ${R2}"
elif [ -n "$R1" ]; then
    read_type="--readFilesIn ${R1}"
else
    echo "Error: No input files found for ${sample}"
    exit 1
fi

STAR --genomeDir ${star_index} \
    ${read_type} \
    --runThreadN ${threads} \
    --outFilterMultimapNmax 10 \
    --alignSJoverhangMin 5 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 2 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --sjdbScore 1 \
    --readFilesCommand zcat \
    --outFileNamePrefix ${in_dir}/12.star/${sample}/star_${sample} \
    --limitBAMsortRAM 40000000000

samtools index ${in_dir}/12.star/${sample}/star_${sample}Aligned.sortedByCoord.out.bam
