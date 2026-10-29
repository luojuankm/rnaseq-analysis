#!/bin/bash
set -e

in_dir=$1
rsem_index=$2
threads=$3
sample=$4

mkdir -p ${in_dir}/06.rsem_expression/${sample}

bam_file=${in_dir}/12.star/${sample}/star_${sample}Aligned.toTranscriptome.out.bam

rsem-calculate-expression --bam \
    --estimate-rspd \
    --calc-ci \
    --paired-end \
    -p ${threads} \
    --no-bam-output \
    --forward-prob 0.5 \
    ${bam_file} \
    ${rsem_index} \
    ${in_dir}/06.rsem_expression/${sample}/rsem_expression_${sample}

genes_results=${in_dir}/06.rsem_expression/${sample}/rsem_expression_${sample}.genes.results
isoforms_results=${in_dir}/06.rsem_expression/${sample}/rsem_expression_${sample}.isoforms.results

awk -F'\t' '{print $1 "," $5}' ${genes_results} > ${in_dir}/06.rsem_expression/${sample}/rsem_expression_count_${sample}.csv
awk -F'\t' '{print $1 "," $6}' ${genes_results} > ${in_dir}/06.rsem_expression/${sample}/rsem_expression_tpm_${sample}.csv
awk -F'\t' '{print $1 "," $7}' ${genes_results} > ${in_dir}/06.rsem_expression/${sample}/rsem_expression_fpkm_${sample}.csv

awk -F'\t' '{print $1 "," $2 "," $5}' ${isoforms_results} > ${in_dir}/06.rsem_expression/${sample}/rsem_transcript_count_${sample}.csv
awk -F'\t' '{print $1 "," $2 "," $6}' ${isoforms_results} > ${in_dir}/06.rsem_expression/${sample}/rsem_transcript_tpm_${sample}.csv
awk -F'\t' '{print $1 "," $2 "," $7}' ${isoforms_results} > ${in_dir}/06.rsem_expression/${sample}/rsem_transcript_fpkm_${sample}.csv
