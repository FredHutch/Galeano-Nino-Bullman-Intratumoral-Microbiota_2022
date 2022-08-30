#!/bin/bash

set -e

# The preprocessing pipeline for single cell NovaSeq data
ml CellRanger/6.1.1
ml GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
ml Python
ml Pysam

raw_data_folder=${CELLRANGER_MKFASTQ_FOLDER} # the folder containing Cellranger mkfastq output folders
root=${WORKING_DIR} # working directory
pathseqdb=${PATHSEQ_DB} # Pathseq database
cellrangerdb=${CELLRANGER_DB}

cd ${root}
mkdir cellranger_count
cd cellranger_count

# Cellranger count processing
for folder in ${raw_data_folder}/*
do
folder=${folder}
folder_name=${folder##*/}
path=${folder}
echo ${path}
cellranger count \
--id=${folder_name} \
--transcriptome=${cellrangerdb} \
--fastqs=${path} \
--sample=${folder_name}
done

outpath=${root}
mkdir ${outpath}
outpath=${outpath}/pathseq
mkdir ${outpath}
# PathSeq process
for folder in *
do
folder_name=${folder##*/}
file=${folder}/outs/possorted_genome_bam.bam
samplename=${folder_name}
echo ${samplename}
gatk --java-options "-Xmx750g" PathSeqPipelineSpark \
    --input ${file} \
    --filter-bwa-image ${pathseqdb}/pathseq_host.fa.img \
    --kmer-file ${pathseqdb}/pathseq_host.bfi \
    --min-clipped-read-length 60 \
    --microbe-fasta ${pathseqdb}/pathseq_microbe.fa \
    --microbe-bwa-image ${pathseqdb}/pathseq_microbe.fa.img \
    --taxonomy-file ${pathseqdb}/pathseq_taxonomy.db \
    --output ${outpath}/${samplename}.pathseq.complete.bam \
    --scores-output ${outpath}/${samplename}.pathseq.complete.csv \
    --is-host-aligned false \
    --filter-duplicates false \
    --min-score-identity .7
done

bam_path=${root}/cellranger_count
pathseq_path=${root}/pathseq
out_path=${root}/python
mkdir ${out_path}
cd ${bam_path}
# generate metadata files using UMI_annotator.py (genus.csv is the metadata matrix)
for each_sample in *
do
echo ${each_sample}
# cell names will be: sample_${each_sample}
python UMI_annotator.py \
${bam_path}/${each_sample}/outs/possorted_genome_bam.bam \
sample_${each_sample} \
${bam_path}/${each_sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
${pathseq_path}/${each_sample}.pathseq.complete.bam \
${pathseq_path}/${each_sample}.pathseq.complete.csv \
${out_path}/${each_sample}.nova.filtered_matrix.readname \
${out_path}/${each_sample}.nova.filtered_matrix.unmap_cbub.bam \
${out_path}/${each_sample}.nova.filtered_matrix.unmap_cbub.fasta \
${out_path}/${each_sample}.nova.filtered_matrix.list \
${out_path}/${each_sample}.nova.raw.filtered_matrix.readnamepath \
${out_path}/${each_sample}.nova.filtered_matrix.genus.cell \
${out_path}/${each_sample}.nova.filtered_matrix.genus.csv \
${out_path}/${each_sample}.nova.filtered_matrix.validate.csv
done

