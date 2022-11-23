#!/bin/bash

set -e

ml GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
ml Krona/2.7.1-GCCcore-9.3.0-Perl-5.30.2
ml Python
ml Pysam

# ROOT is the output directory
root=${ROOT}
# SPACERANGER_FOLDER containing SpaceRanger output folders, named by sample names
bam_path=${SPACERANGER_FOLDER}
# Pathseq database directory
pathseqdb=${PATHSEQ_DB}

cd ${bam_path}
outpath=${root}
mkdir ${outpath}
outpath=${outpath}/pathseq
mkdir ${outpath}
# Pathseq process
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

pathseq_path=${root}/pathseq
out_path=${root}/python
mkdir ${out_path}
cd ${bam_path}

# generate metadata files using UMI_annotator.py (genus.csv is the metadata matrix)
for folder in *
do
each_sample=${folder}
echo ${each_sample}
python UMI_annotator_Visium.py \
${bam_path}/${each_sample}/outs/possorted_genome_bam.bam \
'' \
${bam_path}/${each_sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz \
${pathseq_path}/${each_sample}.pathseq.complete.bam \
${pathseq_path}/${each_sample}.pathseq.complete.csv \
${out_path}/${each_sample}.visium.raw_matrix.readname \
${out_path}/${each_sample}.visium.raw_matrix.unmap_cbub.bam \
${out_path}/${each_sample}.visium.raw_matrix.unmap_cbub.fasta \
${out_path}/${each_sample}.visium.raw_matrix.list \
${out_path}/${each_sample}.visium.raw.raw_matrix.readnamepath \
${out_path}/${each_sample}.visium.raw_matrix.genus.cell \
${out_path}/${each_sample}.visium.raw_matrix.genus.csv \
${out_path}/${each_sample}.visium.raw_matrix.validate.csv
done

