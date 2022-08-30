#!/bin/bash

set -e

ml Trimmomatic/0.39-Java-11
ml BEDTools/2.29.2-GCC-9.3.0
ml SAMtools/1.10-GCCcore-8.3.0
ml FastQC/0.11.9-Java-11
ml picard/2.21.6-Java-11
ml Python
ml Pysam
ml GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8
# ROOT is the output directory
root=${ROOT}
# CELLRANGER_FOLDER containing MiSeq cellranger count output folders, named by sample names
raw_data_folder=${CELLRANGER_MKFASTQ_FOLDER}
# CELLRANGER_FOLDER_NOVA containing GEX cellranger count output folders (not 16s), named by sample names
nova_bam_path=${CELLRANGER_FOLDER_NOVA}
# cellranger database
cellrangerdb=${CELLRANGER_DB}
# Pathseq database directory
pathseqdb=${PATHSEQ_DB}

workdir=${root}
cd ${workdir}

mkdir cellranger_count
cd cellranger_count
# Cellranger process
for folder in ${raw_data_folder}/*
do
folder=${folder}
folder_name=${folder##*/}
path=${folder}
echo ${path}
echo ${folder_name}
cellranger count \
--id=${folder_name} \
--transcriptome=${cellrangerdb} \
--fastqs=${path} \
--sample=${folder_name}
done
# splid the bam file to read 1 and read 2
cd ${workdir}
mkdir split_reads
cd split_reads
for folder in ${raw_data_folder}/*
do
folder_name=${folder##*/}
file=${folder}/outs/possorted_genome_bam.bam
echo ${file}
samplename=${folder_name}
folder_name=${folder##*/}
file=${folder}/outs/possorted_genome_bam.bam
samplename=${folder_name}
bedtools bamtofastq -i ${file} \
                      -fq ${samplename}.r1.fq \
                      -fq2 ${samplename}.r2.fq
done
# fastQC
mkdir ${root}/preqc
fastqc \
-o ${root}/preqc \
-t 36 \
${root}/split_reads/*.fq
# trimmomatic process
cd ${root}/split_reads
mkdir trim
for str in *r1.fq
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE \
-threads 36 \
${str} \
trim/${str}.SE_trim.fq \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:15
done
# fastQC
cd trim
mkdir ${root}/postqc
fastqc \
-o ${root}/postqc \
-t 36 *.SE_trim.fq

mkdir ${workdir}/ubams_r1
cd ${root}/split_reads/trim
# convert fastq files to bam file
for file in *SE_trim.fq
do
java -Xmx750G -jar $EBROOTPICARD/picard.jar FastqToSam \
    FASTQ=${file} \
    OUTPUT=${file}.bam \
    READ_GROUP_NAME=rg1 \
    SAMPLE_NAME=cellranger1 

mv ${file}.bam ${workdir}/ubams_r1
done

folder=${workdir}/ubams_r1
outpath=${workdir}/pathseq_r1
mkdir ${outpath}
cd ${folder}

# PathSeq process
for each_file in *.bam
do
filename="${each_file%.*}"
filename="${filename%.*}"
filename="${filename%.*}"
samplename=$filename
echo $filename
gatk --java-options "-Xmx750g" PathSeqPipelineSpark \
    --input ${each_file} \
    --filter-bwa-image ${pathseqdb}/pathseq_host.fa.img \
    --kmer-file ${pathseqdb}/pathseq_host.bfi \
    --min-clipped-read-length 60 \
    --microbe-fasta ${pathseqdb}/pathseq_microbe.fa \
    --microbe-bwa-image ${pathseqdb}/pathseq_microbe.fa.img \
    --taxonomy-file ${pathseqdb}/pathseq_taxonomy.db \
    --output ${outpath}/${samplename}.pathseq.complete.bam \
    --scores-output ${outpath}/${samplename}.pathseq.complete.txt.csv \
    --is-host-aligned false \
    --filter-duplicates false \
    --min-score-identity .7
done

# generate metadata files using UMI_annotator.py (genus.csv is the metadata matrix)
bam_path=${raw_data_folder}
pathseq_path=${workdir}/pathseq_r1
out_path=${workdir}/python
mkdir ${out_path}
cd ${bam_path}
# cell names will be: sample_${each_sample}
for each_sample in *
do
echo ${each_sample}
python python UMI_annotator.py \
${bam_path}/${each_sample}/outs/possorted_genome_bam.bam \
Sample_${each_sample} \
${nova_bam_path}/${each_sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
${pathseq_path}/${each_sample}.r1.fq.pathseq.complete.bam \
${pathseq_path}/${each_sample}.r1.fq.pathseq.complete.txt.csv \
${out_path}/${each_sample}.mi.filtered_matrix.readname \
${out_path}/${each_sample}.mi.filtered_matrix.unmap_cbub.bam \
${out_path}/${each_sample}.mi.filtered_matrix.unmap_cbub.fasta \
${out_path}/${each_sample}.mi.filtered_matrix.list \
${out_path}/${each_sample}.mi.raw.filtered_matrix.readnamepath \
${out_path}/${each_sample}.mi.filtered_matrix.genus.cell \
${out_path}/${each_sample}.mi.filtered_matrix.genus.csv \
${out_path}/${each_sample}.mi.filtered_matrix.validate.csv
done

