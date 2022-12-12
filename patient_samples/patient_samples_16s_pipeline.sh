#!/bin/bash
ml CellRanger/6.1.1

ml BEDTools/2.29.2-GCC-9.3.0
ml SAMtools/1.16.1-GCC-11.2.0
ml FastQC/0.11.9-Java-11
ml Trimmomatic/0.39-Java-11
ml picard/2.21.6-Java-11


ml GATK/4.1.3.0-GCCcore-8.3.0-Java-1.8

ml Python
ml Pysam

`raw_data_folder` # the folder containing Cellranger mkfastq output folders
`root` # working directory
`pathseqdb` # Pathseq database
`cellrangerdb` # Cellranger database
`gex_bam_path` # barcodes.tsv.gz from GEX is used as a 'whitelist' for real cells

root=${workdir}
# Run cellranger count
cd ${workdir}
mkdir cellranger_count
cd cellranger_count
for folder in ${raw_data_folder}/*
do
folder_name=${folder##*/}
path=${folder}
cellranger count \
--id=${folder_name} \
--transcriptome=${cellrangerdb} \
--fastqs=${path} \
--sample=${folder_name}
done

cd ${workdir}
mkdir split_reads
cd split_reads

# convert cellranger bam file to fastqs
for folder in ${workdir}/cellranger_count/*
do
folder_name=${folder##*/}
file=${folder}/outs/possorted_genome_bam.bam
echo ${file}

samplename=${folder_name}
bedtools bamtofastq -i ${folder}/outs/possorted_genome_bam.bam \
                      -fq ${samplename}.r1.fq \
                      -fq2 ${samplename}.r2.fq
done

# run fastqc before trimmomatic
mkdir ${root}/preqc
fastqc \
-o ${root}/preqc \
${root}/split_reads/*.fq

# run trimmomatic on R1 
cd ${root}/split_reads
mkdir trim
for str in *r1.fq
do
# adjust -threads to number of cores you would like to use
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE \
-threads 36 \
${str} \
trim/${str}.SE_trim.fq \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:15
done

# run fastqc after trimmomatic 
cd trim
mkdir ${root}/postqc
fastqc \
-o ${root}/postqc \
*.SE_trim.fq

mkdir ${workdir}/ubams_r1
cd ${root}/split_reads/trim

# convert R1 to ubam file in order to run Pathseq
for file in *SE_trim.fq
do
java -Xmx700G -jar $EBROOTPICARD/picard.jar FastqToSam \
    FASTQ=${file} \
    OUTPUT=${file}.bam \
    READ_GROUP_NAME=16s \
    SAMPLE_NAME=16s

# move and rename generated ubam files
mv ${file}.bam ${workdir}/ubams_r1
done

ubam_folder=${workdir}/ubams_r1
outpath=${workdir}/pathseq_r1
mkdir ${outpath}

cd ${ubam_folder}

# Pathseq to identify pathogen-associated cells
for each_file in *.bam
do
echo ${each_file}
filename="${each_file%.*}"
filename="${filename%.*}"
filename="${filename%.*}"
samplename=${filename}
echo ${samplename}
# PathSeq process # Please adjust "-Xmx750g" based on the memory you want to use. Adjust --min-score-identity and --min-clipped-read-length based on your samples
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

# Python script to produce a bacteria UMI matrix (based on valid GEX cell)
bam_path=${workdir}/cellranger_count
pathseq_path=${workdir}/pathseq_r1
out_path=${root}/python
mkdir ${out_path}
cd ${bam_path}
# barcodes.tsv.gz from GEX is used as a 'whitelist' for real cells. 
for each_sample in *
do
echo ${each_sample}
echo ${gex_bam_path}
python INVADEseq.py \
${bam_path}/${each_sample}/outs/possorted_genome_bam.bam \
${each_sample} \
${gex_bam_path}/${each_sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
${pathseq_path}/${each_sample}.r1.fq.pathseq.complete.bam \
${pathseq_path}/${each_sample}.r1.fq.pathseq.complete.txt.csv \
${out_path}/${each_sample}.16s.filtered_matrix.readname \
${out_path}/${each_sample}.16s.filtered_matrix.unmap_cbub.bam \
${out_path}/${each_sample}.16s.filtered_matrix.unmap_cbub.fasta \
${out_path}/${each_sample}.16s.filtered_matrix.list \
${out_path}/${each_sample}.16s.raw.filtered_matrix.readnamepath \
${out_path}/${each_sample}.16s.filtered_matrix.genus.cell \
${out_path}/${each_sample}.16s.filtered_matrix.genus.csv \
${out_path}/${each_sample}.16s.filtered_matrix.validate.csv
done



