#!/bin/bash
#
# =============================================================================
# Setup Instructions
# =============================================================================
#
# Grid Engine options (lines prefixed with #$)
#  job name: -N
#  use the current working directory: -cwd
#  number of cores -pe sharedmem
#  runtime limit: -l h_rt
#  memory limit: -l h_vmem
#$ -N STARsolo_LV_3
#$ -cwd
#$ -pe sharedmem 8
#$ -l h_rt=48:00:00
#$ -l h_vmem=20G
#$ -P roslin_macqueen_lab

# Make directory before STAR
#
mkdir /exports/eddie/scratch/s1964742/Chicken_filtered_mat.layer/LV_3
cd /exports/eddie/scratch/s1964742/Chicken_filtered_mat.layer/LV_3

# Module load
#
wget https://github.com/alexdobin/STAR/releases/download/2.7.10a_alpha_220601/STAR_2.7.10a_alpha_220601_Linux_x86_64_static.zip
unzip STAR_2.7.10a_alpha_220601_Linux_x86_64_static.zip

# Running mapping jobs
#
./STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir /exports/eddie/scratch/s1964742/Chicken_filtered_mat.layer \
--sjdbOverhang 99 \
--readFilesIn \
"/exports/eddie/scratch/s1964742/chicken_organoid/LV_230323/outs/fastq_path/AAAYWGVHV/LV_3_S3_L001_R2_001.fastq.gz",\
"/exports/eddie/scratch/s1964742/chicken_organoid/LV_230323/outs/fastq_path/AAAYWGVHV/LV_3_S3_L002_R2_001.fastq.gz",\
"/exports/eddie/scratch/s1964742/chicken_organoid/LV_240323/outs/fastq_path/AAAYW7HHV/LV_3_S3_L001_R2_001.fastq.gz",\
"/exports/eddie/scratch/s1964742/chicken_organoid/LV_240323/outs/fastq_path/AAAYW7HHV/LV_3_S3_L002_R2_001.fastq.gz" \
"/exports/eddie/scratch/s1964742/chicken_organoid/LV_230323/outs/fastq_path/AAAYWGVHV/LV_3_S3_L001_R1_001.fastq.gz",\
"/exports/eddie/scratch/s1964742/chicken_organoid/LV_230323/outs/fastq_path/AAAYWGVHV/LV_3_S3_L002_R1_001.fastq.gz",\
"/exports/eddie/scratch/s1964742/chicken_organoid/LV_240323/outs/fastq_path/AAAYW7HHV/LV_3_S3_L001_R1_001.fastq.gz",\
"/exports/eddie/scratch/s1964742/chicken_organoid/LV_240323/outs/fastq_path/AAAYW7HHV/LV_3_S3_L002_R1_001.fastq.gz" \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMmode Full \
--outSAMattributes GX CB UB \
--soloType CB_UMI_Simple \
--soloCBwhitelist /exports/eddie/scratch/s1964742/3M-february-2018.txt \
--soloFeatures GeneFull_Ex50pAS Velocyto \
--soloUMIdedup 1MM_CR \
--soloCBstart 1 \
--soloCBlen 16 \
--soloUMIstart 17 \
--soloUMIlen 12 \
--soloBarcodeReadLength 0 \
--soloUMIfiltering MultiGeneUMI_CR \
--soloCellFilter EmptyDrops_CR \
--soloMultiMappers EM \
--limitOutSJcollapsed 2500000 \
--soloCellReadStats Standard

rm STAR_2.7.10a_alpha_220601_Linux_x86_64_static.zip
rm STAR

