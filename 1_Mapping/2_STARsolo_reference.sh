#!/bin/bash
#  job name: -N
#  use the current working directory: -cwd
#  number of cores -pe sharedmem
#  runtime limit: -l h_rt
#  memory limit: -l h_vmem
#$ -N STARsolo_reference
#$ -cwd
#$ -pe sharedmem 8
#$ -l h_rt=01:00:00
#$ -l h_vmem=20G
#$ -l rl9=false
#$ -P roslin_macqueen_lab
#$ -e error.txt
#$ -o output.txt
#$ -m beas
#$ -M J.Sun-52@sms.ed.ac.uk

# Make directory before STAR
#
mkdir /exports/eddie/scratch/s1964742/Chicken_filtered_mat.layer
cd /exports/eddie/scratch/s1964742

# Module load
#
source /etc/profile
module load roslin/star/2.7.10a
module load igmm/apps/cellranger/6.0.2

# Download fasta and gtf file
#
wget https://ftp.ensembl.org/pub/release-110/fasta/gallus_gallus_gca016700215v2/dna/Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/gallus_gallus_gca016700215v2/Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.110.gtf.gz
gunzip *gz

# Filtered gtf - protein_coding
#
cellranger mkgtf \
./Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.110.gtf \
./Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.110.filtered.gtf \
--attribute=gene_biotype:protein_coding

# Generating genome indexes
#
STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /exports/eddie/scratch/s1964742/Chicken_filtered_mat.layer/ \
--genomeFastaFiles ./Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.dna.toplevel.fa \
--sjdbGTFfile ./Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.110.filtered.gtf \
--sjdbOverhang 99 \
--genomeChrBinNbits 15

rm ./Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.110.gtf
rm ./Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.dna.toplevel.fa
rm ./Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.110.filtered.gtf
