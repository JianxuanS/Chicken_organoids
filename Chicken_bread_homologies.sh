# The Ensembl 110 layer genome lacks GO information 
# To perform enrichment analysis, we need to map the layer genes to broiler genes. Using the Ensembl rapid release, approximately 95% of the layer genome has corresponding broiler homologs.
wget https://ftp.ensembl.org/pub/rapid-release/species/Gallus_gallus/GCA_016700215.2/ensembl/homology/2022_01/Gallus_gallus-GCA_016700215.2-2022_01-homology.tsv.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/gallus_gallus_gca016700215v2/Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.113.gtf.gz
gunzip *.gz

grep -w "broiler" /exports/eddie/scratch/s1964742/Gallus_gallus-GCA_016700215.2-2022_01-homology.tsv > broiler_homologies.tsv
awk 'BEGIN {OFS="\t"} $7 == "homolog_rbbh" || $7 == "homolog_bbh" {$7 = $6} {print}' broiler_homologies.tsv > modified_broiler_homologies.tsv
awk '{print $5, $6, $7}' modified_broiler_homologies.tsv > broiler_homologies_columns_567.tsv

grep -w "protein_coding" "/exports/eddie/scratch/s1964742/Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.113.gtf" | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[5]"\t"$7}' | sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_biotype "//'| sed 's/gene_name "//' | sed 's/gene_biotype "//' | sed 's/"//g' | sed 's/ //g' | sed '1igene_id\tGeneSymbol\tChromosome\tClass\tStrand' > chicken_annotation_table.tsv
awk 'BEGIN {OFS="\t"} $2 == "gene_sourceensembl" {$2 = $1} {print}' chicken_annotation_table.tsv > modified_chicken_annotation_table.tsv
awk '{print $1,$2}' modified_chicken_annotation_table.tsv > final_chicken_annotation_table.tsv