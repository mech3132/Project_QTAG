#!bin/bash

mkdir -p 03a_assign_taxonomy

conda activate qiime2-2022.8

# Import sequences and tables
qiime tools import \
--input-path 02d_extract_sequences/repset.fasta \
--output-path 03a_assign_taxonomy/repset.qza \
--type 'FeatureData[Sequence]'

## 16S
biom convert --to-hdf5 \
-i 02b_adjust_files_for_processing/otu_adj_16S.txt \
-o 03a_assign_taxonomy/otu_16S.biom \
--table-type 'OTU table'

qiime tools import \
--input-path 03a_assign_taxonomy/otu_16S.biom \
--output-path 03a_assign_taxonomy/otu_16S.qza \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format 

## 18S
biom convert --to-hdf5 \
-i 02b_adjust_files_for_processing/otu_adj_18S.txt \
-o 03a_assign_taxonomy/otu_18S.biom \
--table-type 'OTU table'

qiime tools import \
--input-path 03a_assign_taxonomy/otu_18S.biom \
--output-path 03a_assign_taxonomy/otu_18S.qza \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format 


########## Filter seqs by 16S/18S #############
qiime feature-table filter-seqs \
  --i-data 03a_assign_taxonomy/repset.qza \
  --i-table 03a_assign_taxonomy/otu_16S.qza \
  --o-filtered-data 03a_assign_taxonomy/repset_16S.qza
  
qiime feature-table filter-seqs \
  --i-data 03a_assign_taxonomy/repset.qza \
  --i-table 03a_assign_taxonomy/otu_18S.qza \
  --o-filtered-data 03a_assign_taxonomy/repset_18S.qza



## Get databases
#if [[ -f 03a_assign_taxonomy/silva-138-99-seqs.qza ]]; then
# echo "silva seqs already downloaded"
#else
#wget \
#  -O "03a_assign_taxonomy/silva-138-99-seqs.qza" \
#  "https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza"
#fi

#if [[ -f 03a_assign_taxonomy/silva-138-99-tax.qza ]]; then
# echo "silva seqs already downloaded"
#else
#wget \
#  -O "03a_assign_taxonomy/silva-138-99-tax.qza" \
#  "https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza"
#fi

### Assign taxonomy
####### Done on compute canada cedar

#SBATCH --account=def-haney
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=melissa.c1010@gmail.com
#SBATCH --mail-type=ALL

#module load singularity

##wget \
##  -O "/inputs/silva-138-99-seqs.qza" \
##  "https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza"
##
##wget \
##  -O "/inputs/silva-138-99-tax.qza" \
##  "https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza"
#  
#mkdir -p  outputs/qtag

#singularity exec -B $PWD:/home -B ~/scratch/outputs:/outputs \
#  -B ~/scratch/inputs:/inputs -B ~/scratch/databases:/databases \
# ./qiime2-2022.8.sif \
#  qiime feature-classifier classify-consensus-vsearch \
#--i-query /inputs/qtag/repset_16S.qza \
#--i-reference-reads /databases/silva-138-99-seqs.qza \
#--i-reference-taxonomy /databases/silva-138-99-tax.qza \
#--p-perc-identity 0.8 \
#--p-strand 'both' \
#--p-threads 10 \
#--output-dir /outputs/qtag/vsearch_classify_id_80_taxonomy16

#singularity exec -B $PWD:/home -B ~/scratch/outputs:/outputs \
#  -B ~/scratch/inputs:/inputs -B ~/scratch/databases:/databases \
# ./qiime2-2022.8.sif \
#  qiime feature-classifier classify-consensus-vsearch \
#--i-query /inputs/qtag/repset_18S.qza \
#--i-reference-reads /databases/silva-138-99-seqs.qza \
#--i-reference-taxonomy /databases/silva-138-99-tax.qza \
#--p-perc-identity 0.8 \
#--p-strand 'both' \
#--p-threads 10 \
#--output-dir /outputs/qtag/vsearch_classify_id_80_taxonomy18


qiime tools export \
--input-path 03a_assign_taxonomy/vsearch_classify_id_80_taxonomy16/classification.qza \
--output-path 03a_assign_taxonomy/vsearch_classify_id_80_taxonomy16/export

qiime tools export \
--input-path 03a_assign_taxonomy/vsearch_classify_id_80_taxonomy18/classification.qza \
--output-path 03a_assign_taxonomy/vsearch_classify_id_80_taxonomy18/export


###### Filter out non-chloro and non-mito #########

qiime taxa filter-table \
--i-table 03a_assign_taxonomy/otu_16S.qza \
--i-taxonomy 03a_assign_taxonomy/vsearch_classify_id_80_taxonomy16/classification.qza \
--p-exclude mitochondria,chloroplast,archaea,eukaryote,unassigned \
--o-filtered-table 03a_assign_taxonomy/otu_16S_filtered.qza 

qiime taxa filter-table \
--i-table 03a_assign_taxonomy/otu_18S.qza \
--i-taxonomy 03a_assign_taxonomy/vsearch_classify_id_80_taxonomy18/classification.qza \
--p-exclude bacteria,archaea,unassigned \
--o-filtered-table 03a_assign_taxonomy/otu_18S_filtered.qza 

# Export
qiime tools export \
--input-path 03a_assign_taxonomy/otu_16S_filtered.qza \
--output-path 03a_assign_taxonomy/otu_16S_filtered_exported

qiime tools export \
--input-path 03a_assign_taxonomy/otu_18S_filtered.qza \
--output-path 03a_assign_taxonomy/otu_18S_filtered_exported

# Export raw tables for analysis

qiime tools export \
--input-path 03a_assign_taxonomy/otu_16S.qza \
--output-path 03a_assign_taxonomy/otu_16S_exported

biom convert -i 03a_assign_taxonomy/otu_16S_exported/feature-table.biom --to-tsv \
-o 03a_assign_taxonomy/otu_16S_exported/feature-table.txt

qiime tools export \
--input-path 03a_assign_taxonomy/otu_18S.qza \
--output-path 03a_assign_taxonomy/otu_18S_exported

biom convert -i 03a_assign_taxonomy/otu_18S_exported/feature-table.biom --to-tsv \
-o 03a_assign_taxonomy/otu_18S_exported/feature-table.txt


