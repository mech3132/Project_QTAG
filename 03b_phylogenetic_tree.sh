#!bin/bash

conda activate qiime2-2022.8

mkdir -p 03b_phylogenetic_tree

# Get SILVA tree for insertion
#wget \
#  -O "03_phylogenetic_tree/sepp-refs-silva-128.qza" \
#  "https://data.qiime2.org/2021.11/common/sepp-refs-silva-128.qza"

### Fragment insertion SEPP

## Conducted on compute canada cedar, qiime 2022.8

#mkdir -p 03_phylogenetic_tree/sepp
#qiime fragment-insertion sepp \
#  --i-representative-sequences 03a_assign_taxonomy/repset_16S.qza \
#  --i-reference-database 03b_phylogenetic_tree/sepp-refs-silva-128.qza \
#  --o-tree 03b_phylogenetic_tree/sepp/16S_2022-08-03/insertion-tree_16S.qza \
#  --o-placements 03b_phylogenetic_tree/sepp/16S_2022-08-03/insertion-placements_16S.qza
#  
#qiime fragment-insertion sepp \
#  --i-representative-sequences 03a_assign_taxonomy/repset_18S.qza \
#  --i-reference-database 03_phylogenetic_tree/sepp-refs-silva-128.qza \
#  --o-tree 03_phylogenetic_tree/sepp/insertion-tree_18S.qza \
#  --o-placements 03_phylogenetic_tree/sepp/insertion-placements_18S.qza
#  
#  
# # Filter tree by 16S and 18S
######16S
qiime fragment-insertion filter-features \
--i-tree 03b_phylogenetic_tree/insertion-tree_16S_qtag.qza \
--i-table 03a_assign_taxonomy/otu_16S_filtered.qza \
--output-dir 03b_phylogenetic_tree/filtered_table_16S

qiime phylogeny filter-tree \
--i-tree 03b_phylogenetic_tree/insertion-tree_16S_qtag.qza \
--i-table 03b_phylogenetic_tree/filtered_table_16S/filtered_table.qza \
--output-dir 03b_phylogenetic_tree/filtered-tree_16S

qiime tools export \
--input-path 03b_phylogenetic_tree/filtered-tree_16S/filtered_tree.qza \
--output-path 03b_phylogenetic_tree/filtered-tree_16S/exported_16S_tree

qiime tools export \
--input-path 03b_phylogenetic_tree/filtered_table_16S/filtered_table.qza \
--output-path 03b_phylogenetic_tree/exported_16S_table

biom convert --to-tsv \
-i 03b_phylogenetic_tree/exported_16S_table/feature-table.biom \
-o 03b_phylogenetic_tree/exported_16S_table/feature-table.txt


######18S

qiime fragment-insertion filter-features \
--i-tree 03b_phylogenetic_tree/insertion-tree_18S_qtag.qza \
--i-table 03a_assign_taxonomy/otu_18S_filtered.qza \
--output-dir 03b_phylogenetic_tree/filtered_table_18S

qiime phylogeny filter-tree \
--i-tree 03b_phylogenetic_tree/insertion-tree_18S_qtag.qza \
--i-table 03b_phylogenetic_tree/filtered_table_18S/filtered_table.qza \
--output-dir 03b_phylogenetic_tree/filtered-tree_18S

qiime tools export \
--input-path 03b_phylogenetic_tree/filtered-tree_18S/filtered_tree.qza \
--output-path 03b_phylogenetic_tree/filtered-tree_18S/exported_18S_tree

qiime tools export \
--input-path 03b_phylogenetic_tree/filtered_table_18S/filtered_table.qza \
--output-path 03b_phylogenetic_tree/exported_18S_table

biom convert --to-tsv \
-i 03b_phylogenetic_tree/exported_18S_table/feature-table.biom \
-o 03b_phylogenetic_tree/exported_18S_table/feature-table.txt

######### For troubleshooting

#qiime tools export --input-path 03_phylogenetic_tree/repset_16S.qza --output-path 03_phylogenetic_tree/repset_16S

#qiime tools export --input-path 03_phylogenetic_tree/repset_18S.qza --output-path 03_phylogenetic_tree/repset_18S



