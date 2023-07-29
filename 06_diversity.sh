#!bin/bash

conda activate qiime2-2022.8

mkdir -p 06_diversity


# Coverage based rarefaction

#python3 ../../Code/coverage_based_rarefaction.py \
#-i 04_adjust_files_for_QTAG/split_studies/campbell2013/campbell2013_otu.txt \
#-m 100 \
#-M 300 \
#-b 50 \
#-r 3 \
#-c 0.95 \
#-o 06_diversity/coverage_based_rarefaction

#python3 ../../Code/coverage_based_rarefaction.py \
#-i 04_adjust_files_for_QTAG/otu16.txt \
#-m 100 \
#-M 5000 \
#-b 100 \
#-r 10 \
#-c 0.85 \
#-o 06_diversity/coverage_based_rarefaction_16


python3 ../../Code/coverage_based_rarefaction.py \
-i 04_adjust_files_for_QTAG/otu.txt \
-m 100 \
-M 5000 \
-b 100 \
-r 10 \
-c 0.85 \
-o 06_diversity/coverage_based_rarefaction



# REMOVE NAME
#sed 's/^X/#OTU ID/g' 06_diversity/coverage_based_rarefaction_16/rarefiedTable.txt >  06_diversity/coverage_based_rarefaction_16/rarefiedTable_edit.txt

sed 's/^X/#OTU ID/g' 06_diversity/coverage_based_rarefaction/rarefiedTable.txt >  06_diversity/coverage_based_rarefaction/rarefiedTable_edit.txt


### Diversity matrices


biom convert -i 06_diversity/coverage_based_rarefaction/rarefiedTable_edit.txt \
-o 06_diversity/otu.biom \
--table-type='OTU table' \
--to-hdf5 

qiime tools import \
--input-path 06_diversity/otu.biom \
--output-path 06_diversity/otu.qza \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format 

# separate based on 16S/18S
qiime feature-table filter-samples \
--i-table 06_diversity/otu.qza \
--m-metadata-file 04_adjust_files_for_QTAG/all_meta.txt \
--p-where "target_region=='16S'" \
--o-filtered-table 06_diversity/otu16.qza

qiime feature-table filter-samples \
--i-table 06_diversity/otu.qza \
--m-metadata-file 04_adjust_files_for_QTAG/all_meta.txt \
--p-where "target_region=='18S'" \
--o-filtered-table 06_diversity/otu18.qza



#biom convert -i 04_adjust_files_for_QTAG/otu18.txt \
#--to-hdf5 -o 06_diversity/otu18.biom \
#--table-type 'OTU table'

#qiime tools import \
#--input-path 06_diversity/otu18.biom \
#--output-path 06_diversity/otu18.qza \
#--type 'FeatureTable[Frequency]' \
#--input-format BIOMV210Format 


###### Visualize #####
qiime feature-table summarize --i-table 06_diversity/otu16.qza --o-visualization 06_diversity/summary_16S

qiime feature-table summarize --i-table 06_diversity/otu18.qza --o-visualization 06_diversity/summary_18S

####### DIVERSITY (NONRAREFIED) #######
## WU
qiime diversity beta-phylogenetic \
--i-table 06_diversity/otu16.qza \
--i-phylogeny 03b_phylogenetic_tree/filtered-tree_16S/filtered_tree.qza \
--p-metric 'weighted_unifrac' \
--o-distance-matrix 06_diversity/dm_WU_16S

qiime diversity beta-phylogenetic \
--i-table 06_diversity/otu18.qza \
--i-phylogeny 03b_phylogenetic_tree/filtered-tree_18S/filtered_tree.qza \
--p-metric 'weighted_unifrac' \
--o-distance-matrix 06_diversity/dm_WU_18S

qiime tools export \
--input-path 06_diversity/dm_WU_16S.qza \
--output-path 06_diversity/dm_WU_16S

qiime tools export \
--input-path 06_diversity/dm_WU_18S.qza \
--output-path 06_diversity/dm_WU_18S

## UWU
qiime diversity beta-phylogenetic \
--i-table 06_diversity/otu16.qza \
--i-phylogeny 03b_phylogenetic_tree/filtered-tree_16S/filtered_tree.qza \
--p-metric 'unweighted_unifrac' \
--o-distance-matrix 06_diversity/dm_UWU_16S

qiime diversity beta-phylogenetic \
--i-table 06_diversity/otu18.qza \
--i-phylogeny 03b_phylogenetic_tree/filtered-tree_18S/filtered_tree.qza \
--p-metric 'unweighted_unifrac' \
--o-distance-matrix 06_diversity/dm_UWU_18S

qiime tools export \
--input-path 06_diversity/dm_UWU_16S.qza \
--output-path 06_diversity/dm_UWU_16S

qiime tools export \
--input-path 06_diversity/dm_UWU_18S.qza \
--output-path 06_diversity/dm_UWU_18S



####### DIVERSITY (RAREFIED) #######
## WU
#qiime diversity beta-rarefaction \
#--i-table 06_diversity/otu16.qza \
#--i-phylogeny 03_phylogenetic_tree/2022-08-30_phylogenetictree/insertion-tree_16S.qza \
#--p-metric 'weighted_unifrac' \
#--p-clustering-method 'upgma' \
#--m-metadata-file 04_adjust_files_for_QTAG/all_meta_MANUAL.txt \
#--p-sampling-depth 2000 \
#--output-dir 06_diversity/WU_16_2000

#qiime diversity beta-rarefaction \
#--i-table 06_diversity/otu16.qza \
#--i-phylogeny 03_phylogenetic_tree/2022-08-30_phylogenetictree/insertion-tree_16S.qza \
#--p-metric 'weighted_unifrac' \
#--p-clustering-method 'upgma' \
#--m-metadata-file 04_adjust_files_for_QTAG/all_meta_MANUAL.txt \
#--p-sampling-depth 100 \
#--output-dir 06_diversity/WU_16_100

#qiime diversity beta-phylogenetic \
#--i-table 06_diversity/otu18.qza \
#--i-phylogeny 03_phylogenetic_tree/2022-08-30_phylogenetictree/insertion-tree_18S.qza \
#--p-metric 'weighted_unifrac' \
#--p-clustering-method 'upgma' \
#--p-sampling-depth 2000 \
#--output-dir WU_18_2000

#qiime tools export \
#--input-path 06_diversity/dm_WU_16S.qza \
#--output-path 06_diversity/dm_WU_16S

#qiime tools export \
#--input-path 06_diversity/dm_WU_18S.qza \
#--output-path 06_diversity/dm_WU_18S

### UWU
#qiime diversity beta-rarefaction \
#--i-table 06_diversity/otu16.qza \
#--i-phylogeny 03_phylogenetic_tree/2022-08-30_phylogenetictree/insertion-tree_16S.qza \
#--p-metric 'unweighted_unifrac' \
#--p-clustering-method 'upgma' \
#--p-sampling-depth 2000 \
#--m-metadata-file 04_adjust_files_for_QTAG/all_meta_MANUAL.txt \
#--m-metadata-file \
#--output-dir 06_diversity/UWU_16_2000

#qiime diversity beta-rarefaction \
#--i-table 06_diversity/otu16.qza \
#--i-phylogeny 03_phylogenetic_tree/2022-08-30_phylogenetictree/insertion-tree_16S.qza \
#--p-metric 'unweighted_unifrac' \
#--p-clustering-method 'upgma' \
#--p-sampling-depth 2000 \
#--m-metadata-file 04_adjust_files_for_QTAG/all_meta_MANUAL.txt \
#--m-metadata-file \
#--output-dir 06_diversity/UWU_16_2000


#qiime diversity beta-phylogenetic \
#--i-table 06_diversity/otu18.qza \
#--i-phylogeny 03_phylogenetic_tree/2022-08-30_phylogenetictree/insertion-tree_18S.qza \
#--p-metric 'unweighted_unifrac' \
#--o-distance-matrix 06_diversity/dm_UWU_18S

#qiime tools export \
#--input-path 06_diversity/dm_UWU_16S.qza \
#--output-path 06_diversity/dm_UWU_16S

#qiime tools export \
#--input-path 06_diversity/dm_UWU_18S.qza \
#--output-path 06_diversity/dm_UWU_18S




