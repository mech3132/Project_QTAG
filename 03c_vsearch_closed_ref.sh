#!bin/bash

conda activate qiime2-2021.11

mkdir -p 03_vsearch_closed_ref
mkdir -p 03_vsearch_closed_ref/imported
mkdir -p 03_vsearch_closed_ref/filtered_repsets
mkdir -p 03_vsearch_closed_ref/exported_tables
mkdir -p 03_vsearch_closed_ref/exported_sequences


# Download ref

if [[ -f 03_vsearch_closed_ref/silva-138-99-seqs.qza ]]; then
 echo "silva seqs already downloaded"
else
wget \
  -O "03_vsearch_closed_ref/silva-138-99-seqs.qza" \
  "https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza"
fi
if [[ -f 03_vsearch_closed_ref/silva-138-99-tax.qza ]]; then
 echo "silva seqs already downloaded"
else
wget \
  -O "03_vsearch_closed_ref/silva-138-99-tax.qza" \
  "https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza"
fi


# Import repset
qiime tools import \
--input-path 02_extract_sequences/repset.fasta \
--output-path 03_vsearch_closed_ref/repset.qza \
--type 'FeatureData[Sequence]'


for f in 02_adjust_files_for_processing/studies_split/*; do
	fname=$(echo $f | sed 's/^.*\///g' )
	# Import table
	biom convert -i $f/otu.txt -o 03_vsearch_closed_ref/imported/$fname.biom \
	--to-hdf5 --table-type 'OTU table'

	qiime tools import \
	--input-path 03_vsearch_closed_ref/imported/$fname.biom \
	--output-path 03_vsearch_closed_ref/imported/$fname.qza \
	--type 'FeatureTable[Frequency]' \
	--input-format BIOMV210Format 

	qiime feature-table filter-seqs \
	--i-data 03_vsearch_closed_ref/repset.qza \
	--i-table 03_vsearch_closed_ref/imported/$fname.qza \
	--o-filtered-data 03_vsearch_closed_ref/filtered_repsets/$fname-repset.qza
	
done

# Run through closed ref clustering


## Do vsearch closed ref clustering with both strands for wang and herleman; I think strands are backward
#start=`date +%s`
#qiime vsearch cluster-features-open-reference \
#--i-sequences 03_vsearch_closed_ref/filtered_repsets/wang2021-repset.qza \
#--i-table 03_vsearch_closed_ref/imported/wang2021.qza \
#--i-reference-sequences 03_vsearch_closed_ref/silva-138-99-seqs.qza \
#--p-perc-identity 0.99 \
#--p-strand 'both' \
#--output-dir 03_vsearch_closed_ref/vsearch99_wang2021

#qiime vsearch cluster-features-open-reference \
#--i-sequences 03_vsearch_closed_ref/filtered_repsets/herleman2016-repset.qza \
#--i-table 03_vsearch_closed_ref/imported/herleman2016.qza \
#--i-reference-sequences 03_vsearch_closed_ref/silva-138-99-seqs.qza \
#--p-perc-identity 0.99 \
#--p-strand 'both' \
#--output-dir 03_vsearch_closed_ref/vsearch99_herleman2016
#end=`date +%s`

#runtime=$((end-start))

#echo $runtime


for f in 02_adjust_files_for_processing/studies_split/*; do
fname=$(echo $f | sed 's/^.*\///g' )

if [[ -d 03_vsearch_closed_ref/vsearch99_${fname} ]]; then
echo "$fname already exists; not re-clustering"
else
echo "DOING $fname"
qiime vsearch cluster-features-open-reference \
--i-sequences 03_vsearch_closed_ref/filtered_repsets/$fname-repset.qza \
--i-table 03_vsearch_closed_ref/imported/$fname.qza \
--i-reference-sequences 03_vsearch_closed_ref/silva-138-99-seqs.qza \
--p-perc-identity 0.99 \
--p-strand 'both' \
--output-dir 03_vsearch_closed_ref/vsearch99_${fname}
fi
qiime tools export \
--input-path 03_vsearch_closed_ref/vsearch99_${fname}/clustered_table.qza \
--output-path 03_vsearch_closed_ref/vsearch99_${fname}/clustered_table_exported

biom convert -i 03_vsearch_closed_ref/vsearch99_${fname}/clustered_table_exported/feature-table.biom \
-o 03_vsearch_closed_ref/exported_tables/${fname}_clustered_table.txt \
--to-tsv

done



####### Export
# Export silva database
#qiime tools export \
#--input-path 03_vsearch_closed_ref/silva-138-99-tax.qza \
#--output-path 03_vsearch_closed_ref/silva-138-99-tax

# Merge feature sequences for taxonomic assignment


for f in 02_adjust_files_for_processing/studies_split/*; do
fname=$(echo $f | sed 's/^.*\///g' )
cp 03_vsearch_closed_ref/vsearch99_${fname}/clustered_sequences.qza 03_vsearch_closed_ref/exported_sequences/${fname}_repset.qza
done

qiime feature-table merge-seqs \
--i-data 03_vsearch_closed_ref/exported_sequences/* \
--o-merged-data 03_vsearch_closed_ref/repset_all_clustered99.qza

### RUN THIS ON COMPUTECANADA
#qiime feature-classifier classify-consensus-vsearch \
#--i-query 03_vsearch_closed_ref/repset_all_clustered99.qza \
#--i-reference-reads 03_vsearch_closed_ref/silva-138-99-seqs.qza \
#--i-reference-taxonomy 03_vsearch_closed_ref/silva-138-99-tax.qza \
#--p-perc-identity 0.8 \
#--p-strand 'both' \
#--output-dir 03_vsearch_closed_ref/vsearch_classify_cluster99_id80_taxonomy


qiime tools export \
--input-path 03_vsearch_closed_ref/repset_all_clustered99.qza \
--output-path  03_vsearch_closed_ref/repset_all_clustered99_exported


qiime tools export \
--input-path 03_vsearch_closed_ref/2022-08-30_vsearch_clust99/classification.qza \
--output-path  03_vsearch_closed_ref/2022-08-30_vsearch_clust99/classification_exported


