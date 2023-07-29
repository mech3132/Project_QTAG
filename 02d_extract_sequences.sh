#!bin/bash

mkdir -p 02d_extract_sequences

awk ' NR>1 {print ">"$1"\n"$2}' 02b_adjust_files_for_processing/asvList.txt > 02d_extract_sequences/repset.fasta

for f in 02b_adjust_files_for_processing/studies_split/*/asvList.txt; do
	temppw=$(echo $f | sed 's/asvList.txt/repset.fasta/g')
	awk ' NR>1 {print ">"$1"\n"$2}' $f > $temppw
done
