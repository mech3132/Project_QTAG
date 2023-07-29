#!bin/bash

mkdir -p 05_run_QTAG
mkdir -p 05_run_QTAG_seasons
#mkdir -p 05_run_QTAG_clust99

## Loop through all datasets

for f in 04_adjust_files_for_QTAG/split_studies/*; do
fname=$(echo $f | sed 's/.*\///g')
	if [ -d 05_run_QTAG/$fname ]; then
		echo "$fname already exists; will not overwrite"
	else
		echo ''
		echo "Processing study $f"
		python3 QTAG/QTAG.py \
		-t $f/*_otu.txt \
		-m $f/*_meta.txt \
		-M 'Salinity' \
		-o 05_run_QTAG/$fname \
		--version 'both' \
		--overwrite True
	fi
done

for f in 04_adjust_files_for_QTAG/split_studies_seasons/*; do
fname=$(echo $f | sed 's/.*\///g')
	if [ -d 05_run_QTAG_seasons/$fname ]; then
		echo "$fname already exists; will not overwrite"
	else
		echo ''
		echo "Processing study $f"
		python3 QTAG/QTAG.py \
		-t $f/*_otu.txt \
		-m $f/*_meta.txt \
		-M 'Salinity' \
		-o 05_run_QTAG_seasons/$fname \
		--version 'both' \
		--overwrite True
	fi
done





#for f in 04_adjust_files_for_QTAG/split_studies/*; do
#fname=$(echo $f | sed 's/.*\///g')
#	if [ -d 05_run_QTAG/$fname ]; then
#		echo "$fname already exists; will not overwrite"
#	else
#		echo ''
#		echo "Processing study $f"
#		python3 QTAG/QTAG.py \
#		-t $f/*_otu.txt \
#		-m $f/*_meta.txt \
#		-M 'Salinity' \
#		-o 05_run_QTAG/$fname \
#		--overwrite True
#	fi
#done

#for f in 04_adjust_files_for_QTAG/split_studies/*; do
#fname=$(echo $f | sed 's/.*\///g')
#	if [ -d 05_run_QTAG/${fname}_original ]; then
#		echo "$fname already exists; will not overwrite"
#	else
#		echo ''
#		echo "Processing study $f"
#		python3 QTAG/QTAG_19april2020.py \
#		-t $f/*_otu.txt \
#		-m $f/*_meta.txt \
#		-M 'Salinity' \
#		-o 05_run_QTAG/${fname}_original 
#	fi
#done


#for f in 04_adjust_files_for_QTAG/split_studies_clust99/*; do
#fname=$(echo $f | sed 's/.*\///g')
#fadj=$(echo $f | sed 's/_clust99//g')
#	if [ -d 05_run_QTAG_clust99/$fname ]; then
#		echo "$fname already exists; will not overwrite"
#	else
#		echo ''
#		echo "Processing study $f"
#		python3 QTAG/QTAG.py \
#		-t $f/*_otu.txt \
#		-m $fadj/*_meta.txt \
#		-M 'Salinity' \
#		-o 05_run_QTAG_clust99/$fname \
#		--overwrite True
#	fi
#done

##`


#python3 QTAG/QTAG.py \
#-t 04_adjust_files_for_QTAG/split_studies/chen16/*_otu.txt \
#-m 04_adjust_files_for_QTAG/split_studies/chen16/*_meta.txt \
#-M 'Salinity' \
#-o 05_run_QTAG/TESTING \
#--overwrite True
