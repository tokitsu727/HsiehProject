#!/bin/bash


for i in {01..30}
do
	(( 16<i && i<26)) && continue
	FILES=""
	for FILE in maf_files/*_00${i}*-indels.maf; do
		FILES="${FILES} ${FILE}"
	done
	FILENAME="${i}_indels"
	python3 pysam_filter.py --bulk_graph $FILES -o $FILENAME --reflect_graph
done

for i in {01..30}
do
	(( 16<i && i<26)) && continue
	FILES=""
	for FILE in maf_files/*_00${i}*-snvs.maf; do
		FILES="${FILES} ${FILE}"
	done
	FILENAME="${i}_snvs"
	python3 pysam_filter.py --bulk_graph $FILES -o $FILENAME --reflect_graph
done
