#!/bin/bash


for i in {01..30}
do
	(( 16<i && i<26)) && continue
	FILES=""
	for FILE in /data4/zarkoFiles/vcf2maf/dragenSnvPass/AF_dragen_AllChrom/maf_files/${i}*.maf; do
		FILES="${FILES} ${FILE}"
	done
	FILENAME="${i}_dragen"
	python3 pysam_filter.py --bulk_graph $FILES -o $FILENAME --reflect_graph
done
