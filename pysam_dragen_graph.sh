#!/bin/bash

for i in {01..30}
do
	#(( 16<i && i<26)) && continue
	FILES=""
	#for FILE in /data4/zarkoFiles/vcf2maf/dragenSnvPass/AF_dragen_AllChrom/maf_files/${i}*.maf; do
	for FILE in /home/go/Documents/TimothyOkitsu/dragen_vcf2maf/vcf2maf_sample/maf_files/${i}*.maf; do
		FILES="${FILES} ${FILE}"
	done
	FILENAME="${i}_dragen"
	python3 pysam_filter.py --bulk_graph $FILES -o $FILENAME --one_graph
done
