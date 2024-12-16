#!/bin/bash

for i in {01..30}
do
	(( 16<i && i<26)) && continue
	FILES=""
	for FILE in /data4/TimothyOkitsu/freebayes_exp/vcf2maf_sample/maf_files/*_00${i}_*.maf; do
		FILES="${FILES} ${FILE}"
	done
	FILENAME="${i}freebayes"
	python3 pysam_filter.py --bulk_graph $FILES -o $FILENAME --freebayes_exp
done
