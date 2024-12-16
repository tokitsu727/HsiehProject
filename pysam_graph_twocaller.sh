#!/bin/bash


for i in {01..30}
do
	(( 16<i && i<26)) && continue
	FILES=""
	for FILE in /data4/zarkoFiles/vcf2maf/mutectSnvPass/autosome/maf_files/*_00${i}*.maf; do
		FILES="${FILES} ${FILE}"
	done
	FILENAME="${i}_mutect"
	python3 pysam_filter.py --bulk_graph $FILES -o $FILENAME --one_graph
done
echo "--------SPLIT----------------"
for i in {01..30}
do
	(( 16<i && i<26)) && continue
	FILES=""
	for FILE in /data4/zarkoFiles/vcf2maf/strelkaSnvPass/autosome/maf_files/*_00${i}*.maf; do
		FILES="${FILES} ${FILE}"
	done
	FILENAME="${i}_snvs"
	python3 pysam_filter.py --bulk_graph $FILES -o $FILENAME  --one_graph
done
