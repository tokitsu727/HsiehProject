#!/bin/bash

for i in {01..30}
do
	(( 16<i && i<26)) && continue
	for FILE in /media/go/wgsAnalysis/KUNUSCCLH_00${i}/analysis/*.HC_All.vcf;
	do
		PREFIX=${FILE%%.*}
		NOPATH="${FILE##*/}"
		echo $NOPATH
		FILE_MUTECT="${PREFIX}_MuTect_All_Filtered.vcf"
		FILE_SNVS="${PREFIX}.strelka.passed.somatic.snvs.vcf"
		FILE_INDELS="${PREFIX}.strelka.passed.somatic.indels.vcf"
		python3 pysam_filter.py -i $FILE_MUTECT -m $FILE_SNVS -o "${NOPATH}_MuTect_Strelkasnvs.vcf"
		python3 pysam_filter.py -i $FILE_MUTECT -m $FILE_INDELS -o "${NOPATH}_MuTect_Strelkaindels.vcf"
	done
			
done


