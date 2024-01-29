#!/bin/bash

for i in {27..28}
do
	(( 16<i && i<26)) && continue
	for FILE in /media/go/wgsAnalysis/KgpOut_Echo/KUNUSCCLH_00${i}/analysis/*.HC_All.vcf;
	do
		PREFIX=${FILE%%.*}
		NOPATH="${PREFIX##*/}"
		echo $NOPATH
		FILE_MUTECT="${PREFIX}_MuTect_All_Filtered.vcf"
		FILE_SNVS="${PREFIX}.strelka.passed.somatic.snvs.vcf"
		FILE_INDELS="${PREFIX}.strelka.passed.somatic.indels.vcf"
		python3 pysam_filter.py -i $FILE_SNVS -m $FILE_MUTECT -o "${NOPATH}.MuTect_Strelka-snvs.vcf" &
		python3 pysam_filter.py -i $FILE_INDELS -m $FILE_MUTECT -o "${NOPATH}.MuTect_Strelka-indels.vcf" &
	done
	wait
			
done


