#!/bin/bash

for i in {01..30}
do
	(( 16<i && i<26)) && continue
	for FILE in /media/go/wgsAnalysis/KgpOut_Echo/KUNUSCCLH_00${i}/analysis/*.HC_All.vcf;
	do
		PREFIX=${FILE%%.*}
		NOPATH="${PREFIX##*/}"
		[[ $NOPATH =~ (.*)(\_T[0-9])(.*) ]]
		TUBE=${BASH_REMATCH[2]}
		echo $NOPATH
		TUBE=${TUBE#_}
		DRAGEN=/data3/ZarkoFiles/VcfConverter/crypt_test/dragen_hardfilter/vcf_files/${i}C1_${i}${TUBE}.hard-filtered.vcf
		echo $DRAGEN
		FILE_MUTECT="${PREFIX}_MuTect_All_Filtered.vcf"
		FILE_SNVS="${PREFIX}.strelka.passed.somatic.snvs.vcf"
		FILE_INDELS="${PREFIX}.strelka.passed.somatic.indels.vcf"
		python3 pysam_filter.py -i $DRAGEN -m $FILE_SNVS -o "${NOPATH}.Dragen_StrelkaSnvs.vcf" -d "merged_files_Dragen_strelka" &
		python3 pysam_filter.py -i $DRAGEN -m $FILE_SNVS -o "${NOPATH}.Dragen_StrelkaIndels.vcf" -d "merged_files_Dragen_strelka" &
	done
	wait
			
done


