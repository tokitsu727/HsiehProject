#!/bin/bash
i=16
	(( 16<i && i<26)) && continue
	for j in {1..5}
	do
		echo $i
		PREPREFIX="00${i}_C1_T${j}"
		PREFIX="/data3/EddieFiles/0001/filtered_vcfs/${PREPREFIX}"
		FILE_MUTECT="${PREFIX}.MuTect_All_Filtered.snpEFF.filtered.vcf"
		FILE_SNVS="${PREFIX}.strelka.passed.somatic.snvs.snpEFF.filtered.vcf"
		FILE_INDELS="${PREFIX}.strelka.passed.somatic.indels.snpEFF.filtered.vcf"
		python3 pysam_filter.py -i $FILE_MUTECT -m $FILE_SNVS -o "${PREPREFIX}_MuTect_snvs.vcf"
		python3 pysam_filter.py -i $FILE_MUTECT -m $FILE_INDELS -o "${PREPREFIX}_MuTect_indels.vcf"
	done
