#!/bin/bash

for FILE in /media/go/wgsAnalysis/Analyses/KGP/vcf2maf/Mutect_AllChrom/maf_files/*;
do
	PREFIX=${FILE%%.*}
	NOPATH="${PREFIX##*/}"
	DRAGEN_PREFIX="/media/go/wgsAnalysis/Analyses/Dragen/vcf2maf/AF_dragen_AllChrom_SNVs/maf_files/"
	echo $NOPATH
	if [[ $NOPATH =~ ^KUNUSCCLH_00(.+?)_(.+?)_T(.+?)_(.+?)$ ]]
	then
		FILE_DRAGEN="${BASH_REMATCH[1]:0:2}T${BASH_REMATCH[3]:0:1}.dragen.vcf.PASS.vcf.snv.vcf.AllChrom.maf"
		FILE_DRAGEN="${DRAGEN_PREFIX}${FILE_DRAGEN}"
		FILE_OUT="DragenMutectCsv/${BASH_REMATCH[1]:0:2}T${BASH_REMATCH[3]:0:1}.DragenMutectCombined.csv"
		echo $FILE_OUT
		echo $FILE_DRAGEN
		python3 pysam_filter.py -i $FILE_DRAGEN -m $FILE -o $FILE_OUT -c 
	fi
done
			


