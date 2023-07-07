#!/bin/bash
#Filter input vcf

FILTER="MIN(FMT/GQ) >= 105 && MIN(FMT/DP) > 10"
IFILE=$1
OFILE_TMP="tmp.vcf"
OFILE=$2
OFILE_GZ="$OFILE_TMP".gz

#python3 haploComparator.py $IFILE $OFILE_TMP
bcftools filter -O z -o $OFILE_GZ -i "$FILTER" $IFILE #Change to OFILE_TMP if doing haploComparator step. But, should be unnecessary
#rm $OFILE_TMP
gunzip $OFILE_GZ
python3 AFfilter.py $OFILE_TMP $OFILE
rm $OFILE_TMP

