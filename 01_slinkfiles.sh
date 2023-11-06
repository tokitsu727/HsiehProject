#/bin/bash

cd renamed_vcfs
#FILES='../../Collaborate/KUNUSCCLH_00*/analysis/*vcf'
FILES='/media/go/wgsAnalysis/KUNUSCCLH_00*/analysis/*vcf'

for FILE in $FILES;
do
    FILE2=${FILE/[\.\/a-zA-Z]*\//}
    FILE2=${FILE2/KUNUSCCLH_/}
    FILE2=${FILE2/_01_[a-zA-Z]*_C/_C}
    FILE2=${FILE2/_KHWGSH[_-a-zA-Z]*_T/_T}
    FILE2=${FILE2/_KHWGSH_A[0-9][0-9][0-9][0-9][0-9][\.\_]/.}
    echo "$FILE"
    echo "=> $FILE2"
    echo ""
    ln -s $FILE $FILE2
done
