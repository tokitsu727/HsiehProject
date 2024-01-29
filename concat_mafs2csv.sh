awk 'FNR == 1 && NR!=1 || FNR == 2 && NR!=2{next;}{print}' maf_files/*.maf > Concat_Maf2Csv.csv
