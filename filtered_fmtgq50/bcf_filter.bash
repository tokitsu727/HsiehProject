#!/bin/bash
#Filter input vcf

FILTER="MIN(FMT/GQ) >= 50 && MIN(FMT/DP) > 10" #&& MQ > 50 && MQRankSum > -2.5 && MQRankSum < 2.5 && SOR < 3 && FS < 60"

bcftools filter -O z -o $2 -i "$FILTER" $1

