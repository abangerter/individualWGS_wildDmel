#! /bin/bash


inputDir="/scratch/ab5dr/wildDmel2016/cnvOut/v4outManta"
outputDir="/scratch/ab5dr/wildDmel2016/cnvOut/processMantaV4"

samples=($( ls ${inputDir}/ ))
for ID in "${samples[@]}"; do
	zgrep -e "^#" -e "INS" -e "DEL" -e "DUP" ${inputDir}/${ID}/diploidSV.vcf.gz | grep -e "^#" -e "^X" -e "^2L" -e "^2R" -e "^3R" -e "^3L" > ${outputDir}/${ID}_SVscored_subset.vcf
done

for ID in "${samples[@]}"; do
	grep -v "^#" ${outputDir}/${ID}_SVscored_subset.vcf | awk -v id=${ID} '{OFS="\t";}{IFS="\t";}{print id,$1,$2,$3,$4,$5,$6,$7,$8}' >> ${outputDir}/all_SVscored_subset.txt
done
