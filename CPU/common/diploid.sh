#!/bin/bash

stem=/aidenlab/work/neva/patski_

for i in chrX chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrY
  do
  awk -v chr=$i '$2==chr && $6==chr && $9 >= 10 && $12 >= 10' merged_nodups.txt | ${juiceDir}/scripts/diploid.pl -s ${stem}${i}_chr_pos.txt -o ${stem}${i}_paternal_maternal.txt > diploid_${i}.txt
done
cat diploid_chr1.txt diploid_chr10.txt diploid_chr11.txt diploid_chr12.txt diploid_chr13.txt diploid_chr14.txt diploid_chr15.txt diploid_chr16.txt diploid_chr17.txt diploid_chr18.txt diploid_chr19.txt diploid_chr2.txt diploid_chr3.txt diploid_chr4.txt diploid_chr5.txt diploid_chr6.txt diploid_chr7.txt diploid_chr8.txt diploid_chr9.txt diploid_chrX.txt diploid_chrY.txt > diploid.txt

awk -f ${juiceDir}/scripts/diploid_split.awk diploid.txt 
sort -k2,2d -m maternal.txt maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > mat.txt 
${juiceDir}/scripts/juicer_tools pre mat.txt maternal.hic mm10 
sort -k2,2d -m paternal.txt paternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > pat.txt 
/aidenlab/juicebox pre pat.txt paternal.hic mm10

