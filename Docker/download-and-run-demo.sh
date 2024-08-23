#!/bin/bash
JUICER2DIR=/aidenlab

mkdir -p "$JUICER2DIR/references" 
mkdir -p "$JUICER2DIR/restriction_sites"
mkdir -p $JUICER2DIR/demo/HIC003/fastq 

cd "$JUICER2DIR" || exit 
aria2c -x 16 -s 16 -i "$JUICER2DIR/download-demo.txt"
cd "$JUICER2DIR/demo/HIC003/"
juicer.sh -g hg19 2>&1 | tee output.log
#/opt/juicer/misc/generate_site_positions.py MboI hg19 /aidenlab/references/Homo_sapiens_assembly19.fasta
#awk 'BEGIN{OFS="\t"}{print $1, $NF}' hg19_MboI.txt > hg19.chrom.sizes
#juicer.sh -d "$JUICER2DIR/demo/HIC003" -z "$JUICER2DIR/references/Homo_sapiens_assembly19.fasta" -p "hg19.chrom.sizes" -y "hg19_MboI.txt" -s MboI -t 8 -D "$JUICER2DIR" -a "demo_HIC003" 2>&1 | tee output.log
 

