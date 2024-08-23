#!/bin/bash
# Juicer2 default directory is /aidenlab
JUICER2DIR=/aidenlab

# Create the directory structure that is expected by Juicer2 by default
mkdir -p "$JUICER2DIR/references" 
mkdir -p "$JUICER2DIR/restriction_sites"
# This folder contains the genome fasta file
mkdir -p $JUICER2DIR/demo/HIC003/fastq 

cd "$JUICER2DIR" || exit 
# Download the genome fasta files, reference genome and restriction sites
aria2c -x 16 -s 16 -i "$JUICER2DIR/download-demo.txt"
# Switch to the demo directory. This is where the fastq files are located. Juicer2 expects the fastq files to be in a folder called "fastq"
cd "$JUICER2DIR/demo/HIC003/"
# Run the Juicer2 pipeline with default parameters
juicer.sh -g hg19 2>&1 | tee output.log

exit 0

# To run the pipeline with custom parameters using only the genome .fasta and .fastq files, you can use the following commands

# Index a genome fasta file
bwa index /aidenlab/references/Homo_sapiens_assembly19.fasta

# create a restriction site file for MboI
/opt/juicer/misc/generate_site_positions.py MboI hg19 /aidenlab/references/Homo_sapiens_assembly19.fasta

# create the chrom.sizes file using the restriction site
awk 'BEGIN{OFS="\t"}{print $1, $NF}' hg19_MboI.txt > hg19.chrom.sizes

# Run the juicer Pipeline with specific parameters 
juicer.sh -d "$JUICER2DIR/demo/HIC003" -z "$JUICER2DIR/references/Homo_sapiens_assembly19.fasta" -p "hg19.chrom.sizes" -y "hg19_MboI.txt" -s MboI -t 8 -D "$JUICER2DIR" -a "demo_HIC003" 2>&1 | tee output.log
 

