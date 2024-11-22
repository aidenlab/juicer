#!/bin/bash -l
# Juicer2 default directory is ./juicer2_demo
JUICER_DEMO_DIR=$(pwd)/juicer2_demo

# Create the directory structure that is expected by Juicer2 by default
mkdir -p "$JUICER_DEMO_DIR/references" 
mkdir -p "$JUICER_DEMO_DIR/restriction_sites"
# This folder contains the genome fasta file
mkdir -p $JUICER_DEMO_DIR/HIC003/fastq 

cd "$JUICER_DEMO_DIR" || exit 
# Download the genome fasta files, reference genome and restriction sites
aria2c --check-integrity --auto-file-renaming=false --conditional-get=true -x 16 -s 16 -i "/aidenlab/download-demo.txt"
# Switch to the demo directory. This is where the fastq files are located. Juicer2 expects the fastq files to be in a folder called "fastq"
cd "$JUICER_DEMO_DIR/HIC003/"
# Run the Juicer2 pipeline with default parameters
#juicer.sh -g hg19 2>&1 | tee output.log

# create the chrom.sizes file using the restriction site
awk 'BEGIN{OFS="\t"}{print $1, $NF}' $JUICER_DEMO_DIR/restriction_sites/hg19_MboI.txt > $JUICER_DEMO_DIR/restriction_sites/hg19.chrom.sizes


juicer.sh -g hg19 -z "$JUICER_DEMO_DIR/references/Homo_sapiens_assembly19.fasta" \
    -p "$JUICER_DEMO_DIR/restriction_sites/hg19.chrom.sizes" -y "$JUICER_DEMO_DIR/restriction_sites/hg19_MboI.txt" \
    -s MboI --qc 2>&1 | tee output.log

exit $?

# To run the pipeline with custom parameters using only the genome .fasta and .fastq files, you can use the following commands

# Index a genome fasta file
bwa index $JUICER_DEMO_DIR/references/Homo_sapiens_assembly19.fasta

# create a restriction site file for MboI
/opt/juicer/misc/generate_site_positions.py MboI hg19 $JUICER_DEMO_DIR/references/Homo_sapiens_assembly19.fasta


# Run the juicer Pipeline with specific parameters 
juicer.sh -d "$JUICER_DEMO_DIR/HIC003" -z "$JUICER_DEMO_DIR/references/Homo_sapiens_assembly19.fasta" \
    -p "$JUICER_DEMO_DIR /restriction_sites/hg19.chrom.sizes" -y "$JUICER_DEMO_DIR /restriction_sites/hg19_MboI.txt" \
    -s MboI -t 8 -D "$JUICER_DEMO_DIR" -a "demo_HIC003" 2>&1 | tee output.log
 

