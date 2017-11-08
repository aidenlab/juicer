#!/bin/bash
# Converts bam file to merged_nodups.txt file
# Calls samtools to filter duplicates then 
# chimeric_blacklist to associate read pair group
# Then assigns fragment and sorts
# Sorting should be optional
# Usage: bamtotxt.sh <bam_file> <site_file> <output_file> [sort]
set -e

if [ $# -lt 3 ] || [ $# -gt 4 ]
then
    echo "Usage: $0 <bam_file> <site_file> <output_file> [sort]"
    exit 1
else
    bamfile=$1
    sitefile=$2 
    outputfile=$3
fi
if [ $# -eq 4 ]
then
    sort=1
fi

samtools view -F 1024 "$bamfile" | awk -v fname=$bamfile -f "$( dirname "${BASH_SOURCE[0]}" )"/chimeric_blacklist.awk
"$( dirname "${BASH_SOURCE[0]}" )"/fragment.pl "${bamfile}_norm.txt" "${bamfile}.frag.txt" "$sitefile"
if [ -n $sort ]
then
    sort -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n "${bamfile}.frag.txt" > "$outputfile"
else
    mv "${bamfile}.frag.txt" "$outputfile"
fi
rm "${bamfile}_norm.txt" "${bamfile}.frag.txt"
