#!/bin/bash

inputFile=$1

awk '{chr1=$2;chr2=$6;gsub(/ /, "", chr1);gsub(/ /, "", chr2);chr1=tolower(chr1);chr2=tolower(chr2);gsub(/chr/, "", chr1);gsub(/chr/, "", chr2);chr1=toupper(chr1);chr2=toupper(chr2); print $0 >> "merged_nodups.txt_"chr1"_"chr2;}' "$inputFile"
