#!/bin/bash

inputFile=$1

awk 'BEGIN{bytecounter=0}{
	if (bytecounter==0) {
		chr1=$2;
		chr2=$6;
		gsub(/ /, "", chr1);
		gsub(/ /, "", chr2);
		chr1=tolower(chr1);
		chr2=tolower(chr2);
		gsub(/chr/, "", chr1);
		gsub(/chr/, "", chr2);	
		chr1=toupper(chr1);
		chr2=toupper(chr2);
		print chr1"_"chr2 "," bytecounter;
		currentChr1 = chr1;
		currentChr2 = chr2;
	} 
	chr1=$2;
        chr2=$6;
        gsub(/ /, "", chr1);
        gsub(/ /, "", chr2);
        chr1=tolower(chr1);
        chr2=tolower(chr2);
        gsub(/chr/, "", chr1);
        gsub(/chr/, "", chr2);
        chr1=toupper(chr1);
        chr2=toupper(chr2);
	if (chr1!=currentChr1 || chr2!=currentChr2) {
		print chr1"_"chr2 "," bytecounter;
		currentChr1 = chr1;
		currentChr2 = chr2;
		bytecounter += length($0) + 1;
	} else {
		bytecounter += length($0) + 1;
	}
	  
}' "$inputFile"
