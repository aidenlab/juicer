#!/bin/bash

inputFile=$1
chunkSize=$2

awk -v chunkSize=$chunkSize 'BEGIN{bytecounter=0; numLines=0; chunkCounter=0;}{
	if (bytecounter==0) {
		chr1=$2;
		chr2=$6;
		gsub(/ /, "", chr1);
		gsub(/ /, "", chr2);
		#chr1=tolower(chr1);
		#chr2=tolower(chr2);
		#gsub(/chr/, "", chr1);
		#gsub(/chr/, "", chr2);	
		#chr1=toupper(chr1);
		#chr2=toupper(chr2);
		#print chr1"-"chr2 "," bytecounter;
		startPos = bytecounter;
		numLines++;
		currentChr1 = chr1;
		currentChr2 = chr2;
	} 
	chr1=$2;
        chr2=$6;
        gsub(/ /, "", chr1);
        gsub(/ /, "", chr2);
        #chr1=tolower(chr1);
        #chr2=tolower(chr2);
        #gsub(/chr/, "", chr1);
        #gsub(/chr/, "", chr2);
        #chr1=toupper(chr1);
        #chr2=toupper(chr2);
	if (chr1!=currentChr1 || chr2!=currentChr2) {
		print currentChr1"-"currentChr2 "," chunkCounter "," startPos "," bytecounter-startPos;
		currentChr1 = chr1;
		currentChr2 = chr2;
		startPos = bytecounter;
		chunkCounter=0;
		numLines=1;
		bytecounter += length($0) + 1;
	} else if (numLines==chunkSize) {
		print currentChr1"-"currentChr2 "," chunkCounter "," startPos "," bytecounter-startPos;
		startPos = bytecounter;
		chunkCounter++;
		numLines=1;
		bytecounter += length($0) + 1;
	} else {
		numLines++
		bytecounter += length($0) + 1;
	}
	  
}END {
	print currentChr1"-"currentChr2 "," chunkCounter "," startPos "," bytecounter-startPos;
}' "$inputFile"
