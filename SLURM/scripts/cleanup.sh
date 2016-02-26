#!/bin/bash
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########

# Script to clean up big repetitive files and zip fastqs. Run after you are 
# sure the pipeline ran successfully.  Run from top directory (HIC001 e.g.).

juiceDir="/home/keagen/programs/juicer/SLURM"
# top level directory, can also be set in options
topDir=$(pwd)
# unique name for jobs in this run
groupname=$(basename $topDir)

#output messages
outDir="$topDir/debug"

total=`ls -l juicer_aligned/merged_sort.txt | awk '{print $5}'`
total2=`ls -l juicer_aligned/merged_nodups.txt juicer_aligned/dups.txt juicer_aligned/opt_dups.txt | awk '{sum = sum + $5}END{print sum}'`
if [ $total -eq $total2 ]; then rm juicer_aligned/merged_sort.txt; rm -r juicer_splits; else echo "Problem: The sum of merged_nodups and the dups files is not the same size as merged_sort.txt"; fi

jid=`sbatch <<- CLEANUP | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p kornberg,owners,normal
	#SBATCH -o $outDir/cleanup-%j.out
	#SBATCH -e $outDir/cleanup-%j.err
	#SBATCH -t 16:00:00
	#SBATCH -n 1
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=250G
	#SBATCH -J "${groupname}_cleanup"

	module load pbzip2
	pbzip2 -p16 fastq/*.fastq
CLEANUP`