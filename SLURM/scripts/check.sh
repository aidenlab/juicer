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
#
# Sanity check once pipeline is complete to be sure number of lines adds up, deduping
# proceeded appropriately, and that files were created

# top level directory, can also be set in options
topDir=$(pwd)
# unique name for jobs in this run
groupname=$(basename $topDir)

splitdir=${topDir}"/juicer_splits"
donesplitdir=$topDir"/juicer_done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/juicer_output"
tmpdir=${topDir}"/juicer_HIC_tmp"


if [ -z $ARG1 ]
then
# Start by checking the statistics to see if they add up 
    res1=`cat ${splitdir}/*norm*res*`
    check1=`cat ${splitdir}/*norm*res* | awk '{s2+=$2; s3+=$3; s4+=$4; s5+=$5; s6+=$6}END{if (s2 != s3+s4+s5+s6){print 0}else{print 1}}'`
    if [ $check1 -eq 0 ] || [ -z "$res1" ]
    then
        echo "***! Error! The statistics do not add up. Alignment likely failed to complete on one or more files. Run relaunch_prep.sh"
        exit 100
    fi
fi

# Check the sizes of merged_sort versus the dups/no dups files to be sure
# no reads were lost
total=1
total2=0
total=`ls -l ${outputdir}/merged_sort.txt | awk '{print $5}'`
total2=`ls -l ${outputdir}/merged_nodups.txt ${outputdir}/dups.txt ${outputdir}/opt_dups.txt | awk '{sum = sum + $5}END{print sum}'`

if [ -z $total ] || [ -z $total2 ] || [ $total -ne $total2 ]
then
    echo "***! Error! The sorted file and dups/no dups files do not add up, or were empty. Merge or dedupping likely failed, restart pipeline with -S merge or -S dedup"
    exit 100
fi

wctotal=`cat ${splitdir}/*_linecount.txt | awk '{sum+=$1}END{print sum/4}'`
check2=`cat ${splitdir}/*norm*res* | awk '{s2+=$2;}END{print s2}'`
check4=`cat ${splitdir}/*norm*res* | awk '{s4+=$4;}END{print s4}'`
 
if [ -z $ARG1 ] 
then
    if [ $wctotal -ne $check2 ]
    then
        echo "***! Error! The number of reads in the fastqs (${wctotal}) is not the same as the number of reads reported in the stats (${check2}), likely due to a failure during alignment"
        exit 100
    fi
else
   # Alternate aligner check
   unmapamb=`awk 'NR%4==1{split($1,a,"_");split(a[1],b,"@"); c[b[2]]++}END{print length(c)}' ${splitdir}/*.fq`
   calctotal=`expr $check4 + $unmapamb`
   if [ $wctotal -ne $calctotal ]
   then
        echo "***! Error! The number of reads in the fastqs (${wctotal}) is not the same as the number of alignable reads reported in the stats (${check4}) plus unmapped/ambiguous (${unmapamb}), likely due to a failure during alignment. Run relaunch_prep.sh"
        exit 100
   else
       echo "Total: $wctotal"
       echo "Unalignable: $unmapamb"
   fi
fi

if [ -f ${outputdir}/${groupname}_inter.hic ] && [ -s ${outputdir}/${groupname}_inter.hic ] && [ -f ${outputdir}/${groupname}_inter_30.hic ] && [ -s ${outputdir}/${groupname}_inter_30.hic ]
then
    echo "(-: Pipeline successfully completed (-:";
    echo "Run cleanup.sh to remove the splits directory";
else
    echo "***! Error! Either inter.hic or inter_30.hic were not created"
    exit 100
fi
