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
# MegaMap script. 
#
#                                                                       
# [topDir]		    - Should contain the results of all base experiments
#
# From the top-level directory, the following two directory is created:
#                                                                              
# [topDir]/mega     - Location of result of processing the mega map
#
# Juicer version 1.5
juicer_version="1.5" 
## Set the following variables to work with your system

## use cluster load commands:
#usePath=""
#load_java=""
#load_cluster=""
#load_coreutils=""
#load_gpu=""
# Juicer directory, contains scripts/ and restriction_sites/
juiceDir="/opt/juicer"
# default queue, can also be set in options
queue="Stat"
# default long queue, can also be set in options
long_queue="Merge"
# unique name for jobs in this run
groupname="a$(date +%s)"

## Default options, overridden by command line arguments

# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
site="DpnII"
# genome ID, default to human, can also be set in options
genomeID="hg19"

## Read arguments                                                     
usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-s site] [-hx]"
genomeHelp="   genomeID must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \"$genomeID\")"
dirHelp="   [topDir] is the top level directory (default \"$topDir\") and must contain links to all merged_nodups files underneath it"
siteHelp="   [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" (default \"$site\"); alternatively, this can be the restriction site file"
excludeHelp="   -x: exclude fragment-delimited maps from Hi-C mega map (will run much faster)"
helpHelp="   -h: print this help and exit"

printHelpAndExit() {
    echo "$usageHelp"
    echo "$genomeHelp"
    echo "$dirHelp"
    echo "$siteHelp"
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:hxs:" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	s) site=$OPTARG ;;
  x) exclude=1 ;;
	[?]) printHelpAndExit 1;;
    esac
done

## Set ligation junction based on restriction enzyme
case $site in
	HindIII) ligation="AAGCTAGCTT";;
	DpnII) ligation="GATCGATC";;
	MboI) ligation="GATCGATC";;
  none) ligation="XXXX";;
	*)  ligation="XXXX"
      site_file=$site
      echo "$site not listed as recognized enzyme, so trying it as site file."
      echo "Ligation junction is undefined";;
esac

if [ -z "$site_file" ]
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [ ! -e "$site_file" ] && [ "$site" != "none" ]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    echo "The site file is used for statistics even if fragment delimited maps are excluded"
    exit 100
fi

## Directories to be created and regex strings for listing files
megadir=${topDir}"/mega"
outputdir=${megadir}"/aligned"
tmpdir=${megadir}"/HIC_tmp"
export TMPDIR=${megadir}"/HIC_tmp"
outfile=${megadir}/lsf.out
touchfile1=${megadir}/touch1

## Check for existing merge_nodups files:

merged_count=`find -L ${topDir} | grep merged_nodups.txt | wc -l`
if [ "$merged_count" -lt "1" ]
then
	echo "***! Failed to find at least one merged_nodups files under ${topDir}"
	exit 100
fi

merged_names=$(find -L ${topDir} | grep merged_nodups.txt | tr '\n' ' ')
inter_names=$(find -L ${topDir} | grep inter.txt | tr '\n' ' ')

## Create output directory, exit if already exists
if [[ -d "${outputdir}" ]] 
then
    echo "***! Move or remove directory \"${outputdir}\" before proceeding."
	exit 101			
else
	mkdir -p ${outputdir}
fi

## Create temporary directory
if [ ! -d "$tmpdir" ]; then
    mkdir $tmpdir
    chmod 777 $tmpdir
fi

#Prepare an empty log file:
touch ${outfile}
chmod 777 ${outfile}

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

#source $usePath
#$load_cluster

jid1=`bsub -o ${megadir}/lsf.out -q "${queue}" -J "${groupname}_topstats" <<-TOPSTATS

export LC_ALL=C
if ! awk -f ${juiceDir}/scripts/makemega_addstats.awk ${inter_names} > ${outputdir}/inter.txt
then  
    echo "***! Some problems occurred somewhere in creating top stats files."
    exit 100
else
    cp ${outputdir}/inter.txt ${outputdir}/inter_30.txt
fi
touch $touchfile1
TOPSTATS`

touchfile2=${megadir}/touch2
# Merge all merged_nodups.txt files found under current dir
jid2=`bsub -o ${megadir}/lsf.out -q "${long_queue}" -J ${groupname}_merge -R "rusage[mem=16000]" -w "done(${groupname}_topstats)" <<- MRGSRT
if [ ! -f "${touchfile1}" ]
then
   echo "***! Top stats job failed, type bjobs -l $jid1 to see what happened."
   exit 100;
fi
if ! sort -T ${tmpdir} -m -k2,2d -k6,6d ${merged_names} > ${outputdir}/merged_nodups.txt
then 
echo "***! Some problems occurred somewhere in creating sorted merged_nodups files."
    exit 100
else
echo "(-: Finished sorting all merged_nodups files into a single merge."
  rm -r ${tmpdir}
fi
touch $touchfile2
MRGSRT`

touchfile3=${megadir}/touch3  

# Create statistics files for MQ > 0
jid3=`bsub -o ${megadir}/lsf.out -q "${queue}" -J "${groupname}_inter0" -w "done(${groupname}_merge)" <<- INTER0
if [ ! -f "${touchfile2}" ]
then
   echo "***! Sort job failed."
   exit 100;
fi

${juiceDir}/scripts/statistics.pl -q 1 -o${outputdir}/inter.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt
touch $touchfile3
INTER0`

touchfile4=${megadir}/touch4

# Create statistics files for MQ > 30
jid4=`bsub -o ${megadir}/lsf.out -q "${queue}" -J "${groupname}_inter30" -w "done(${groupname}_merge)" <<- INTER30
if [ ! -f "${touchfile2}" ]
then
   echo "***! Sort job failed."
   exit 100;
fi

${juiceDir}/scripts/statistics.pl -q 30 -o${outputdir}/inter_30.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt 
touch $touchfile4
INTER30`
touchfile5=${megadir}/touch5

# Create HIC maps file for MQ > 0
jid5=`bsub -o ${megadir}/lsf.out -q "${queue}" -J "${groupname}_hic0" -w "done(${groupname}_inter30)" -R "rusage[mem=16000]" <<- HIC0
#source $usePath
#$load_java
if [ ! -f "${touchfile3}" ]
then
   echo "***! Statistics q=1 job failed."
   exit 100;
fi
exitcode=-999
if [ -z "$exclude" ]
then
	${juiceDir}/scripts/juicebox pre -f ${site_file} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
    exitcode=\$?
else
	${juiceDir}/scripts/juicebox pre -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
    exitcode=\$?
fi
if [ "\${exitcode}" -eq 0 ]
then
    touch $touchfile5
fi
HIC0`

touchfile6=${megadir}/touch6
# Create HIC maps file for MQ > 30
jid6=`bsub -o ${megadir}/lsf.out -q "${queue}" -J "${groupname}_hic30" -w "done(${groupname}_inter30)" -R "rusage[mem=16000]" <<- HIC30
#source $usePath
#$load_java	
if [ ! -f "${touchfile4}" ]
then
   echo "***! Statistics q=30 job failed."
   exit 100;
fi
exitcode=-999
if [ -z "${exclude}" ]
then
	${juiceDir}/scripts/juicebox pre -f ${site_file} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}
   exitcode=\$?
else
	${juiceDir}/scripts/juicebox pre -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}
   exitcode=\$?
fi
if [ "\${exitcode}" -eq 0 ]
then
touch $touchfile6
fi
HIC30`

touchfile7=${megadir}/touch7
# Create loop and domain lists file for MQ > 30
jid7=`bsub  -o ${megadir}/lsf.out -q "${queue}" -J "${groupname}_postproc" -w "done(${groupname}_hic30)" <<- POSTPROC
#source $usePath
#$load_java	
export _JAVA_OPTIONS=-Xmx16384m;
export LC_ALL=C
if [ ! -f "${touchfile6}" ]
then
   echo "***! Failed to create domain and loop lists"
   exit 100;
fi
${juiceDir}/scripts/juicer_postprocessing.sh -j ${juiceDir}/scripts/juicebox -i ${outputdir}/inter_30.hic -m ${juiceDir}/references/motif -g ${genomeID}
touch $touchfile7
POSTPROC`

jid8=`bsub -o ${megadir}/lsf.out -j y -q "${queue}" -J "${groupname}_done" -w "done(${groupname}_postproc)" <<- FINAL
if [ ! -f "${touchfile5}" ]
then
   echo "***! Failed to make inter.hic."   
   exit 100;
fi
if [ ! -f "${touchfile6}" ]
then
   echo "***! Failed to make inter_30.hic."   
   exit 100;
fi
if [ ! -f "${touchfile7}" ]
then
   echo "***! Failed in postprocessing"   
   exit 100;
fi

rm $touchfile1 $touchfile2 $touchfile3 $touchfile4 $touchfile5 $touchfile6 $touchfile7
echo "(-: Successfully completed making mega map. Done. :-)"
FINAL`

echo "(-: Finished adding all jobs... please wait while processing."
