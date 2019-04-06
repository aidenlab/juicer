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
usePath=/broad/software/scripts/useuse
load_java="use Java-1.7"
load_cluster="use UGER"
# Juicer directory, contains scripts/ and restriction_sites/
juiceDir="/broad/aidenlab"
#default queue
queue="long"
# unique name for jobs in this run
groupname="a$(date +%s)"

## Default options, overridden by command line arguments

# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
site="MboI"
# genome ID, default to human, can also be set in options
genomeID="hg19"

## Read arguments                                                     
usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-s site] [-r resolutions] [-hx]"
genomeHelp="   genomeID must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \"$genomeID\")"
dirHelp="   [topDir] is the top level directory (default \"$topDir\") and must contain links to all merged_nodups files underneath it"
siteHelp="   [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" (default \"$site\"); alternatively, this can be the restriction site file"
resolutionsHelp="   [resolutions] is a comma-delimited list of resolutions, such as 10000,5000,1000,5f (default is 2.5M,1M,500K,250K,100K,50K,25K,10K,5K in base pair and 500f,250f,100f,50f,25f,10f,5f,2f,1f)"
excludeHelp="   -x: exclude fragment-delimited maps from Hi-C mega map (will run much faster)"
helpHelp="   -h: print this help and exit"

printHelpAndExit() {
    echo "$usageHelp"
    echo "$genomeHelp"
    echo "$dirHelp"
    echo "$siteHelp"
    echo "$resolutionsHelp"
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:r:hxs:" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	s) site=$OPTARG ;;
	x) exclude=1 ;;
	r) resolutions=$OPTARG ;;
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
    exit 1
fi

## Directories to be created and regex strings for listing files
megadir=${topDir}"/mega"
outputdir=${megadir}"/aligned"
tmpdir=${megadir}"/HIC_tmp"
# set global tmpdir so no problems with /var/tmp
export TMPDIR=${tmpdir}
#output messages
logdir=${megadir}/debug

## Check for existing merge_nodups files:

merged_count=`find -L ${topDir} | grep merged_nodups.txt | wc -l`
if [ "$merged_count" -lt "1" ]
then
	echo "***! Failed to find at least one merged_nodups files under ${topDir}"
	exit 1
fi

merged_names1=$(find -L ${topDir} | grep merged_nodups.txt)
merged_names=$(echo $merged_names1 | tr '\n' ' ')
inter_names=$(find -L ${topDir} | grep inter.txt | tr '\n' ' ')

if [[ $merged_names == *".txt.gz"* ]]
then
    gzipped=1
    echo "***! Mega map of gzipped files not yet supported, please unzip before running."
    exit 1
    # we need to unzip here
    for i in $merged_names1
    do
	if [[ $i != *".txt.gz"* ]]
	then
	    echo "***! Mixture of gzipped and unzipped merged_nodups files"
	    echo "Ensure that the merged_nodups are all either unzipped or gzipped then rerun"
	    echo "Files: $merged_names"
	    exit 1
	fi
    done
fi

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

## Create log directory
if [ ! -d "$logdir" ]; then
    mkdir $logdir
    chmod 777 $logdir
fi

if [ -n "$resolutions" ]; then
    resolutions="-r $resolutions"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline
qsub -o ${logdir}/header.out -j y -r y -N ${groupname}cmd -cwd <<-EOF
  date
  echo "Juicer version:$juicer_version"
  echo "$0 $@"
EOF

source $usePath
$load_cluster
touchfile1=${megadir}/touch1
jid1=`qsub -terse -o ${logdir}/topstats.out -e ${logdir}/topstats.err -r y -N ${groupname}_topstats <<-TOPSTATS
export LC_ALL=en_US.UTF-8
if ! awk -f ${juiceDir}/scripts/makemega_addstats.awk ${inter_names} > ${outputdir}/inter.txt
then  
    echo "***! Some problems occurred somewhere in creating top stats files."
    exit 1
else
    echo "(-: Finished creating top stats files."
    cp ${outputdir}/inter.txt ${outputdir}/inter_30.txt
fi
touch $touchfile1
TOPSTATS`
echo "topstats job id: $jid1" > $logdir/jobs.out

touchfile2=${megadir}/touch2
# Merge all merged_nodups.txt files found under current dir
jid2=`qsub -terse -o ${logdir}/merge.out -e ${logdir}/merge.err -q ${queue} -r y -N ${groupname}_merge -l h_vmem=16g -hold_jid ${groupname}_topstats <<- MRGSRT
if [ -n "$gzipped" ]
then 
  # This code didn't work, doesn't sort appropriately unfortunately.
  if ! gunzip -c $merged_names | sort -T ${tmpdir} -m -k2,2d -k6,6d > ${outputdir}/merged_nodups.txt 
  then
    echo "***! Some problems occurred somewhere in merging sorted merged_nodups files."
    exit 1
  else
    echo "(-: Finished sorting all merged_nodups files into a single merge."
    rm -r ${tmpdir}
  fi
elif ! sort -T ${tmpdir} -m -k2,2d -k6,6d ${merged_names} > ${outputdir}/merged_nodups.txt
then 
  echo "***! Some problems occurred somewhere in merging sorted merged_nodups files."
  exit 1
else
  echo "(-: Finished sorting all merged_nodups files into a single merge."
  rm -r ${tmpdir}
fi
touch $touchfile2
MRGSRT`
echo "merge job id: $jid2" >> $logdir/jobs.out

touchfile3=${megadir}/touch3
holdjobs1="-hold_jid ${groupname}_merge,${groupname}_topstats";    

# Create statistics files for MQ > 0
jid3=`qsub -terse -o ${logdir}/inter0.out -e ${logdir}/inter0.err -q ${queue} -r y -N ${groupname}_inter0 $holdjobs1 <<- INTER0
if [ ! -f "${touchfile1}" ]
then
   echo "***! Top stats job failed, type qacct -j $jid1 to see what happened."
   exit 1
fi
if [ ! -f "${touchfile2}" ]
then
   echo "***! Sort job failed, type qacct -j $jid2 to see what happened."
   exit 1
fi

${juiceDir}/scripts/statistics.pl -q 1 -o${outputdir}/inter.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt
touch $touchfile3
INTER0`
echo "statistics mapq0 job id: $jid3" >> $logdir/jobs.out

touchfile4=${megadir}/touch4
holdjobs2="-hold_jid ${groupname}_inter0";    

# Create statistics files for MQ > 30
jid4=`qsub -terse -o ${logdir}/inter30.out -e ${logdir}/inter30.err -q ${queue} -r y -N ${groupname}_inter30 $holdjobs1 <<- INTER30
if [ ! -f "${touchfile1}" ]
then
   echo "***! Top stats job failed, type qacct -j $jid1 to see what happened."
   exit 1
fi
if [ ! -f "${touchfile2}" ]
then
   echo "***! Sort job failed, type qacct -j $jid2 to see what happened."
   exit 1
fi

${juiceDir}/scripts/statistics.pl -q 30 -o${outputdir}/inter_30.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt 
touch $touchfile4
INTER30`
echo "statistics mapq30 job id: $jid4" >> $logdir/jobs.out

holdjobs3="-hold_jid ${groupname}_inter30";
touchfile5=${megadir}/touch5

# Create HIC maps file for MQ > 0
jid5=`qsub -terse -o ${logdir}/hic0.out -e ${logdir}/hic0.err -q ${queue} -r y -N ${groupname}_hic0  -l h_vmem=16g $holdjobs2 <<- HIC0
source $usePath
$load_java
if [ ! -f "${touchfile3}" ]
then
   echo "***! Statistics q=1 job failed, type qacct -j $jid3 to see what happened."
   exit 1
fi
myexitcode=1
if [ -z "$exclude" ]
then
    echo "Launching ${juiceDir}/scripts/juicer_tools pre $resolutions -f ${site_file} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}"
    if ${juiceDir}/scripts/juicer_tools pre $resolutions -f ${site_file} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
    then
       touch $touchfile5
    fi
else
    echo "Launching ${juiceDir}/scripts/juicer_tools pre $resolutions -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}"
    if ${juiceDir}/scripts/juicer_tools pre $resolutions -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
    then
       touch $touchfile5
    fi
fi
HIC0`
echo "hic mapq0 job id: $jid5" >> $logdir/jobs.out

touchfile6=${megadir}/touch6
# Create HIC maps file for MQ > 30
jid6=`qsub -terse -o ${logdir}/hic30.out -e ${logdir}/hic30.err -q ${queue} -r y -N ${groupname}_hic30 -l h_vmem=16g $holdjobs3 <<- HIC30
source $usePath
$load_java	
if [ ! -f "${touchfile4}" ]
then
   echo "***! Statistics q=30 job failed, type qacct -j $jid4 to see what happened."
   exit 1
fi
myexitcode=1
if [ -z "${exclude}" ]
then
    echo "Launching ${juiceDir}/scripts/juicer_tools pre $resolutions -f ${site_file} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}"
    if ${juiceDir}/scripts/juicer_tools pre $resolutions -f ${site_file} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}
    then
       touch $touchfile6
    fi
else
    echo "Launching ${juiceDir}/scripts/juicer_tools pre $resolutions -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}"
    if ${juiceDir}/scripts/juicer_tools pre $resolutions -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}
    then
       touch $touchfile6
    fi
fi
HIC30`
echo "hic mapq30 job id: $jid6" >> $logdir/jobs.out

touchfile7=${megadir}/touch7
jid7=`qsub -terse -o ${logdir}/hiccups.out -e ${logdir}/hiccups.err -q ${queue} -r y -N ${groupname}_hiccups  -l h_vmem=16g -hold_jid ${groupname}_hic30 <<- HICCUPS
source $usePath;
$load_java;
export _JAVA_OPTIONS=-Xmx16384m;
export LC_ALL=en_US.UTF-8;
if [ ! -f "${touchfile6}" ]
then
   echo "***! Failed to make inter_30.hic, type qacct -j $jid6 to see what happened."
   exit 1
fi
${juiceDir}/scripts/juicer_hiccups.sh -j ${juiceDir}/scripts/juicer_tools -i ${outputdir}/inter_30.hic -m ${juiceDir}/references/motif -g ${genomeID}
touch $touchfile7
HICCUPS`
echo "hiccups job id: $jid7" >> $logdir/jobs.out

touchfile8=${megadir}/touch8
jid8=`qsub -terse -o ${logdir}/arrowhead.out -e ${logdir}/arrowhead.err -q ${queue} -r y -N ${groupname}_arrowhead  -l h_vmem=4g -hold_jid ${groupname}_hic30 <<- ARROWHEAD
source $usePath;
$load_java;
export _JAVA_OPTIONS=-Xmx16384m;
export LC_ALL=en_US.UTF-8;
if [ ! -f "${touchfile6}" ]
then
   echo "***! Failed to make inter_30.hic, type qacct -j $jid6 to see what happened."
   exit 1
fi
${juiceDir}/scripts/juicer_arrowhead.sh -j ${juiceDir}/scripts/juicer_tools -i ${outputdir}/inter_30.hic
touch $touchfile8
ARROWHEAD`
echo "arrowhead job id: $jid8" >> $logdir/jobs.out

holdjobs="-hold_jid ${groupname}_hic0,${groupname}_hic30,${groupname}_hiccups,${groupname}_arrowhead";

qsub -o ${logdir}/done.out -e ${logdir}/done.err -r y -N ${groupname}_done $holdjobs <<- FINAL
if [ ! -f "${touchfile5}" ]
then
   echo "***! Failed to make inter.hic, type qacct -j $jid5 to see what happened."   
   exit 1
fi
if [ ! -f "${touchfile6}" ]
then
   echo "***! Failed to make inter_30.hic, type qacct -j $jid6 to see what happened."   
   exit 1
fi
if [ ! -f "${touchfile7}" ]
then
   echo "***! Failed to create loop lists, type qacct -j $jid7 to see what happened."   
   exit 1
fi
if [ ! -f "${touchfile8}" ]
then
   echo "***! Failed to create contact domain lists, type qacct -j $jid8 to see what happened."   
   exit 1
fi
rm -r ${tmpdir}
rm $touchfile1 $touchfile2 $touchfile3 $touchfile4 $touchfile5 $touchfile6 $touchfile7 $touchfile8
echo "(-: Successfully completed making mega map. Done. :-)"
FINAL

echo "(-: Finished adding all jobs... please wait while processing."
