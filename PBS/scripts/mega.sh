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
load_java="module load java/jdk1.8.0_131"
load_cuda="module load cuda/7.5.18/gcc/4.4.7"
# Juicer directory, contains scripts/ and restriction_sites/
juiceDir="/lustre1/mzhibo/hic/apps/juicer"
# default queue and time
#queue="batch"
#walltime="walltime=12:00:00"

# unique name for jobs in this run
groupname="C$(date +"%s"|cut -c 6-11)"


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

while getopts "d:g:r:h:x:s" opt; do
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
echo $site
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
echo $ligation
if [ -z "$site_file" ]
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi
echo $site_file
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
# set global tmpdir so no problems with /var/tmp
export TMPDIR=${tmpdir}
#output messages
logdir=${megadir}"/debug"

## Check for existing merge_nodups files:

merged_count=`find -L ${topDir} | grep merged_nodups.txt | wc -l`
if [ "$merged_count" -lt "1" ]
then
	echo "***! Failed to find at least one merged_nodups files under ${topDir}"
	exit 100
fi

merged_names1=$(find -L ${topDir} | grep merged_nodups.txt)
merged_names=$(echo $merged_names1 | tr '\n' ' ')
inter_names=$(find -L ${topDir} | grep inter.txt | tr '\n' ' ')

if [[ $merged_names == *".txt.gz"* ]]
then
    gzipped=1
    echo "***! Mega map of gzipped files not yet supported, please unzip before running."
    exit 100
    # we need to unzip here
    for i in $merged_names1
    do
	if [[ $i != *".txt.gz"* ]]
	then
	    echo "***! Mixture of gzipped and unzipped merged_nodups files"
	    echo "Ensure that the merged_nodups are all either unzipped or gzipped then rerun"
	    echo "Files: $merged_names"
	    exit 100
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
qsub -o ${logdir}/header.log -j oe -q batch -N ${groupname}_cmd <<-EOF
  date
  echo "Juicer version:$juicer_version"
  echo "$0 $@"
EOF

jid1=$(qsub -o ${logdir}/topstats.log -j oe -N ${groupname}_Tstats -l mem=20gb -l walltime=24:00:00 -l nodes=1:ppn=1:AMD -q batch <<-TOPSTATS
export LC_ALL=C
if ! awk -f ${juiceDir}/scripts/makemega_addstats.awk ${inter_names} > ${outputdir}/inter.txt
then  
    echo "***! Some problems occurred somewhere in creating top stats files."
    exit 100
else
    cp ${outputdir}/inter.txt ${outputdir}/inter_30.txt
fi
TOPSTATS
)
jobIDstr=${jid1}
# Merge all merged_nodups.txt files found under current dir
jid2=$(qsub -o ${logdir}/merge.log -j oe -q batch -N ${groupname}_merge -l mem=20gb -l walltime=24:00:00 -l nodes=1:ppn=1:AMD <<- MRGSRT
if ! sort -T ${tmpdir} -m -k2,2d -k6,6d ${merged_names} > ${outputdir}/merged_nodups.txt
then 
echo "***! Some problems occurred somewhere in merging sorted merged_nodups files."
    exit 100
else
echo "Finished sorting all merged_nodups files into a single merge."
  rm -r ${tmpdir}
fi
MRGSRT
)
jobIDstr="${jobIDstr}:${jid2}"

# Create statistics files for MQ > 0
jid3=$(qsub -o ${logdir}/inter0.log -j oe -q batch -N ${groupname}_inter0 -l mem=20gb -l walltime=24:00:00 -l nodes=1:ppn=1:AMD -W depend=afterok:${jid1}:${jid2} <<- INTER0
${juiceDir}/scripts/statistics.pl -q 1 -o${outputdir}/inter.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt
INTER0
)
jobIDstr="${jobIDstr}:${jid3}"

# Create statistics files for MQ > 30
jid4=$(qsub -o ${logdir}/inter30.log -j oe -q batch -N ${groupname}_inter30 -l mem=20gb -l walltime=24:00:00 -l nodes=1:ppn=1:AMD -W depend=afterok:${jid1}:${jid2}  <<- INTER30
${juiceDir}/scripts/statistics.pl -q 30 -o${outputdir}/inter_30.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt 
INTER30
)
jobIDstr="${jobIDstr}:${jid4}"

# Create HIC maps file for MQ > 0
jid5=$(qsub -o ${logdir}/hic0_${groupname}.log -j oe -q batch -M ${EMAIL} -m ae -N ${groupname}_hic0 -l mem=40gb -l walltime=168:00:00 -l nodes=1:ppn=1:AMD -W depend=afterok:${jid4} <<- HIC0
$load_java
if [ -z "$exclude" ]
then
  echo "Launching ${juiceDir}/scripts/juicer_tools pre ${resolutions} -f ${site_file} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}"
  ${juiceDir}/scripts/juicer_tools pre ${resolutions} -f ${site_file} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
else
  echo "Launching ${juiceDir}/scripts/juicer_tools pre ${resolutions} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}"
  ${juiceDir}/scripts/juicer_tools pre ${resolutions} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
fi
HIC0
)
jobIDstr="${jobIDstr}:${jid5}"
# Create HIC maps file for MQ > 30
jid6=$(qsub -o ${logdir}/hic30_${groupname}.log -j oe -q batch -M ${EMAIL} -m ae -N ${groupname}_hic30 -l mem=60gb -l walltime=168:00:00 -l nodes=1:ppn=1:AMD -W depend=afterok:${jid4} <<- HIC30
$load_java
if [ -z "${exclude}" ]
then
    echo "Launching ${juiceDir}/scripts/juicer_tools pre ${resolutions} -f ${site_file} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}"
    ${juiceDir}/scripts/juicer_tools pre ${resolutions} -f ${site_file} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}
else
   echo "Launching ${juiceDir}/scripts/juicer_tools pre ${resolutions} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}"
   ${juiceDir}/scripts/juicer_tools pre ${resolutions} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}
fi
HIC30
)
jobIDstr="${jobIDstr}:${jid6}"

# Create loop and domain lists file for MQ > 30
jid7=$(qsub  -o ${logdir}/hiccups.log -j oe -q batch -N ${groupname}hiccups -l mem=60gb -l walltime=100:00:00 -l nodes=1:ppn=1:gpus=1:GPU -W depend=afterok:${jid6} <<- HICCUPS
$load_java
$load_cuda
echo $PBS_GPUFILE
export _JAVA_OPTIONS=-Xmx16384m;
export LC_ALL=C
${juiceDir}/scripts/juicer_hiccups.sh -j ${juiceDir}/scripts/juicer_tools -i ${outputdir}/inter_30.hic -m ${juiceDir}/references/motif -g ${genomeID}
B
HICCUPS

)

jobIDstr="${jobIDstr}:${jid7}"

jid8=$(qsub -o ${logdir}/arrowhead.log -j oe -q batch -N ${groupname}_ArwHead -l mem=60gb -l walltime=100:00:00 -l nodes=1:ppn=1:gpus=1:GPU -W depend=afterok:${jid6} <<- ARROWHEAD
$load_java
$load_cuda
echo $PBS_GPUFILE
export _JAVA_OPTIONS=-Xmx16384m;
export LC_ALL=C
${juiceDir}/scripts/juicer_arrowhead.sh -j ${juiceDir}/scripts/juicer_tools -i ${outputdir}/inter_30.hic
ARROWHEAD
)

qsub -o ${logdir}/done.log -j oe -q batch -N ${groupname}_done -W depend=afterok:${jobIDstr} <<- FINAL
echo "All jobs finished processing!"
FINAL
qsub -o ${logdir}/done.out -j oe -q batch -N ${groupname}_fail -W depend=afternotok:${jobIDstr} <<- FINAL
echo "Error occored in placing the jobs. Please check err file of each step to find out"
FINAL
