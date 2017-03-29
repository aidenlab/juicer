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
load_java="module load Java/7.1.2.0"
#load_cluster=""
#load_coreutils=""
load_gpu="module load gcccuda/2.6.10;export PATH=\$PATH:/usr/local/cuda/bin/;"
# Juicer directory, contains scripts/ and restriction_sites/
juiceDir="/poscratch/aidenlab/juicer"
# default queue, can also be set in options
queue="commons"
# default long queue, can also be set in options
long_queue="commons"

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
usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-s site] [-h]"
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
export TMPDIR=${tmpdir}
outfile=${megadir}/lsf.out
#output messages
logdir="$megadir/debug"
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

## Create output directory, used for reporting commands output
if [ ! -d "$logdir" ]; then
        mkdir "$logdir"
        chmod 777 "$logdir"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline
#source $usePath
#$load_cluster

# Create top statistics file from all inter.txt files found under current dir
touchfile1=${megadir}/touch1
jid1=`sbatch <<- TOPSTATS | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p $queue
#SBATCH -t 1440
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH -o $logdir/topstats-%j.out
#SBATCH -e $logdir/topstats-%j.err
#SBATCH -J "${groupname}_topstats"
export LC_COLLATE=C
if ! awk -f ${juiceDir}/scripts/makemega_addstats.awk ${inter_names} > ${outputdir}/inter.txt
then  
    echo "***! Some problems occurred somewhere in creating top stats files."
    exit 100
else
    echo "(-: Finished creating top stats files."
    cp ${outputdir}/inter.txt ${outputdir}/inter_30.txt
fi
touch $touchfile1
TOPSTATS`
dependtopstats="afterok:$jid1"

# Merge all merged_nodups.txt files found under current dir
touchfile2=${megadir}/touch2
jid2=`sbatch <<- MRGSRT | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p ${long_queue}
#SBATCH -t 1440
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16384
#SBATCH -o $logdir/merge-%j.out
#SBATCH -e $logdir/merge-%j.err
#SBATCH -J "${groupname}_merge"
#SBATCH -d "${dependtopstats}"
if [ ! -f "${touchfile1}" ]
then
    echo "***! Top stats job failed, type \"scontrol show job $jid1\" to see what happened."
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
dependmerge="afterok:$jid2"

# Create statistics files for MQ > 0
touchfile3=${megadir}/touch3  
jid3=`sbatch <<- INTER0 | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p ${long_queue}
#SBATCH -t 1440
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH -o $logdir/inter0-%j.out
#SBATCH -e $logdir/inter0-%j.err
#SBATCH -J "${groupname}_inter0"
#SBATCH -d "${dependmerge}"
if [ ! -f "${touchfile2}" ]
then
   echo "***! Sort job failed."
   exit 100;
fi
${juiceDir}/scripts/statistics.pl -q 1 -o${outputdir}/inter.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt
touch $touchfile3
INTER0`
dependinter0="afterok:$jid3"

# Create statistics files for MQ > 30
touchfile4=${megadir}/touch4
jid4=`sbatch <<- INTER30 | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p ${long_queue}
#SBATCH -t 1440
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH -o $logdir/inter30-%j.out
#SBATCH -e $logdir/inter30-%j.err
#SBATCH -J "${groupname}_inter30"
#SBATCH -d "${dependmerge}"
if [ ! -f "${touchfile2}" ]
then
   echo "***! Sort job failed."
   exit 100;
fi
${juiceDir}/scripts/statistics.pl -q 30 -o${outputdir}/inter_30.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt 
touch $touchfile4
INTER30`
dependinter30="afterok:$jid4"

# Create HIC maps file for MQ > 0
touchfile5=${megadir}/touch5
jid5=`sbatch <<- HIC0 | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p ${long_queue}
#SBATCH -t 1440
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH -o $logdir/hic0-%j.out
#SBATCH -e $logdir/hic0-%j.err
#SBATCH -J "${groupname}_hic0"
#SBATCH -d "${dependinter0}"
#SBATCH --mem-per-cpu=32G
#source $usePath
$load_java
if [ ! -f "${touchfile3}" ]
then
   echo "***! Statistics q=1 job failed."
   exit 100;
fi
if [ -z "$exclude" ]
then
	${juiceDir}/scripts/juicer_tools pre -f ${site_file} -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
    exitcode=\$?
else
	${juiceDir}/scripts/juicer_tools pre -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
    exitcode=\$?
fi
if [ "\${exitcode}" -eq 0 ]
then
    touch $touchfile5
fi
HIC0`
dependhic0="afterok:$jid5"

# Create HIC maps file for MQ > 30
touchfile6=${megadir}/touch6
jid6=`sbatch <<- HIC30 | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p ${long_queue}
#SBATCH -t 1440
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH -o $logdir/hic30-%j.out
#SBATCH -e $logdir/hic30-%j.err
#SBATCH -J "${groupname}_hic30"
#SBATCH -d "${dependinter30}"
#SBATCH --mem-per-cpu=32G
#source $usePath
$load_java	
if [ ! -f "${touchfile4}" ]
then
   echo "***! Statistics q=30 job failed."
   exit 100;
fi
if [ -z "${exclude}" ]
then
	${juiceDir}/scripts/juicer_tools pre -f ${site_file} -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}
   exitcode=\$?
else
	${juiceDir}/scripts/juicer_tools pre -s ${outputdir}/inter_30.txt -g ${outputdir}/inter_30_hists.m -q 30 ${outputdir}/merged_nodups.txt ${outputdir}/inter_30.hic ${genomeID}
   exitcode=\$?
fi
if [ "\${exitcode}" -eq 0 ]
then
    touch $touchfile6
fi
HIC30`
dependhic30="${dependhic0}:$jid6"
dependhic30only="afterok:$jid6"

# Create loop lists file for MQ > 30
touchfile7=${megadir}/touch7
jid7=`sbatch <<- HICCUPS | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p ${long_queue}
#SBATCH -t 1440
#SBATCH -c 2
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G 
#SBATCH --gres=gpu:kepler:1
#SBATCH -o $logdir/hiccups-%j.out
#SBATCH -e $logdir/hiccups-%j.err
#SBATCH -J "${groupname}_hiccups"
#SBATCH -d "${dependhic30only}"
#source $usePath
$load_java	
if [ ! -f "${touchfile6}" ]
then
   echo "***! HIC maps q=30 job failed."
   exit 100;
fi
${load_gpu}
${juiceDir}/scripts/juicer_hiccups.sh -j ${juiceDir}/scripts/juicer_tools -i $outputdir/inter_30.hic -m ${juiceDir}/references/motif -g $genomeID
touch $touchfile7
HICCUPS`
dependhiccups="afterok:$jid7"

touchfile8=${megadir}/touch8
# Create domain lists for MQ > 30
jid8=`sbatch <<- ARROWHEAD | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p ${long_queue}
#SBATCH -t 1440
#SBATCH -c 2
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G 
#SBATCH -o $logdir/arrowhead-%j.out
#SBATCH -e $logdir/arrowhead-%j.err
#SBATCH -J "${groupname}_arrowhead"
#SBATCH -d "${dependhic30only}"
#source $usePath
$load_java	
if [ ! -f "${touchfile6}" ]
then
   echo "***! HIC maps q=30 job failed."
   exit 100;
fi
${juiceDir}/scripts/juicer_arrowhead.sh -j ${juiceDir}/scripts/juicer_tools -i $outputdir/inter_30.hic
touch $touchfile8
ARROWHEAD`
dependarrows="${dependhic0}:$jid8"
dependarrowsonly="afterok:$jid8"

# Final checks
jid9=`sbatch <<- FINAL | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p ${long_queue}
#SBATCH -t 100
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH -o $logdir/done-%j.out
#SBATCH -e $logdir/done-%j.err
#SBATCH -J "${groupname}_done"
#SBATCH -d "${dependarrows}"

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
   echo "***! Failed to create loop lists."   
   exit 100;
fi
if [ ! -f "${touchfile8}" ]
then
   echo "***! Failed to create domain lists."   
   exit 100;
fi
rmdir ${tmpdir}
rm $touchfile1 $touchfile2 $touchfile3 $touchfile4 $touchfile5 $touchfile6 $touchfile7 $touchfile8
echo "(-: Successfully completed making mega map. Done. :-)"
echo $topDir, $site, $genomeID, $genomePath | mail -r MegaJuicer@rice.edu -s \"Mega Juicer pipeline finished successfully @ Rice" -t $USER@rice.edu;
FINAL`
dependfinal="afterok:$jid9"
echo "(-: Finished adding all jobs... please wait while processing."
