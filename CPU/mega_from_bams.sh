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
#  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# MegaMap script. 
#
#                                                                       
# [topDir]    - Should contain the results of all base experiments
#
# From the top-level directory, the following two directories are created:
#
# [topDir]/mega     - Location of result of processing the mega map
#
# Juicer version 2.0
juicer_version="2.0"
# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
site="none"
# genome ID, default to human, can also be set in options
genomeID="hg19"
# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir="/aidenlab"
# by default exclude fragment delimited maps
exclude=1
# single-end input, default no
singleend=0

usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-s site] [-S stage] [-D Juicer scripts directory] [-T threadsHic] [-y site_file] [-f] [-h]"
genomeHelp="   genomeID is either defined in the script, e.g. \"hg19\" or \"mm10\" or the path to the chrom.sizes file"
dirHelp="   [topDir] is the top level directory (default \"$topDir\") and must contain links to all bam files underneath it"
siteHelp="   [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" (default \"$site\"); alternatively, this can be the restriction site file"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
threadsHicHelp="* [threads for hic file creation]: number of threads when building hic file"
stageHelp="* [stage]: must be one of \"final\", \"postproc\", or \"early\".\n    -Use \"final\" when the reads have been combined into merged but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the stats and\n     hic files"
singleEndHelp="* -u: Single end alignment"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
excludeHelp="   -f: include fragment-delimited maps from Hi-C mega map (will run slower)"
helpHelp="   -h: print this help and exit"

printHelpAndExit() {
    echo "Juicer Version: ${juicer_version}"
    echo "$usageHelp"
    echo "$genomeHelp"
    echo "$dirHelp"
    echo "$siteHelp"
    echo "$siteFileHelp"
    echo "$stageHelp"
    echo "$scriptDirHelp"
    echo "$threadsHicHelp"
    echo "$excludeHelp"
    echo "$singleEndHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:hfs:S:D:y:T:u:" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	s) site=$OPTARG ;;
	f) exclude=0 ;;
	y) site_file=$OPTARG ;;
	S) stage=$OPTARG ;;
	D) juiceDir=$OPTARG ;;
	T) threadsHic=$OPTARG ;;
  u) singleend=1 ;;
	[?]) printHelpAndExit 1;;
    esac
done

## If DNAse-type experiment, no fragment maps; or way to get around site file
if [[ "$site" == "none" ]]
then
    exclude=1;
fi

if [ -z "$site_file" ]
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number
if [[ ! -e "$site_file" ]] && [[ "$site" != "none" ]] && [[ ! "$site_file" =~ "none" ]]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    echo "The site file is used for statistics even if fragment delimited maps are excluded"
    exit 1
elif [[ "$site" != "none" ]] && [[ ! "$site_file" =~ "none" ]]
then
    echo  "Using $site_file as site file"
fi

resolutionsToBuildString="-r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100,50,20,10"

if [ "$exclude" -eq 1 ]
then
  buildFragmentMapString=""
else
  buildFragmentMapString="-f $site_file"
fi

if [ -n "$stage" ]
then
    case $stage in
        final) final=1 ;;
	early) early=1 ;;
	postproc) postproc=1 ;;
        *)  echo "$usageHelp"
	    echo "$stageHelp"
	    exit 1
    esac
fi

## Directories to be created and regex strings for listing files
megaDir=${topDir}"/mega"
outputDir=${megaDir}"/aligned"
tmpdir=${megaDir}"/HIC_tmp"
export TMPDIR=${tmpdir}
tempdirPre=${outputDir}"/HIC_tmp"
tempdirPre30=${outputDir}"/HIC30_tmp"

if [ -z "$threadsHic" ]
then
  threadsHic=1
	threadHicString=""
	threadHic30String=""
	threadNormString=""
else
  threadHicString="--threads $threadsHic -i ${outputDir}/merged1_index.txt -t ${tempdirPre}"
	threadHic30String="--threads $threadsHic -i ${outputDir}/merged30_index.txt -t ${tempdirPre30}"
	threadNormString="--threads $threadsHic"
fi

cThreads="$(getconf _NPROCESSORS_ONLN)"
cThreadString="-@ $cThreads"

## Check for existing merged files:
merged_count=$(find -L "${topDir}" | grep -c merged_dedup.bam)
if [ "$merged_count" -lt "2" ]
then
    echo "***! Failed to find at least two merged_dedup.bam files under ${topDir}"
    exit 1
fi

bams_to_merge=$(find -L "${topDir}" | grep merged_dedup.bam | tr '\n' ' ')

## Create output directory, exit if already exists
if [[ -d "${outputDir}" ]] && [ -z $final ] && [ -z $postproc ]
then
    echo "***! Move or remove directory \"${outputDir}\" before proceeding."
    exit 1
else
    mkdir -p "${outputDir}"
fi

## Create temporary directory
if [ ! -d "$tmpdir" ]; then
    mkdir "$tmpdir"
    chmod 777 "$tmpdir"
fi



## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

# Not in final or postproc
if [ -z $final ] && [ -z $postproc ]
then
    samtools merge -c -t cb "$cThreadString" -o "${outputDir}"/mega_merged_dedup.bam "${bams_to_merge}"
    samtools view "$cThreadString" -F 1024 -O sam "${outputDir}"/mega_merged_dedup.bam | awk -v mapq=1 -f "${juiceDir}"/scripts/common/sam_to_pre.awk > "${outputDir}"/merged1.txt
    samtools view "$cThreadString" -F 1024 -O sam "${outputDir}"/mega_merged_dedup.bam | awk -v mapq=30 -f "${juiceDir}"/scripts/common/sam_to_pre.awk > "${outputDir}"/merged30.txt

    # Create statistics file
    if [ $singleend -eq 1 ]
    then
        ret=$(samtools view "$cThreadString" -f 1024 -F 256 "${outputDir}"/mega_merged_dedup.bam | awk '{if ($0~/rt:A:7/){singdup++}else{dup++}}END{print dup,singdup}')
        dups=$(echo "$ret" | awk '{print $1}')
        singdups=$(echo "$ret" | awk '{print $2}')
        cat "$splitdir"/*.res.txt | awk -v dups="$dups" -v singdups="$singdups" -v ligation="XXXX" -v singleend=1 -f "${juiceDir}"/scripts/common/stats_sub.awk >> "$outputDir"/inter.txt
    else
        dups=$(samtools view -c -f 1089 -F 256 "$cThreadString" "${outputDir}"/mega_merged_dedup.bam)
        cat "$splitdir"/*.res.txt | awk -v dups="$dups" -v ligation="XXXX" -f "${juiceDir}"/scripts/common/stats_sub.awk >> "$outputDir"/inter.txt
    fi
    cp "$outputDir"/inter.txt "$outputDir"/inter_30.txt

    echo "(-: Finished creating top stats files."
    cp "${outputDir}"/inter.txt "${outputDir}"/inter_30.txt
    echo "(-: Finished sorting all files into a single merge."
    "${juiceDir}"/scripts/common/juicer_tools statistics "$site_file" "$outputDir"/inter.txt "$outputDir"/merged1.txt "$genomeID"
    "${juiceDir}"/scripts/common/juicer_tools statistics "$site_file" "$outputDir"/inter_30.txt "$outputDir"/merged30.txt "$genomeID"

    mkdir "${tempdirPre}"
	  if [[ $threadsHic -gt 1 ]] && [[ ! -s "${outputDir}"/merged1_index.txt ]]
	  then
	      "${juiceDir}"/scripts/common/index_by_chr.awk "${outputDir}"/merged1.txt 500000 > "${outputDir}"/merged1_index.txt
	  fi
	  "${juiceDir}"/scripts/common/juicer_tools pre -n -s "$outputDir"/inter.txt -g "$outputDir"/inter_hists.m -q 1 "$resolutionsToBuildString" "$buildFragmentMapString" "$threadHicString" "$outputDir"/merged1.txt "$outputDir"/inter.hic "$genomeID"
	  "${juiceDir}"/scripts/common/juicer_tools addNorm "$threadNormString" "${outputDir}"/inter.hic
	  rm -Rf "${tempdirPre}"

	  mkdir "${tempdirPre30}"
	  if [[ $threadsHic -gt 1 ]] && [[ ! -s "${outputDir}"/merged30_index.txt ]]
	  then
	      "${juiceDir}"/scripts/common/index_by_chr.awk "${outputDir}"/merged30.txt 500000 > "${outputDir}"/merged30_index.txt
	  fi
    "${juiceDir}"/scripts/common/juicer_tools pre -n -s "$outputDir"/inter_30.txt -g "$outputDir"/inter_30_hists.m -q 30 "$resolutionsToBuildString" "$buildFragmentMapString" "$threadHic30String" "$outputDir"/merged30.txt "$outputDir"/inter_30.hic "$genomeID"
	  "${juiceDir}"/scripts/common/juicer_tools addNorm "$threadNormString" "${outputDir}"/inter_30.hic
	  rm -Rf "${tempdirPre30}"
fi
if [ -z $early ]
then
    # Create loop lists file for MQ > 30
    "${juiceDir}"/scripts/common/juicer_hiccups.sh -j "${juiceDir}"/scripts/common/juicer_tools -i "$outputDir"/inter_30.hic -m "${juiceDir}"/references/motif -g "$genomeID"
    "${juiceDir}"/scripts/common/juicer_arrowhead.sh -j "${juiceDir}"/scripts/common/juicer_tools -i "$outputDir"/inter_30.hic
fi

if [ -s "${outputDir}"/inter.hic ] && [ -s "${outputDir}"/inter_30.hic ]
then
   rm -fr "${tmpdir}"
   echo "(-: Successfully completed making mega map. Done. :-)"
else
   echo "!*** Error: one or both hic files are empty. Check debug directory for hic logs"
fi