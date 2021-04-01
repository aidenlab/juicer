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
# Single CPU version of Juicer. 
#
# Alignment script. Sets the reference genome and genome ID based on the input
# arguments (default human, none). Optional arguments are the description for stats file, 
# stage to relaunch at, paths to various files if needed,
# path to scripts directory, and the top-level directory (default 
# current directory). In lieu of setting the genome ID, you can instead set the
# reference sequence and the chrom.sizes file path, but the directory 
# containing the reference sequence must also contain the BWA index files.
#
# Aligns the fastq files, handles chimeric reads, sorts, and merges. 
#
# If all is successful, takes the final merged bam file, marks duplicates,
# creates hic contact maps, normalizes, and annotates features. 
# Final product will be hic file, stats file, dedup bam in the aligned directory.
#                                                                       
# [topDir]/fastq  - Should contain the fastq files. This code assumes that
#                   there is an "R" in the appropriate files, i.e. *R*.fastq
# From the top-level directory, the following two directories are created:
#                                                                              
# [topDir]/splits  - Where to write the scratch split files (fastq files and
#                    intermediate SAM files). This can be deleted after 
#                    execution.
# [topDir]/aligned - Where to write the final output files.
#
# The following globals should be set correctly before proceeding:
#
# splitsize - The number of lines that each split fastq should contain. Larger
#             means fewer files and longer overall, but too small means there
#             are so many jobs that the cluster won't run them. This can be
#             set with the -C command as well
# read1str  - portion of fastq filename that indicates this is the "read 1"
#             file; used to loop over only the read 1 and within that loop,
#             also align read 2 and merge.  If this is not set correctly,
#             script will not work. The error will often manifest itself
#             through a "*" in the name because the wildcard was not able to
#             match any files with the read1str.   
# Juicer version 2.0
shopt -s extglob
export LC_ALL=C
juicer_version="2.0"
### LOAD BWA AND SAMTOOLS 
bwa_cmd="bwa"
# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value 
read1str="_R1"
read2str="_R2"

## Default options, overridden by command line arguments 
# Juicer directory, contains scripts/, references/, and restriction_sites/ 
# can also be set in options via -D   
juiceDir="/aidenlab"

# size to split fastqs. adjust to match your needs. 4000000=1M reads per split
# can also be changed via the -C flag
splitsize=90000000

# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
# default: not set
site="none"
# genome ID, default to human, can also be set in options
genomeID="hg19"
# description, default empty
about=""
# do not include fragment delimited maps by default
nofrag=1
# use wobble for dedupping by default (not just exact matches)
justexact=0

## TODO Change usage to be correct and print nicely 
## Read arguments                                                     
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-s site]\n                 [-a about] [-S stage] [-p chrom.sizes path]\n                 [-y restriction site file] [-z reference genome file]\n                 [-D Juicer scripts directory]\n                 [-b ligation] [-t threads] [-T threadsHic]\n                 [-e] [-h] [-f] [-j]"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \n  \"$genomeID\"); alternatively, it can be defined using the -z command"
dirHelp="* [topDir] is the top level directory (default\n  \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
siteHelp="* [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" \n  (default \"$site\")"
aboutHelp="* [about]: enter description of experiment, enclosed in single quotes"
stageHelp="* [stage]: must be one of \"chimeric\", \"merge\", \"dedup\", \"final\", \"postproc\", or \"early\".\n    -Use \"merge\" when alignment has finished but the merged_sort file has not\n     yet been created.\n    -Use \"dedup\" when the files have been merged into merged_sort but\n     merged_nodups has not yet been created.\n    -Use \"final\" when the reads have been deduped into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the hic files\n    Can also use -e flag to exit early"
pathHelp="* [chrom.sizes path]: enter path for chrom.sizes file"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
refSeqHelp="* [reference genome file]: enter path for reference sequence file, BWA index\n  files must be in same directory"
ligationHelp="* [ligation junction]: use this string when counting ligation junctions"
threadsHelp="* [threads]: number of threads when running BWA alignment"
threadsHicHelp="* [threads for hic file creation]: number of threads when building hic file"
excludeHelp="* -f: include fragment-delimited maps in hic file creation"
justHelp="* -j: just exact duplicates excluded at dedupping step"
earlyexitHelp="* -e: Use for an early exit, before the final creation of the hic files"
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$genomeHelp"
    echo -e "$dirHelp"
    echo -e "$siteHelp"
    echo -e "$aboutHelp"
    echo -e "$stageHelp"
    echo -e "$pathHelp"
    echo -e "$siteFileHelp"
    echo -e "$refSeqHelp"
    echo -e "$scriptDirHelp"
    echo -e "$ligationHelp"
    echo -e "$threadsHelp"
    echo -e "$threadsHicHelp"
    echo -e "$userHelp"
    echo -e "$earlyexitHelp"
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:a:hs:p:y:z:S:D:b:t:jfecT:1:2:i:" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	s) site=$OPTARG ;;
	a) about=$OPTARG ;;
	p) genomePath=$OPTARG ;;  
	y) site_file=$OPTARG ;;
	z) refSeq=$OPTARG ;;
	S) stage=$OPTARG ;;
	D) juiceDir=$OPTARG ;;
	f) nofrag=0 ;;
	b) ligation=$OPTARG ;;
	t) threads=$OPTARG ;;
	1) read1files=$OPTARG ;;
	2) read2files=$OPTARG ;;
	j) justexact=1 ;;
	e) earlyexit=1 ;;
	T) threadsHic=$OPTARG ;;
	i) sampleName=$OPTARG ;;
	[?]) printHelpAndExit 1;;
    esac
done

if [ ! -z "$stage" ]
then
    case $stage in
	chimeric) chimeric=1 ;;
        merge) merge=1 ;;
        dedup) dedup=1 ;;
        early) earlyexit=1 ;;
        final) final=1 ;;
	postproc) postproc=1 ;; 
	alignonly) alignonly=1 ;;
	chimericonly) chimericonly=1 ;;
	deduponly) deduponly=1 ;;
        *)  echo "$usageHelp"
	    echo "$stageHelp"
	    exit 1
    esac
fi

## Set reference sequence based on genome ID
if [ -z "$refSeq" ]
then 
    case $genomeID in
	mm9)	refSeq="${juiceDir}/references/Mus_musculus_assembly9_norandom.fasta";;
	mm10)	refSeq="${juiceDir}/references/Mus_musculus_assembly10/v0/Mus_musculus_assembly10.fasta";;
	hg38)	refSeq="${juiceDir}/references/hg38/hg38.fa";;
	hg19)	refSeq="${juiceDir}/references/Homo_sapiens_assembly19.fasta";;
	hg18)	refSeq="${juiceDir}/references/hg18.fasta";;
	*)	echo "$usageHelp"
	    echo "$genomeHelp"
	    exit 1
    esac
else
    ## Reference sequence passed in, so genomePath must be set for the .hic 
    ## file to be properly created
    if [[ -z "$genomePath" ]] && [[ -z $earlyexit ]] && [ -z "$alignonly" ]
    then
        echo "***! You must define a chrom.sizes file or a standard genome ID via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq; you may use \"-p hg19\" or other standard genomes";
        exit 1;
    fi
fi

## Check that refSeq exists 
if [ ! -e "$refSeq" ]; then
    echo "***! Reference sequence $refSeq does not exist";
    exit 1;
fi

## Check that index for refSeq exists
if [[ ! -e "${refSeq}.bwt" ]] 
then
    echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run bwa index on this file before running juicer.";
    exit 1;
fi

## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]; then
    case $site in
	HindIII) ligation="AAGCTAGCTT";;
	MseI)  ligation="TTATAA";;
	DpnII) ligation="GATCGATC";;
	MboI) ligation="GATCGATC";;
        NcoI) ligation="CCATGCATGG";;
	none) ligation="XXXX";;
	Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
	*)  ligation="XXXX"
	    echo "$site not listed as recognized enzyme."
	    echo "Ligation junction is undefined"
    esac
fi

## If DNAse-type experiment, no fragment maps; or way to get around site file
if [[ "$site" == "none" ]] 
then
    nofrag=1;
fi

if [ -z "$site_file" ]
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [[ ! -e "$site_file" ]] && [[ "$site" != "none" ]] &&  [[ ! "$site_file" =~ "none" ]]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    exit 1
elif [[ "$site" != "none" ]] && [[ ! "$site_file" =~ "none" ]]
then
    echo  "Using $site_file as site file"
fi

## Set threads for sending appropriate parameters to cluster and string for BWA call
if [ -z "$threads" ]
then
    threads="$(getconf _NPROCESSORS_ONLN)"
    threadstring="-t $threads"
    sthreadstring="-@ $threads"
else
    threadstring="-t $threads"
    sthreadstring="-@ $threads"
fi

if [ -n "$read2files" ] && [ -z "$read1files" ]
then
    echo "***! When fastqs for read2 are specified with \"-2\", corresponding read1 fastqs must be specified with \"-1\" "
    exit 1
fi

### If not set in options, set samplename based on directory name
if [ -z "$sampleName" ]
then
    sampleName="NULL"
fi

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/aligned"

if [ -z "$threadsHic" ]
then
    threadsHic=1
    threadHicString=""
    threadHic30String=""
    threadNormString=""
else
    threadHicString="--threads $threadsHic -i ${outputdir}/merged0_index.txt -t ${outputdir}/HIC_tmp"
    threadHic30String="--threads $threadsHic -i ${outputdir}/merged30_index.txt -t ${outputdir}/HIC30_tmp"
    threadNormString="--threads $threadsHic"
fi

## Create output directory, only if not in postproc, dedup or final stages
if [[ -d "$outputdir" && -z "$final" && -z "$dedup" && -z "$postproc" && -z "$deduponly" ]]
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 1			
else
    if [[ -z "$final" && -z "$dedup" && -z "$postproc" && -z "$deduponly" ]]; then
        mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 1; } 
    fi
fi

## Create split directory
if [ -d "$splitdir" ]; then
    splitdirexists=1
else
    mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 1; }
fi

if [ -z "$read1files" ]
then
## Check that fastq directory exists and has proper fastq files; only if necessary
    if [[ -z "$final" && -z "$dedup" && -z "$postproc" && -z "$deduponly" && -z "$merge" && -z "$mergeonly" ]]; then
	if [ ! -d "$topDir/fastq" ]; then
	    echo "Directory \"$topDir/fastq\" does not exist."
	    echo "Create \"$topDir/fastq\" and put fastq files to be aligned there."
	    echo "Type \"juicer.sh -h\" for help"
	    exit 1
	else
	    if stat -t ${fastqdir} >/dev/null 2>&1
	    then
		echo "(-: Looking for fastq files...fastq files exist"
		if [ ! $splitdirexists ]
		then
		    echo "(-: Created $splitdir."
		    ln -s ${fastqdir} ${splitdir}/.
		else
		    echo -e "---  Using already created files in $splitdir\n"
		fi
		testname=$(ls -lgG ${fastqdir} | awk 'NR==1{print $7}')
		if [ "${testname: -3}" == ".gz" ]
		then
		    read1=${splitdir}"/*${read1str}*.fastq.gz"
		    gzipped=1
		else
		    read1=${splitdir}"/*${read1str}*.fastq"
		fi
	    else
		if [ ! -d "$splitdir" ]; then
		    echo "***! Failed to find any files matching ${fastqdir}"
		    echo "***! Type \"juicer.sh -h \" for help"
		    exit
		fi
	    fi
	fi
    fi
    read1files=()
    read2files=()
    for i in ${read1}
    do
	ext=${i#*$read1str}
	name=${i%$read1str*}
        # these names have to be right or it'll break   
	name1=${name}${read1str}
	name2=${name}${read2str}
	read1filesstr+=$name1$ext" "
	read2filesstr+=$name2$ext" "
    done
    read1files=( $read1filesstr )
    read2files=( $read2filesstr )
else
    if [ -z "$read2files" ]
    then
	echo "***! When fastqs for read1 are specified with \"-1\", corresponding read2 fastqs must be specified with \"-2\" "
	exit 1
    else
	 # replace commas with spaces for iteration, put in array 
	read1files=($(echo $read1files | tr ',' ' '))
	read2files=($(echo $read2files | tr ',' ' '))
    fi
fi

if [ "${#read1files[@]}" -ne "${#read2files[@]}" ]
then
    echo "***! The number of read1 fastqs specified (${#read1files[@]}) is not equal to the number of read2 fastqs specified (${#read2files[@]})"
    exit 1
fi

## Arguments have been checked and directories created. Now begins 
## the real work of the pipeline
headfile=${outputdir}/header
date > $headfile
# Experiment description
if [ -n "${about}" ]
then
    echo -ne 'Experiment description: ${about}; ' >> $headfile
else
    echo -ne 'Experiment description: ' >> $headfile
fi
echo -ne "Sample name $sampleName;"  >> $headfile
# Get version numbers of all software   
echo -ne " Juicer version $juicer_version;" >> $headfile
$bwa_cmd 2>&1 | awk '$1=="Version:"{printf(" BWA %s; ", $2)}' >> $headfile
echo -ne "$threads threads; " >> $headfile
java -version 2>&1 | awk 'NR==1{printf("%s; ", $0);}' >> $headfile
${juiceDir}/scripts/common/juicer_tools -V 2>&1 | awk '$1=="Juicer" && $2=="Tools"{printf("%s; ", $0);}' >> $headfile
echo "$0 $@" >> $headfile


## ALIGN FASTQ IN SINGLE END MODE, SORT BY READNAME, HANDLE CHIMERIC READS     

## Not in merge, dedup, final, or postproc stage, i.e. need to align files. 
if [ -z $merge ] && [ -z $mergeonly ] && [ -z $final ] && [ -z $dedup ] && [ -z $deduponly ] && [ -z $postproc ]
then
    if [ "$nofrag" -eq 0 ]
    then
	echo -e "(-: Aligning files matching $fastqdir\n to genome $genomeID with site file $site_file"
    else
        echo -e "(-: Aligning files matching $fastqdir\n to genome $genomeID with no fragment delimited maps."
    fi

    for ((i = 0; i < ${#read1files[@]}; ++i)); do
	usegzip=0
	file1=${read1files[$i]}
	file2=${read2files[$i]}

	ext=${file1#*$read1str}
	name=${file1%$read1str*} 

	# these names have to be right or it'll break
	name1=${name}${read1str}
	name2=${name}${read2str}	
	jname=$(basename "$name")${ext}
	rgtag="@RG\tID:${sampleName}\tLB:NULL\tPL:ILLUMINA\tSM:NULL"
	
	if [ ${file1: -3} == ".gz" ]
	then
	    usegzip=1
	fi

	source ${juiceDir}/scripts/common/countligations.sh
	if [ -z "$chimeric" ]
	then
	    echo "Running command $bwa_cmd mem -SP5M $threadstring $refSeq $file1 $file2 > $name$ext.sam"
	    $bwa_cmd mem -SP5M -R "${rgtag}" $threadstring $refSeq $file1 $file2 > $name$ext.sam
	    if [ $? -ne 0 ]
	    then
		echo "***! Alignment of $file1 $file2 failed."
		exit 1
	    else
		echo "(-:  Align of $name$ext.sam done successfully"
	    fi
	else
	    echo "Using already aligned reads $name$ext.sam";
        fi

	# call chimeric script to deal with chimeric reads; sorted file is sorted by read name at this point		
	if [ "$site" != "none" ] && [ -e "$site_file" ] 
	then		
	    awk -v stem=${name}${ext}_norm -v site_file=$site_file -f $juiceDir/scripts/common/chimeric_sam.awk $name$ext.sam | samtools sort -t cb -n $sthreadstring >  ${name}${ext}.bam
	else
	    awk -v stem=${name}${ext}_norm -f $juiceDir/scripts/common/chimeric_sam.awk $name$ext.sam > $name$ext.sam2
	    awk -v avgInsertFile=${name}${ext}_norm.txt.res.txt -f $juiceDir/scripts/common/adjust_insert_size.awk $name$ext.sam2 | samtools sort -t cb -n $sthreadstring >  ${name}${ext}.bam
	fi

	if [ $? -ne 0 ]
	then
	    echo "***! Failure during chimera handling of $name${ext}"
	    exit 1 
	else  
	    rm $name$ext.sam* 
	fi  
    done # done looping over all fastq split files
fi  # Not in merge, dedup,  or final stage, i.e. need to split and align files.

if [ -n "$alignonly" ]
then
    exit 0
fi

#MERGE SORTED AND ALIGNED FILES
# Not in final, dedup, or postproc
if [ -z $final ] && [ -z $dedup ] && [ -z $deduponly ] && [ -z $postproc ]
then
    if [ -d $donesplitdir ]
    then
	mv $donesplitdir/* $splitdir/.
    fi
    
    if ! samtools merge -t cb -n $sthreadstring $outputdir/merged_sort.bam  $splitdir/*.bam
    then
	echo "***! Some problems occurred somewhere in creating sorted align files."
	exit 1
    else
	echo "(-: Finished sorting all sorted files into a single merge."
    fi
fi

if [ -n "$mergeonly" ]
then
    exit 0
fi

#REMOVE DUPLICATES 
if [ -z $final ] && [ -z $postproc ]
then
    if [ $justexact -eq 1 ]
    then
	samtools view $sthreadstring -h $outputdir/merged_sort.bam | awk -f $juiceDir/scripts/common/dups_sam.awk -v nowobble=1 -v fname=$outputdir/dup_count > $outputdir/merged_dedup.sam
    else
	samtools view $sthreadstring -h $outputdir/merged_sort.bam | awk -f $juiceDir/scripts/common/dups_sam.awk -v fname=$outputdir/dup_count > $outputdir/merged_dedup.sam
    fi
    if [ $? -ne 0 ]
    then
	echo "***! Mark duplicates of $outputdir/merged_sort.bam failed."
	exit 1
    else
	echo "(-:  Mark duplicates done successfully"

    fi
fi

if [ -n "$deduponly" ]
then
    exit 0
fi

#CREATE HIC FILES 
if [ -z "$genomePath" ]
then
    #If no path to genome is give, use genome ID as default.
    genomePath=$genomeID
fi

#Skip if post-processing only is required
if [ -z $postproc ]
    then

    # Check that dedupping worked properly
    # in ideal world, we would check this in split_rmdups and not remove before we know they are correct
    size1=$(samtools view -h ${outputdir}/merged_sort.bam | wc -l | awk '{print $1}')
    size2=$(wc -l ${outputdir}/merged_dedup.sam | awk '{print $1}')
	
    if [ $size1 -ne $size2 ]
    then
	echo "***! Error! The sorted file and dups/no dups files do not add up, or were empty."
	exit 1
    fi

    samtools view $sthreadstring -F 1024 -O sam ${outputdir}/merged_dedup.sam | awk -v mapq=1 -f ${juiceDir}/scripts/common/sam_to_pre.awk > ${outputdir}/merged0.txt
    samtools view $sthreadstring -F 1024 -O sam ${outputdir}/merged_dedup.sam | awk -v mapq=30 -f ${juiceDir}/scripts/common/sam_to_pre.awk > ${outputdir}/merged30.txt

    if [[ $threadsHic -gt 1 ]]
    then
	time ${juiceDir}/scripts/common/index_by_chr.awk ${outputdir}/merged0.txt 500000 > ${outputdir}/merged0_index.txt
	time ${juiceDir}/scripts/common/index_by_chr.awk ${outputdir}/merged30.txt 500000 > ${outputdir}/merged30_index.txt
    fi

    if samtools view -b $sthreadstring ${outputdir}/merged_dedup.sam > ${outputdir}/merged_dedup.bam
    then
	rm ${outputdir}/merged_dedup.sam
	rm ${outputdir}/merged_sort.bam
    fi

    export IBM_JAVA_OPTIONS="-Xmx1024m -Xgcthreads1"
    export _JAVA_OPTIONS="-Xmx1024m -Xms1024m"
    tail -n1 $headfile | awk '{printf"%-1000s\n", $0}' > $outputdir/inter.txt
    cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/common/stats_sub.awk >> $outputdir/inter.txt
    awk 'NR==FNR{if ($1 ~ /Alignable/){split($0,a,":"); split(a[2],b); gsub(",","", b[1]); tot=int(b[1]);}}FNR!=NR{sum+=$1}END{print "Unique Reads: ",int(tot-sum); print "Duplicates:",sum}' $outputdir/inter.txt $outputdir/*_count >> $outputdir/inter.txt
    cp $outputdir/inter.txt $outputdir/inter_30.txt
    perl ${juiceDir}/scripts/common/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt $outputdir/merged0.txt
    perl ${juiceDir}/scripts/common/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt $outputdir/merged30.txt

    # if early exit, we stop here, once the stats are calculated
    if [ ! -z "$earlyexit" ]
    then
	export splitdir=${splitdir}; export outputdir=${outputdir}; export early=1; ${juiceDir}/scripts/common/check.sh
    fi
    export IBM_JAVA_OPTIONS="-Xmx150000m -Xgcthreads1"
    export _JAVA_OPTIONS="-Xmx150000m -Xms150000m"
    mkdir ${outputdir}"/HIC_tmp"
    if [ "$nofrag" -eq 1 ]
    then 
	time ${juiceDir}/scripts/common/juicer_tools pre -n -s $outputdir/inter.txt -g $outputdir/inter_hists.m -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 $threadHicString $outputdir/merged0.txt $outputdir/inter.hic $genomePath
    else
	time ${juiceDir}/scripts/common/juicer_tools pre -n -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 $threadHicString $outputdir/merged0.txt $outputdir/inter.hic $genomePath 
    fi
    if [[ $threadsHic -gt 1 ]]
    then 
	time ${juiceDir}/scripts/common/juicer_tools addNorm $threadNormString ${outputdir}/inter.hic
    fi
    rm -R ${outputdir}"/HIC_tmp"
    date
    mkdir ${outputdir}"/HIC30_tmp"
    if [ "$nofrag" -eq 1 ]
    then 
	time ${juiceDir}/scripts/common/juicer_tools pre -n -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 $threadHic30String $outputdir/merged30.txt $outputdir/inter_30.hic $genomePath
    else
	time ${juiceDir}/scripts/common/juicer_tools pre -n -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100 $threadHic30String $outputdir/merged30.txt $outputdir/inter_30.hic $genomePath
    fi
    if [[ $threadsHic -gt 1 ]]
    then 
	time ${juiceDir}/scripts/common/juicer_tools addNorm $threadNormString ${outputdir}/inter_30.hic
    fi
    rm -R ${outputdir}"/HIC30_tmp"
    # POSTPROCESSING 
    ${juiceDir}/scripts/common/juicer_postprocessing.sh -j ${juiceDir}/scripts/common/juicer_tools -i ${outputdir}/inter_30.hic -m ${juiceDir}/references/motif -g ${genomeID}  
fi
#CHECK THAT PIPELINE WAS SUCCESSFUL
export early=$earlyexit
export splitdir=$splitdir
source ${juiceDir}/scripts/common/check.sh
