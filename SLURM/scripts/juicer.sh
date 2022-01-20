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
# Alignment script. Sets the reference genome and genome ID based on the input
# arguments (default human, none). Optional arguments are the queue for the 
# alignment, description for stats file, 
# stage to relaunch at, paths to various files if needed,
# chunk size, path to scripts directory, and the top-level directory (default 
# current directory). In lieu of setting the genome ID, you can instead set the
# reference sequence and the chrom.sizes file path, but the directory 
# containing the reference sequence must also contain the BWA index files.
#
# Splits the fastq files, creates jobs to align them, creates merge jobs that
# wait for the alignment to finish, and creates a final merge job.
#
# Also creates "cleanup" jobs that at each stage, deletes jobs off the cluster
# if any one of them fails.
#
# If all is successful, takes the final merged file, removes name duplicates,
# removes PCR duplicates, and creates the hic job and stats job.  Final
# product will be hic file and stats file in the aligned directory.
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
juicer_version="2.0"
## Set the following variables to work with your system

# Aiden Lab specific check
isRice=$(host $(hostname) | awk 'NR==1{if ($1~/rice/) {print 1} else {print 0}; exit}') #'
isBCM=$(host $(hostname) | awk 'NR==1{if ($1~/bcm/) {print 1} else {print 0}; exit}') #'
isVoltron=0
## path additionals, make sure paths are correct for your system
## use cluster load commands
if [ $isRice -eq 1 ] 
then
    export PATH=$HOME/bin:$PATH
    isNots=$(host $(hostname) | awk 'NR==1{if ($1~/nots/){print 1}else {print 0}}') #'
    if [ $isNots -eq 1 ]
    then
	load_bwa="module load  GCCcore/7.3.0 BWA/0.7.17"
	load_java="module load Java/1.8.0_162" 
	load_gpu="module load gcccuda/2016a;module load CUDA/8.0.44;" 
	load_samtools=""
    else
	load_bwa="export PATH=/home/ncd4/bwa:$PATH"
	load_java="module load Java/8.0.3.22" 
	load_gpu="module load gcccuda/2016a;module load CUDA/8.0.54;" 
	load_samtools=""
    fi
    
    # Juicer directory, contains scripts/, references/, and restriction_sites/
    # can also be set in options via -D
    juiceDir="/projects/ea14/juicer" ### RICE
    # default queue, can also be set in options via -q
    queue="commons"
    queue_time="24:00:00"
    # default long queue, can also be set in options via -l
    long_queue="commons"
    long_queue_time="24:00:00"
elif [ $isBCM -eq 1 ]
then    
    # Juicer directory, contains scripts/, references/, and restriction_sites/
    # can also be set in options via -D
    # XXX verify / change me
    juiceDir="/gpfs0/juicer2"
    load_bwa="unset PYTHONPATH; conda activate"
    # default queue, can also be set in options via -q
    queue="mhgcp"
    queue_time="1200"
    # default long queue, can also be set in options via -l
    long_queue="mhgcp"
    long_queue_time="3600"
else
    isVoltron=1
    #export PATH=/gpfs0/biobuild/biobuilds-2016.11/bin:$PATH 
    # unset MALLOC_ARENA_MAX # on IBM platform this parameter has significant speed efect but may result in memory leaks
    load_bwa="spack load bwa@0.7.17 arch=\`spack arch\`"
    load_awk="spack load gawk@4.1.4 arch=\`spack arch\`"
    load_gpu="spack load cuda@8.0.61 arch=\`spack arch\` && CUDA_VISIBLE_DEVICES=0,1,2,3"
    load_samtools="spack load samtools@1.13 arch=\`spack arch\`"
    call_bwameth="/gpfs0/home/neva/bwa-meth/bwameth.py"
    # Juicer directory, contains scripts/, references/, and restriction_sites/
    # can also be set in options via -D
    juiceDir="/gpfs0/juicer2"
    # default queue, can also be set in options
    queue="commons"
    queue_time="2880"
    # default long queue, can also be set in options
    long_queue="long"
    long_queue_time="7200"
fi

# size to split fastqs. adjust to match your needs. 4000000=1M reads per split
# can also be changed via the -C flag
splitsize=90000000

# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1" 
read2str="_R2" 
ostem="output" # potentially make a flag for this 

# unique name for jobs in this run
groupname="a$(date +%s)"

## Default options, overridden by command line arguments

# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
# default: not set
site="none"
# description, default empty
about=""
# do not include fragment delimited maps by default
nofrag=1
# use wobble for dedupping by default (not just exact matches)
justexact=0
wobbleDist=4
# assembly mode, produce old merged_nodups, early exit
assembly=0
# force cleanup
cleanup=0
# qc apa
qc_apa=0
# single-end input, default no
singleend=0
# sample name for RG tag
sampleName="HiC_sample"
# library name for RG tag
libraryName="HiC_library"

## Read arguments                                                     
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]\n                 [-a about] [-S stage] [-p chrom.sizes path] [-C chunk size]\n                 [-y restriction site file] [-z reference genome file]\n                 [-D Juicer scripts parent dir] [-Q queue time limit]\n                 [-L long queue time limit] [-b ligation] [-t threads]\n                 [-T threadsHic] [-A account] [-i sample] [-k library] [-w wobble]\n                 [-e] [-h] [-f] [-j] [-u] [-m] [--assembly] [--cleanup] [--qc_apa] [--qc] [--in-situ]"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\"; alternatively, it can be defined using the -z command"
dirHelp="* [topDir] is the top level directory (default\n  \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
queueHelp="* [queue] is the queue for running alignments (default \"$queue\")"
longQueueHelp="* [long queue] is the queue for running longer jobs such as the hic file\n  creation (default \"$long_queue\")"
siteHelp="* [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" \n  (default \"$site\")"
aboutHelp="* [about]: enter description of experiment, enclosed in single quotes"
stageHelp="* [stage]: must be one of \"chimeric\", \"merge\", \"dedup\", \"afterdedup\", \"final\", \"postproc\", or \"early\".\n    -Use \"chimeric\" when alignment has finished or to start from previously\n     aligned files\n    -Use \"merge\" when chimeric handling has finished but the merged_sort file\n     has not yet been created.\n    -Use \"dedup\" when the files have been merged into merged_sort but\n     merged_dedup has not yet been created.\n    -Use \"afterdedup\" when dedup is complete but statistics haven't been run\n    -Use \"final\" when the reads have been deduped into merged_dedup but the\n     final hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the hic files\n    Can also use -e flag to exit early"
pathHelp="* [chrom.sizes path]: enter path for chrom.sizes file; can also use canonical\n  genome name here such as hg38"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
chunkHelp="* [chunk size]: number of lines in split files, must be multiple of 4\n  (default ${splitsize}, which equals $(awk -v ss=${splitsize} 'BEGIN{print ss/4000000}') million reads)"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
refSeqHelp="* [reference genome file]: enter path for reference sequence file, BWA index\n  files must be in same directory"
queueTimeHelp="* [queue time limit]: time limit for queue, i.e. -W 12:00 is 12 hours\n  (default ${queue_time})"
longQueueTimeHelp="* [long queue time limit]: time limit for long queue, i.e. -W 168:00 is one week\n  (default ${long_queue_time})"
ligationHelp="* [ligation junction]: use this string when counting ligation junctions"
threadsHelp="* [threads]: number of threads when running BWA alignment"
threadsHicHelp="* [threads for hic file creation]: number of threads when building hic file"
userHelp="* [account name]: user account name on cluster"
sampleHelp="* [sample name]: will be added to the SM portion of the read group (RG) tag"
libraryHelp="* [library name]: will be added to the LB portion of the read group (RG) tag"
wobbleHelp="* [wobble dist]: adjust wobble for deduping (default 4)"
excludeHelp="* -f: include fragment-delimited maps in hic file creation"
justHelp="* -j: just exact duplicates excluded at dedupping step"
earlyexitHelp="* -e: Use for an early exit, before the final creation of the hic files"
singleEndHelp="* -u: Single end alignment"
methylationHelp="* -m: Methylation + Hi-C library"
assemblyHelp="* --assembly: For use before 3D-DNA; early exit and create old style merged_nodups"
cleanupHelp="* --cleanup: Automatically clean up files if pipeline successfully completes"
qcapaHelp="* --qc_apa: Run QC APA"
qcHelp="* --qc: Only build map down to 1000bp, no feature annotation" 
insituHelp="* --in-situ: Only build map down to 1000bp"
helpHelp="* -h, --help: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$genomeHelp"
    echo -e "$dirHelp"
    echo -e "$queueHelp"
    echo -e "$longQueueHelp"
    echo -e "$siteHelp"
    echo -e "$aboutHelp"
    echo -e "$stageHelp"
    echo -e "$pathHelp"
    echo -e "$siteFileHelp"
    echo -e "$refSeqHelp"
    echo -e "$chunkHelp"
    echo -e "$scriptDirHelp"
    echo -e "$queueTimeHelp"
    echo -e "$longQueueTimeHelp"
    echo -e "$ligationHelp"
    echo -e "$threadsHelp"
    echo -e "$threadsHicHelp"
    echo -e "$userHelp"
    echo -e "$sampleHelp"
    echo -e "$libraryHelp"
    echo -e "$wobbleHelp"
    echo -e "$justHelp"
    echo -e "$earlyexitHelp"
    echo -e "$excludeHelp"
    echo -e "$singleEndHelp"
    echo -e "$methylationHelp"
    echo -e "$assemblyHelp"
    echo -e "$cleanupHelp"
    echo -e "$qcapaHelp"
    echo -e "$qcHelp"
    echo -e "$insituHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:a:hq:s:p:l:y:z:S:C:D:Q:L:b:A:i:t:jfuec-:T:w:k:m" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	l) long_queue=$OPTARG ;;
	q) queue=$OPTARG ;;
	s) site=$OPTARG ;;
	a) about=$OPTARG ;;
	p) genomePath=$OPTARG ;;  
	y) site_file=$OPTARG ;;
	z) refSeq=$OPTARG ;;
	S) stage=$OPTARG ;;
	C) splitsize=$OPTARG; splitme=1 ;;
	D) juiceDir=$OPTARG ;;
	Q) queue_time=$OPTARG ;;
	L) long_queue_time=$OPTARG ;;
	f) nofrag=0 ;;
	b) ligation=$OPTARG ;;
	t) threads=$OPTARG ;;
	A) user=$OPTARG ;;
	j) justexact=1 ;;
	e) earlyexit=1 ;;
	T) threadsHic=$OPTARG ;;
	i) sampleName=$OPTARG ;;
	u) singleend=1 ;;
	w) wobbleDist=$OPTARG ;;
	k) libraryName=$OPTARG ;;
	m) methylation=1 ;;
	-) case "${OPTARG}" in 
	    assembly) earlyexit=1; assembly=1 ;;
	    cleanup)  cleanup=1 ;;
	    qc_apa)   qc_apa=1 ;;
	    qc) qc=1 ;;
	    "help")   printHelpAndExit 0;;
	    in-situ) insitu=1 ;;
	    *) echo "Unknown argument --${OPTARG}"; 
	       printHelpAndExit 1;;
           esac;;
	[?]) printHelpAndExit 1;;
    esac
done

if [ ! -z "$stage" ]
then
    case $stage in
	chimeric) chimeric=1 ;;
        merge) merge=1 ;;
        dedup) dedup=1 ;;
	afterdedup) afterdedup=1 ;;
        early) earlyexit=1 ;;
        final) final=1 ;;
	postproc) postproc=1 ;; 
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
	GRCh38) 
	    refSeq="${juiceDir}/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
	    site_file="${juiceDir}/restriction_sites/ENCFF132WAM.txt"
	    genomeID="hg38"
	    ;;
	hg19)	refSeq="${juiceDir}/references/Homo_sapiens_assembly19.fasta";;
	hg18)	refSeq="${juiceDir}/references/hg18.fasta";;
	*)	echo "$usageHelp"
	    echo "$genomeHelp"
	    exit 1
    esac
else
    ## Reference sequence passed in, so genomePath must be set for the .hic 
    ## file to be properly created
    if [[ -z "$genomePath" ]] && [[ -z $earlyexit ]]
    then
        echo "***! You must define a chrom.sizes file or a standard genome ID via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq; you may use \"-p hg19\" or other standard genomes";
        exit 1;
    fi
fi

## Alignment checks; not necessary if later stages
if [[ -z "$chimeric" && -z "$merge" &&  -z "$final" && -z "$dedup" && -z "$postproc" && -z "$afterdedup" ]] 
then
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
fi

## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]; then
    case $site in
	HindIII) ligation="AAGCTAGCTT";;
	MseI)  ligation="TTATAA";;
	DpnII) ligation="GATCGATC";;
	MboI) ligation="GATCGATC";;
        MboI) ligation="GATCGATC";;
	NlaIII) ligation="CATG";;
        NcoI) ligation="CCATGCATGG";;
	MspI) ligation="CCGCGG";;
	HinP1I) ligation="GCGCGC";;
	StyD4I) ligation="CCNGGCCNGG";;
	SaII) ligation="GTCGATCGAC";;
	NheI) ligation="GCTAGCTAGC";;
	StyI) ligation="CCWWGCWWGG";;
	XhoI) ligation="CTCGATCGAG";;
	BglII) ligation="AGATCGATCT";;
	CviJI) ligation="'(AGCC|GGCT|AGCT|GGCC|GGGG|GGGA|AGGG|CCCT|CCCC|TCCC|AGAG|CTCT)'";;
	MboI+MseI) ligation="'(GATCGATC|TTATAA)'";;
	CviQI+MseI) ligation="'(GTATAC|TTATAA|GTATAA|TTATAC)'";;
	MseI+CviAII) ligation="'(CATATG|CATTAA|TTAATG|TTATAA)'";;
	Csp6I+MseI) ligation="'(GTATAC|TTATAA|GTATAA|TTATAC)'";;
	MseI+Csp6I) ligation="'(GTATAC|TTATAA|GTATAA|TTATAC)'";;
	MluCI) ligation="AATTAATT";;
	CviAII) ligation="CATATG";;
        Csp6I) ligation="GTATAC";;
        Csp6I+MseI) ligation="'(GTATAC|TTATAA|GTATAA|TTATAC)'";;
        MseI+Csp6I) ligation="'(GTATAC|TTATAA|GTATAA|TTATAC)'";;
        MseI+CviAII) ligation="'(CATATG|CATTAA|TTAATG|TTATAA)'";;
        CviAII+MseI) ligation="'(CATATG|CATTAA|TTAATG|TTATAA)'";;
        CviAII+Csp6I) ligation="'(CATATG|GTATAC|CATTAC|GTAATG)'";;
        Csp6I+CviAII) ligation="'(CATATG|GTATAC|CATTAC|GTAATG)'";;
        Csp6I+CviAII+MseI) ligation="'(CATATG|CATTAA|CATTAC|TTAATG|TTATAA|TTATAC|GTATAC|GTATAA|GTAATG)'";;
	Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
	none) ligation="XXXX";;
	*)  ligation="XXXX"
	    echo "$site not listed as recognized enzyme."
	    echo "Ligation junction is undefined"
    esac
fi

if [ "$methylation" = 1 ]
then
    ligation=$(echo $ligation | awk '{printf "'\''%s'\'' ", gensub("C","[CT]",$0)}')
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
    # default is 8 threads; may need to adjust
    if [ $isRice -eq 1 ]
    then
	threads=8
	sortthreads=8
	threadstring="-t $threads"
	sthreadstring="-@ $threads"
    elif [ $isVoltron -eq 1 ]
    then
	threads=8 
	## On voltron with 8 thread per core Power8 CPU bwa can use more threads
	threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"
	sthreadstring="-@ \$SLURM_JOB_CPUS_PER_NODE"
	sortthreads=8
    else
	# only one thread if undefined; also if isBCM
	threads=1
	sortthreads=1
	threadstring="-t $threads"
	sthreadstring="-@ $threads"

    fi
else
    if [ $isVoltron -eq 1 ]
    then
	threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"
	sthreadstring="-@ \$SLURM_JOB_CPUS_PER_NODE"
	sortthreads=8
    else
	threadstring="-t $threads"
	sthreadstring="-@ $threads"
	sortthreads=$threads
    fi
fi

alloc_mem=$(($threads * 8000))

if [ $alloc_mem -gt 80000 ]
then
    alloc_mem=80000
fi

if [ $isBCM -eq 1 ] || [ $isRice -eq 1 ]
then
    alloc_mem=50000
fi

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/HIC_tmp"
debugdir=${topDir}"/debug"

if [ -z "$threadsHic" ]
then
    threadsHic=1
    threadHicString=""
    threadHic30String=""
    threadNormString=""
else
    threadHicString="--threads $threadsHic -i ${outputdir}/merged1_index.txt -t ${outputdir}/HIC_tmp"
    threadHic30String="--threads $threadsHic -i ${outputdir}/merged30_index.txt -t ${outputdir}/HIC30_tmp"
    threadNormString="--threads $threadsHic"
fi

## Alignment checks; not necessary if later stages
if [[ -z "$chimeric" && -z "$merge" &&  -z "$final" && -z "$dedup" && -z "$afterdedup" && -z "$postproc" ]] 
then
    ## Check that fastq directory exists and has proper fastq files
    if [ ! -d "$topDir/fastq" ]; then
	echo "Directory \"$topDir/fastq\" does not exist."
	echo "Create \"$topDir/fastq\" and put fastq files to be aligned there."
	echo "Type \"juicer.sh -h\" for help"
	exit 1
    else 
	if stat -t ${fastqdir} >/dev/null 2>&1
	then
	    echo "(-: Looking for fastq files...fastq files exist"
	else
	    if [ ! -d "$splitdir" ]; then 
		echo "***! Failed to find any files matching ${fastqdir}"
		echo "***! Type \"juicer.sh -h \" for help"
		exit 1		
	    fi
	fi
    fi

    testname=$(ls -lgG ${fastqdir} | awk 'NR==1{print $7}')
    if [ "${testname: -3}" == ".gz" ]
    then
	read1=${splitdir}"/*${read1str}*.fastq.gz"
	gzipped=1
    else
	read1=${splitdir}"/*${read1str}*.fastq"
    fi
elif [[ -n "$chimeric" ]]
then
    read1=${splitdir}"/*.sam"
fi

## Create split directory
if [ -d "$splitdir" ]; then
    splitdirexists=1
elif  [[ -n "$chimeric" ]]; then
    echo "***! Chimeric stage requires already aligned files in ${splitdir}"  
    exit 1 
else
    mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 1; }
fi


## Create output directory, only if not in postproc, dedup or final stages
if [[ -d "$outputdir" && -z "$final" && -z "$afterdedup" && -z "$dedup" && -z "$postproc" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 1			
else
    if [[ -z "$final" && -z "$dedup" && -z "$postproc" && -z "$afterdedup" ]]; then
        mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 1; } 
    fi
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$dedup" ] && [ -z "$postproc" ] && [ -z "$afterdedup" ] ; then
    mkdir "$tmpdir"
    chmod 777 "$tmpdir"
fi

## Create output directory, used for reporting commands output
if [ ! -d "$debugdir" ]; then
    mkdir "$debugdir"
    chmod 777 "$debugdir"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline
# If chunk size sent in, split. Otherwise check size before splitting
if [ $isVoltron -ne 1 ]
then
    if [ -z $splitme ]
    then
	fastqsize=$(ls -lgGL  ${fastqdir} | awk '{sum+=$3}END{print sum}')
	if [ "$fastqsize" -gt "2592410750" ]
	then
	    splitme=1
	fi
    fi
fi

if [ -z "$user" ]
then
    userstring=""
else
    userstring="#SBATCH -A $user"
fi

# Add header containing command executed and timestamp
if [ "$methylation" = 1 ]
then
    queuestring="#SBATCH -p weka"
else
    queuestring="#SBATCH -p $queue"
fi
jid=`sbatch <<- HEADER | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l 
        $userstring
	$queuestring
	#SBATCH -t 2
	#SBATCH -c 1
	#SBATCH -o $debugdir/head-%j.out
	#SBATCH -e $debugdir/head-%j.err
	#SBATCH -J "${groupname}_cmd"
	date
	${load_bwa}
	${load_java}
	${load_awk}

	# Experiment description
	if [ -n "${about}" ]
	then
		echo -ne 'Experiment description: ${about}; '
	else
		echo -ne 'Experiment description: '
	fi

	# Get version numbers of all software
	echo -ne "Juicer version $juicer_version;" 
	bwa 2>&1 | awk '\\\$1=="Version:"{printf(" BWA %s; ", \\\$2)}'
	if [ "$methylation" = 1 ]
	then
		activate conda
		$call_bwameth --version 2>&1 | awk '{printf("%s; ", \\\$0)}'
	fi
	echo -ne "$threads threads; "
	if [ -n "$splitme" ]
	then
		echo -ne "splitsize $splitsize; "
	fi  
	java -version 2>&1 | awk 'NR==1{printf("%s; ", \\\$0);}'
	${juiceDir}/scripts/juicer_tools -V 2>&1 | awk '\\\$1=="Juicer" && \\\$2=="Tools"{printf("%s; ", \\\$0);}'
	
	echo "$0 $@"
HEADER`
headfile="${debugdir}/head-${jid}.out"

## Record if we failed while aligning, so we don't waste time on other jobs
## Remove file if we're relaunching Juicer 
errorfile=${debugdir}/${groupname}_alignfail
if [ -f $errorfile ]
then
    rm $errorfile
fi

# Not in merge, dedup,  or final stage, i.e. need to split and align files.
if [ -z $merge ] && [ -z $final ] && [ -z $dedup ] && [ -z $postproc ] && [ -z $afterdedup ]
then
    if [ "$nofrag" -eq 0 ]
    then
	echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $refSeq with site file $site_file"
    else
        echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $refSeq with no fragment delimited maps."
    fi
    
    ## Split fastq files into smaller portions for parallelizing alignment 
    ## Do this by creating a text script file for the job on STDIN and then 
    ## sending it to the cluster	
    dependsplit="afterok"
    if [ ! $splitdirexists ]
    then
	echo "(-: Created $splitdir and $outputdir."
	if [ -n "$splitme" ]
        then
            for i in ${fastqdir}
            do
		filename=$(basename $i)
		filename=${filename%.*}      
                if [ -z "$gzipped" ]
                then	
		    jid=`sbatch <<- SPLITEND | egrep -o -e "\b[0-9]+$"
			#!/bin/bash -l
                        #SBATCH -p $queue
			#SBATCH -t $queue_time
			#SBATCH -c 1
			#SBATCH --mem=5G
			#SBATCH -o $debugdir/split-%j.out
			#SBATCH -e $debugdir/split-%j.err
			#SBATCH -J "${groupname}_split_${i}"
                        $userstring			
			date
			echo "Split file: $filename"
			split -a 3 -l $splitsize -d --additional-suffix=.fastq $i $splitdir/$filename
			date
SPLITEND`
		else
		    jid=`sbatch <<- SPLITEND | egrep -o -e "\b[0-9]+$"
			#!/bin/bash -l
			#SBATCH -p $queue
			#SBATCH -t $queue_time
			#SBATCH -c 1
			#SBATCH --mem=5G
			#SBATCH -o $debugdir/split-%j.out
			#SBATCH -e $debugdir/split-%j.err
			#SBATCH -J "${groupname}_split_${i}"
                        $userstring			
			date
			echo "Split file: $filename"
			zcat $i | split -a 3 -l $splitsize -d --additional-suffix=.fastq - $splitdir/$filename
			date
SPLITEND`
		fi
		dependsplit="$dependsplit:$jid"
                # if we split files, the splits are named .fastq
                read1=${splitdir}"/*${read1str}*.fastq"
	    done
	    
	    srun -c 1 -p "$queue" -t 1 -o $debugdir/wait-%j.out -e $debugdir/wait-%j.err -d $dependsplit -J "${groupname}_wait" sleep 1
        else
            cp -rs ${fastqdir} ${splitdir}
            wait
        fi
    else
        ## No need to re-split fastqs if they already exist
        echo -e "---  Using already created files in $splitdir\n"
	# unzipped files will have .fastq extension, softlinked gz 
        testname=$(ls -lgG ${splitdir} | awk '$7~/fastq$/||$7~/gz$/{print $7; exit}')

	if [[ -z "$chimeric" ]]
	then
            if [[ ${testname: -3} == ".gz" ]]
            then
		read1=${splitdir}"/*${read1str}*.fastq.gz"
            else
		read1=${splitdir}"/*${read1str}*.fastq"
            fi
	fi
    fi
    
    ## Launch job. Once split/move is done, set the parameters for the launch. 
    echo "(-: Starting job to launch other jobs once splitting is complete"
    
    ## Loop over all read1/read2 fastq files and create jobs for aligning.
    ## Then call chimeric script on aligned, sort individual
    ## Wait for splits to be individually sorted, then do a big merge sort.
    ## ARRAY holds the names of the jobs as they are submitted
    countjobs=0
    declare -a ARRAY
    declare -a JIDS
    declare -a TOUCH

    dependmerge="afterok"

    for i in ${read1}
    do
	ext=${i#*$read1str}
	name=${i%$read1str*} 
	# these names have to be right or it'll break
	name1=${name}${read1str}
	name2=${name}${read2str}	
	jname=$(basename "$name")${ext}


	# RG group; ID derived from paired-end name, sample and library can be user set
	if [ $singleend -eq 1 ]
	then
	    rg="@RG\\tID:${jname%%.fastq*}\\tSM:${sampleName}\\tPL:LS454\\tLB:${libraryName}"
	else
	    rg="@RG\\tID:${jname%%.fastq*}\\tSM:${sampleName}\\tPL:ILM\\tLB:${libraryName}"
	fi
	touchfile=${tmpdir}/${jname}
	if [ -z "$chimeric" ]
	then
            usegzip=0
            if [ "${ext: -3}" == ".gz" ]
            then
		usegzip=1
	    fi

	    # count ligations
	    jid=`sbatch <<- CNTLIG |  egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -t $queue_time
		#SBATCH -c 1
		#SBATCH -o $debugdir/count_ligation-%j.out
		#SBATCH -e $debugdir/count_ligation-%j.err
		#SBATCH -J "${groupname}_${jname}_Count_Ligation"
		#SBATCH --mem=5G
                $userstring			

		date
		export usegzip=${usegzip}; export name=${name}; export name1=${name1}; export name2=${name2}; export ext=${ext}; export ligation=${ligation}; export singleend=${singleend}; ${juiceDir}/scripts/countligations.sh
		date
CNTLIG`
	    dependcount="$jid"

	    # align fastqs
	    jid=`sbatch <<- ALGNR1 | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		$queuestring
		#SBATCH -o $debugdir/align1-%j.out
		#SBATCH -e $debugdir/align1-%j.err
		#SBATCH -t $queue_time
		#SBATCH -n 1
		#SBATCH -c $threads
		#SBATCH --ntasks=1
		#SBATCH --mem=$alloc_mem
		#SBATCH -J "${groupname}_align1_${jname}"
		#SBATCH --threads-per-core=1		
                $userstring			

		${load_bwa}
		# Align reads
		date
		if [ "$methylation" = 1 ]
		then
		   conda activate
		fi

		if [ $singleend -eq 1 ]
		then
		   if [ "$methylation" = 1 ]
		   then
		      # The -M flag is already used
		      echo "Running bwameth.py $threadstring -5 --do-not-penalize-chimeras --reference ${refSeq} --read-group $rg $name1$ext > $name$ext.sam"
		      $call_bwameth $threadstring -5 --do-not-penalize-chimeras --reference ${refSeq} --read-group '$rg' $name1$ext > $name$ext.sam
		   else
			echo "Running command bwa mem -5M $threadstring -R $rg $refSeq $name1$ext > $name$ext.sam"
			bwa mem -5M $threadstring -R '$rg' $refSeq $name1$ext > $name$ext.sam 
                   fi
		else
		   if [ "$methylation" = 1 ]
		   then
		      # --read-group '$rg' seems to cause problems
		      # The -M flag is already used by bwameth
		      echo "Running bwameth.py $threadstring -5SP --do-not-penalize-chimeras --read-group '$rg'  --reference ${refSeq} $name1$ext $name2$ext > $name$ext.sam"
		      $call_bwameth $threadstring -5SP --do-not-penalize-chimeras --read-group '$rg' --reference ${refSeq} $name1$ext $name2$ext > $name$ext.sam      
		   else 
			echo "Running command bwa mem -SP5M $threadstring -R $rg $refSeq $name1$ext $name2$ext > $name$ext.sam" 
			bwa mem -SP5M $threadstring -R '$rg' $refSeq $name1$ext $name2$ext > $name$ext.sam
		   fi
		fi
		if [ \$? -ne 0 ]
		then  
		       touch $errorfile
		       exit 1
		else
		       echo "(-: Mem align of $name$ext.sam done successfully"
		fi
		date
ALGNR1`
	    dependalign="afterok:$jid:$dependcount"
	else
	    name=${i%.sam}
	    ext=""
	    jid=`sbatch <<- CNTLINE |  egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -t $queue_time
		#SBATCH -c 1
		#SBATCH -o $debugdir/count_line-%j.out
		#SBATCH -e $debugdir/count_line-%j.err
		#SBATCH -J "${groupname}_${jname}_Count_Line"
		#SBATCH --mem=5G
		${load_awk}
		${load_samtools}
                $userstring			

		date
		echo -ne "0 " > ${name}${ext}_norm.txt.res.txt
		samtools flagstat $i | awk '\\\$0~/paired in sequencing/{print \\\$1*2; exit}' > ${i}_linecount.txt
		date
CNTLINE`

	    dependalign="afterok:$jid"
	fi

	# wait for alignment, chimeric read handling
	if [ "$site" != "none" ] && [ -e "$site_file" ] 
	then		
	    jid=`sbatch <<- MRGALL | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $long_queue
		#SBATCH -o $debugdir/merge-%j.out
		#SBATCH -e $debugdir/merge-%j.err
		#SBATCH --mem=10G
		#SBATCH -t $long_queue_time
		#SBATCH -c 1
		#SBATCH --ntasks=1
		#SBATCH -d $dependalign
		#SBATCH -J "${groupname}_merge_${jname}"
                #SBATCH --threads-per-core=1
                $userstring
		${load_awk}

		date
		if [ $singleend -eq 1 ]
		then
		    time awk -v stem=${name}${ext}_norm -v site_file=$site_file -v singleend=$singleend -f $juiceDir/scripts/chimeric_sam.awk $name$ext.sam > $name$ext.sam3
		else
		    time awk -v stem=${name}${ext}_norm -v site_file=$site_file -f $juiceDir/scripts/chimeric_sam.awk $name$ext.sam > $name$ext.sam3
		fi
		date
MRGALL`
	    dependalign="afterok:$jid"
	else
	    if [ $singleend -eq 1 ]
	    then
		jid=`sbatch <<- MRGALL1 | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $long_queue
		#SBATCH -o $debugdir/merge1-%j.out
		#SBATCH -e $debugdir/merge1-%j.err
		#SBATCH --mem=10G
		#SBATCH -t $long_queue_time
		#SBATCH -c 1
		#SBATCH --ntasks=1
		#SBATCH -d $dependalign
		#SBATCH -J "${groupname}_merge_${jname}"
                #SBATCH --threads-per-core=1
                $userstring
		${load_awk}
		#time awk -v maxcount=1000000 -f $juiceDir/scripts/calculate_insert_size.awk $name$ext.sam > $name$ext.insert_size
		#will need to combine chimeric_sam and adjust_insert_size 
		time awk -v stem=${name}${ext}_norm -v singleend=$singleend -f $juiceDir/scripts/chimeric_sam.awk $name$ext.sam > $name$ext.sam3
MRGALL1`
		dependalign="afterok:$jid"
	    else
		jid=`sbatch <<- MRGALL1 | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $long_queue
		#SBATCH -o $debugdir/merge1-%j.out
		#SBATCH -e $debugdir/merge1-%j.err
		#SBATCH --mem=10G
		#SBATCH -t $long_queue_time
		#SBATCH -c 1
		#SBATCH --ntasks=1
		#SBATCH -d $dependalign
		#SBATCH -J "${groupname}_merge_${jname}"
                #SBATCH --threads-per-core=1
                $userstring
		${load_awk}
		#time awk -v maxcount=1000000 -f $juiceDir/scripts/calculate_insert_size.awk $name$ext.sam > $name$ext.insert_size
		#will need to combine chimeric_sam and adjust_insert_size 
		time awk -v stem=${name}${ext}_norm -f $juiceDir/scripts/chimeric_sam.awk $name$ext.sam > $name$ext.sam2
MRGALL1`
		dependalign="afterok:$jid"
		jid=`sbatch <<- MRGALL3 | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $long_queue
		#SBATCH -o $debugdir/merge2-%j.out
		#SBATCH -e $debugdir/merge2-%j.err
		#SBATCH --mem=10G
		#SBATCH -t $long_queue_time
		#SBATCH -c 1
		#SBATCH --ntasks=1
		#SBATCH -d $dependalign
		#SBATCH -J "${groupname}_merge_${jname}"
                #SBATCH --threads-per-core=1
                $userstring
		${load_awk}

		time awk -v avgInsertFile=${name}${ext}_norm.txt.res.txt -f $juiceDir/scripts/adjust_insert_size.awk $name$ext.sam2 > $name$ext.sam3
MRGALL3`
		dependalign="afterok:$jid"
	    fi
	fi
	
	jid2=`sbatch <<- MRGALL2 | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p $long_queue
#SBATCH -o $debugdir/mergesort-%j.out
#SBATCH -e $debugdir/mergesort-%j.err
#SBATCH --mem-per-cpu=2G
#SBATCH -t $long_queue_time
#SBATCH -c $sortthreads
#SBATCH --ntasks=1
#SBATCH -d $dependalign 
#SBATCH -J "${groupname}_mergesort_${jname}"
#SBATCH --threads-per-core=1 
$userstring
${load_samtools}
#we should probably set the -m based on memory / num of threads
if time samtools sort -t cb -n -O SAM -@ $sortthreads -l 0 -m 2G $name$ext.sam3 >  ${name}${ext}.sam
then
   rm -f $name$ext.sam2 $name$ext.sam3
   touch $touchfile
else
   echo "***! Failure during chimera handling of $name${ext}"
   touch $errorfile
   exit 1
fi
MRGALL2`
	dependmerge="${dependmerge}:${jid2}"
	ARRAY[countjobs]="${groupname}_mergesort_${jname}"
	JIDS[countjobs]="${jid}"
	TOUCH[countjobs]="$touchfile"
        countjobs=$(( $countjobs + 1 ))

    done # done looping over all fastq split files
    
    # list of all jobs. print errors if failed    
    for (( i=0; i < $countjobs; i++ ))
    do
	f=${TOUCH[$i]}
	msg="***! Error in job ${ARRAY[$i]}  Type squeue -j ${JIDS[$i]} to see what happened"
	
	# check that alignment finished successfully
	jid=`sbatch <<- EOF
		#!/bin/bash -l
		#SBATCH -o $debugdir/aligncheck-%j.out
		#SBATCH -e $debugdir/aligncheck-%j.err
		#SBATCH -t $queue_time
		#SBATCH -p $queue
		#SBATCH -J "${groupname}_check"
		#SBATCH -d $dependmerge
                $userstring			

		date
		echo "Checking $f"
		if [ ! -e $f ]
		then
			echo $msg
			touch $errorfile
		fi
		date
EOF`
	jid=$(echo $jid | egrep -o -e "\b[0-9]+$")
	dependmergecheck="${dependmerge}:${jid}"
    done
fi  # Not in merge, dedup,  or final stage, i.e. need to split and align files.

# Not in final, dedup, or postproc
if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ] && [ -z $afterdedup ]
then
    if [ -z $merge ]
    then
	sbatch_wait="#SBATCH -d $dependmergecheck"
    else
        sbatch_wait=""
    fi
    
    # merge the sorted files into one giant file that is also sorted. jid=`sbatch <<- MRGSRT | egrep -o -e "\b[0-9]+$"
    
    if [ $isVoltron -eq 1 ]
    then  
	sbatch_time="#SBATCH -t 10080"
    else
	sbatch_time="#SBATCH -t 1440"
    fi
    if [ $isBCM -eq 1 ]
    then
	sbatch_cpu_alloc="#SBATCH -c 1"
	sbatch_mem_alloc="#SBATCH --mem=80G"
    else
	sbatch_cpu_alloc="#SBATCH -c 8"
	sbatch_mem_alloc="#SBATCH --mem=64G"
    fi

    jid=`sbatch <<- EOF
		#!/bin/bash -l
		#SBATCH -o $debugdir/fragmerge-%j.out
		#SBATCH -e $debugdir/fragmerge-%j.err
		${sbatch_mem_alloc}
		${sbatch_time}
		#SBATCH -p $long_queue
		${sbatch_cpu_alloc}
		#SBATCH -J "${groupname}_fragmerge"
		${sbatch_wait}
                $userstring			

		date
		if [ -f "${errorfile}" ]
		then
			echo "***! Found errorfile. Exiting." 
			exit 1 
		fi
		export LC_COLLATE=C
		if [ -d $donesplitdir ]
		then
			mv $donesplitdir/* $splitdir/.
		fi

		if ! samtools merge -c -t cb -n $sthreadstring -O SAM $outputdir/merged_sort.sam  $splitdir/*.sam
		then
			echo "***! Some problems occurred somewhere in creating sorted align files."
			touch $errorfile
			exit 1
		else
			echo "(-: Finished sorting all sorted files into a single merge."
		fi
		date
EOF`

    jid=$(echo $jid | egrep -o -e "\b[0-9]+$")
    dependmrgsrt="afterok:$jid"
fi

# Remove the duplicates from the big sorted file
if [ -z $final ] && [ -z $postproc ] && [ -z $afterdedup ]
then
    if [ -z $dedup ]
    then
        sbatch_wait="#SBATCH -d $dependmrgsrt"
    else
        sbatch_wait=""
    fi
    # Guard job for dedup. this job is a placeholder to hold any job submitted after dedup.
    # We keep the ID of this guard, so we can later alter dependencies of inner dedupping phase.
    # After dedup is done, this job will be released. 
    guardjid=`sbatch <<- DEDUPGUARD | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/dedupguard-%j.out
	#SBATCH -e $debugdir/dedupguard-%j.err
	#SBATCH -t 10
	#SBATCH -c 1
	#SBATCH -H
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_dedup_guard"
	${sbatch_wait}
        $userstring			

	date
DEDUPGUARD`

    dependguard="afterok:$guardjid"

    # if jobs succeeded, kill the cleanup job, remove the duplicates from the big sorted file
    jid=`sbatch <<- DEDUP | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $debugdir/dedup-%j.out
	#SBATCH -e $debugdir/dedup-%j.err
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_dedup"
	${sbatch_wait}
        $userstring
	
	${load_awk}
	date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
	squeue -u $USER -o "%A %T %j %E %R" | column -t
	awk -v queue=$long_queue -v groupname=$groupname -v debugdir=$debugdir -v dir=$outputdir -v topDir=$topDir -v juicedir=$juiceDir -v site=$site -v genomeID=$genomeID -v genomePath=$genomePath -v user=$USER -v guardjid=$guardjid -v justexact=$justexact -v wobbleDist=$wobbleDist -f $juiceDir/scripts/split_rmdups_sam.awk $outputdir/merged_sort.sam

	##Schedule new job to run after last dedup part:
	##Push guard to run after last dedup is completed:
	##srun --ntasks=1 -c 1 -p "$queue" -t 1 -o ${debugdir}/dedup_requeue-%j.out -e ${debugdir}/dedup-requeue-%j.err -J "$groupname_msplit0" -d singleton echo ID: $ echo "\${!SLURM_JOB_ID}"; scontrol update JobID=$guardjid dependency=afterok:\$SLURM_JOB_ID
	squeue -u $USER -o "%A %T %j %E %R" | column -t
	date
	
	scontrol release $guardjid
DEDUP`

    dependosplit="afterok:$jid"

    #Push dedup guard to run only after dedup is complete:
    scontrol update JobID=$guardjid dependency=afterok:$jid

    #Wait for all parts of split_rmdups to complete:
    jid=`sbatch <<- MSPLITWAIT | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/post_dedup-%j.out
	#SBATCH -e $debugdir/post_dedup-%j.err
	#SBATCH -t 100
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_post_dedup"
	#SBATCH -d ${dependguard}
        $userstring			

	date
	rm -Rf $tmpdir;
	find $debugdir -type f -size 0 | xargs rm
	squeue -u $USER -o "%A %T %j %E %R" | column -t
	date
MSPLITWAIT`

    dependmsplit="afterok:$jid"
    sbatch_wait="#SBATCH -d $dependmsplit"
else
    sbatch_wait=""
fi

if [ -z "$genomePath" ]
then
    #If no path to genome is give, use genome ID as default.
    genomePath=$genomeID
fi

#Skip if only final or post-processing only is required
if [ -z $postproc ] && [ -z $final ]
    then
    # Check that dedupping worked properly
    awkscript='BEGIN{sscriptname = sprintf("%s/.%s_rmsplit.slurm", debugdir, groupname);}NR==1{if (NF == 2 && $1 == $2 ){print "Sorted and dups/no dups files add up"; printf("#!/bin/bash -l\n#SBATCH -o %s/dup-rm.out\n#SBATCH -e %s/dup-rm.err\n#SBATCH -p %s\n#SBATCH -J %s_msplit0\n#SBATCH -d singleton\n#SBATCH -t 1440\n#SBATCH -c 1\n#SBATCH --ntasks=1\ndate;\nrm %s/*_msplit*; rm %s/split*;\n rm %s/merged_sort.sam;\ndate\n", debugdir, debugdir, queue, groupname, dir, dir, dir) > sscriptname; sysstring = sprintf("sbatch %s", sscriptname); system(sysstring);close(sscriptname); }else{print "Problem"; print "***! Error! The sorted file and dups/no dups files do not add up, or were empty."; exit 1}}'
    jid=`sbatch <<- DUPCHECK | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/dupcheck-%j.out
	#SBATCH -e $debugdir/dupcheck-%j.err
	#SBATCH -t $queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=8G
	#SBATCH -J "${groupname}_dupcheck"
	${sbatch_wait}
        $userstring			
	${load_awk}
	date 
	wc -l ${outputdir}/merged_sort.sam |  awk '{printf("%s ", \\\$1)}' > $debugdir/dupcheck-${groupname}
	wc -l ${outputdir}/merged_dedup.sam | awk '{printf("%s ", \\\$1)}' >> $debugdir/dupcheck-${groupname}
	cat $debugdir/dupcheck-${groupname}
	awk -v debugdir=$debugdir -v queue=$queue -v groupname=$groupname -v dir=$outputdir '$awkscript' $debugdir/dupcheck-${groupname}
DUPCHECK`

    sbatch_wait="#SBATCH -d afterok:$jid"
    jid1=`sbatch <<- MERGED1 | egrep -o -e "\b[0-9]+$" 
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/merged1-%j.out
	#SBATCH -e $debugdir/merged1-%j.err
	#SBATCH -t $queue_time
	${sbatch_cpu_alloc}
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=10G
	#SBATCH -J "${groupname}_merged1"
	${sbatch_wait}
	$userstring
	${load_samtools}

	samtools view -F 1024 -O sam $sthreadstring ${outputdir}/merged_dedup.sam | awk -v mapq=1 -f ${juiceDir}/scripts/sam_to_pre.awk > ${outputdir}/merged1.txt
        date
MERGED1`

    sbatch_wait1="#SBATCH -d afterok:$jid1"
    jid2=`sbatch <<- MERGED30 | egrep -o -e "\b[0-9]+$" 
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/merged30-%j.out
	#SBATCH -e $debugdir/merged30-%j.err
	#SBATCH -t $queue_time
	${sbatch_cpu_alloc}
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=10G
	#SBATCH -J "${groupname}_merged30"
	${sbatch_wait}
	$userstring
	${load_samtools}

	samtools view -F 1024 -O sam $sthreadstring ${outputdir}/merged_dedup.sam | awk -v mapq=30 -f ${juiceDir}/scripts/sam_to_pre.awk > ${outputdir}/merged30.txt 
        date
MERGED30`
    sbatch_wait2="#SBATCH -d afterok:$jid2"

    sbatch_wait0="#SBATCH -d afterok:$jid1:$jid2"

    jid=`sbatch <<- PRESTATS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/prestats-%j.out
	#SBATCH -e $debugdir/prestats-%j.err
	#SBATCH -t $queue_time
	${sbatch_cpu_alloc} 
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=1G
	#SBATCH -J "${groupname}_prestats"
	${sbatch_wait}
	$userstring
	${load_awk}
        date
        ${load_java}
	${load_samtools}
        export IBM_JAVA_OPTIONS="-Xmx1024m -Xgcthreads1"
        export _JAVA_OPTIONS="-Xmx1024m -Xms1024m"
        tail -n1 $headfile | awk '{printf"%-1000s\n", \\\$0}' > $outputdir/inter.txt

	# count duplicates via samtools view -c
	# for paired end, count only first in pair to avoid double counting
	if [ $singleend -eq 1 ] 
	then 
		if [ -f $outputdir/merged_dedup.bam ]
		then
			samtools view $sthreadstring -f 1024 -F 256 $outputdir/merged_dedup.bam | awk '{if (\\\$0~/rt:A:7/){singdup++}else{dup++}}END{print dup,singdup}' > $outputdir/tmp
		else
			awk 'and(\\\$2,1024)>0 && and(\\\$2,256)==0{if (\\\$0~/rt:A:7/){singdup++}else{dup++}}END{print dup,singdup}' $outputdir/merged_dedup.sam > $outputdir/tmp 
		fi
		cat $splitdir/*.res.txt | awk -v fname=$outputdir/tmp -v ligation=$ligation -v singleend=1 -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt
	else 
		if [ -f $outputdir/merged_dedup.bam ] 
		then
			samtools view $sthreadstring -c -f 1089 -F 256 $outputdir/merged_dedup.bam > $outputdir/tmp
		else
		        awk 'and(\\\$2,1024)>0 && and(\\\$2,1)>0 && and(\\\$2,64)>0 && and(\\\$2,256)==0{dup++}END{print dup}' $outputdir/merged_dedup.sam > $outputdir/tmp
		fi
		cat $splitdir/*.res.txt | awk -v fname=$outputdir/tmp -v ligation=$ligation -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt
	fi
        cp $outputdir/inter.txt $outputdir/inter_30.txt

        date
PRESTATS`
    sbatch_wait000="${sbatch_wait1}:$jid"

    jid=`sbatch <<- BAMRM  | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $debugdir/bamrm-%j.out
	#SBATCH -e $debugdir/bamrm-%j.err
	#SBATCH -t $queue_time
	${sbatch_cpu_alloc}
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=10G
	#SBATCH -J "${groupname}_bamrm"
	${sbatch_wait0}
	$userstring
	${load_samtools}
	if samtools view -b $sthreadstring ${outputdir}/merged_dedup.sam > ${outputdir}/merged_dedup.bam
	then
		rm ${outputdir}/merged_dedup.sam
	fi
BAMRM`

    if [ "$methylation" = 1 ]
    then
	tmpj=`sbatch <<- METH | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p commons
		#SBATCH -o $debugdir/meth-%j.out
		#SBATCH -e $debugdir/meth-%j.err
		#SBATCH -t $queue_time
		${sbatch_cpu_alloc}
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=10G
		#SBATCH -J "${groupname}_meth"
		#SBATCH -d afterok:$jid
		$userstring
		${load_samtools}                            
		export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs0/home/neva/lib
		samtools sort $sthreadstring ${outputdir}/merged_dedup.bam > ${outputdir}/merged_dedup_sort.bam
		/gpfs0/home/neva/bin/MethylDackel extract -F 1024 --keepSingleton --keepDiscordant $refSeq ${outputdir}/merged_dedup_sort.bam 
		/gpfs0/home/neva/bin/MethylDackel extract -F 1024 --keepSingleton --keepDiscordant --cytosine_report --CHH --CHG $refSeq ${outputdir}/merged_dedup_sort.bam
		${juiceDir}/scripts/conversion.sh ${outputdir}/merged_dedup_sort.cytosine_report.txt > ${outputdir}/conversion_report.txt
		rm ${outputdir}/merged_dedup_sort.bam ${outputdir}/merged_dedup_sort.cytosine_report.txt*
METH`
    fi

    sbatch_wait00="${sbatch_wait2}:$jid"
    jid=`sbatch <<- STATS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/stats-%j.out
	#SBATCH -e $debugdir/stats-%j.err
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem=25G
	#SBATCH -J "${groupname}_stats"
	${sbatch_wait000}
        $userstring			

	date
	if [ -f "${errorfile}" ]
	then 
		echo "***! Found errorfile. Exiting." 
		exit 1 
	fi
	if [ $assembly -eq 1 ]
        then
		${juiceDir}/scripts/juicer_tools statistics $site_file $outputdir/inter.txt $outputdir/merged1.txt none
	else
		${juiceDir}/scripts/juicer_tools statistics $site_file $outputdir/inter.txt $outputdir/merged1.txt $genomePath
	fi
	date
STATS`
    sbatch_wait1="#SBATCH -d afterok:$jid"

    dependstats="afterok:$jid"
    jid=`sbatch <<- STATS30 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/stats30-%j.out
	#SBATCH -e $debugdir/stats30-%j.err
	#SBATCH -t $long_queue_time
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem=25G
	#SBATCH -J "${groupname}_stats30"
	${sbatch_wait00}
	$userstring

	date
	if [ $assembly -eq 1 ]
        then
		${juiceDir}/scripts/juicer_tools statistics $site_file $outputdir/inter_30.txt $outputdir/merged30.txt none
	else
	${juiceDir}/scripts/juicer_tools statistics $site_file $outputdir/inter_30.txt $outputdir/merged30.txt $genomePath
	fi
	date
STATS30`

    dependstats30="afterok:$jid"
    sbatch_wait1="${sbatch_wait1}:$jid"
    sbatch_waitstats="#SBATCH -d $dependstats"
    sbatch_waitstats30="#SBATCH -d $dependstats30"
else
    sbatch_wait1=""
    sbatch_waitstats=""
    sbatch_waitstats30=""
fi

if [ -z $postproc ] 
then
    # if early exit, we stop here, once the stats are calculated
    if [ ! -z "$earlyexit" ]
    then
	if [ $assembly -eq 1 ]
	then
	    jid=`sbatch <<- MND | egrep -o -e "\b[0-9]+$" 
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem=2G
	#SBATCH -o $debugdir/mnd-%j.out
	#SBATCH -e $debugdir/mnd-%j.err
	#SBATCH -t 1200
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_mnd"     
	${sbatch_wait1}
        $userstring	   
	${load_samtools}
	date

	samtools view $sthreadstring -O SAM -F 1024 $outputdir/merged_dedup.*am | awk -v mnd=1 -f ${juiceDir}/scripts/sam_to_pre.awk > ${outputdir}/merged_nodups.txt 
	date
MND`
	    sbatch_wait1="#SBATCH -d afterok:$jid"
	fi

	jid=`sbatch <<- FINCLN1 | egrep -o -e "\b[0-9]+$" 
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem=2G
	#SBATCH -o $debugdir/fincln1-%j.out
	#SBATCH -e $debugdir/fincln1-%j.err
	#SBATCH -t 1200
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_prep_done"     
        #SBATCH --mail-type=END,FAIL
	${sbatch_wait1}
        $userstring	   
	date
	export splitdir=${splitdir}; export outputdir=${outputdir}; export early=1; 
	if ${juiceDir}/scripts/check.sh 
	then 
           if [ "$cleanup" = 1 ] 
           then
	      ${juiceDir}/scripts/cleanup.sh
           fi
	fi
	date
FINCLN1`
	echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id $jid"
	exit 0
    fi
    
    if [ "$qc" = 1 ] || [ "$insitu" = 1 ]
    then
        resstr="-r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000"
    else
        resstr="-r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100"
    fi
    if [ "$nofrag" -eq 1 ]
    then
	fragstr=""
    else
        fragstr="-f $site_file"
    fi

    jid=`sbatch <<- HIC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/hic-%j.out
	#SBATCH -e $debugdir/hic-%j.err	
	#SBATCH -t $long_queue_time
	#SBATCH -c $threadsHic
	#SBATCH --ntasks=1
	#SBATCH --mem=150G
	#SBATCH -J "${groupname}_hic"
	${sbatch_waitstats}
        $userstring			

	${load_java}
	export IBM_JAVA_OPTIONS="-Xmx150000m -Xgcthreads1"
	export _JAVA_OPTIONS="-Xmx150000m -Xms150000m"
	date
	if [ -f "${errorfile}" ]
	then 
		echo "***! Found errorfile. Exiting." 
		exit 1 
	fi 
	mkdir ${outputdir}"/HIC_tmp"

	# multithreaded and index doesn't exist yet
	if [[ $threadsHic -gt 1 ]] && [[ ! -s ${outputdir}/merged1_index.txt ]] 
	then
	  time ${juiceDir}/scripts/index_by_chr.awk ${outputdir}/merged1.txt 500000 > ${outputdir}/merged1_index.txt
	fi

	time ${juiceDir}/scripts/juicer_tools pre -n -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $resstr $fragstr $threadHicString $outputdir/merged1.txt $outputdir/inter.hic $genomePath
	time ${juiceDir}/scripts/juicer_tools addNorm $threadNormString ${outputdir}/inter.hic 
	rm -Rf ${outputdir}"/HIC_tmp"
	date
HIC`

    dependhic="afterok:$jid"

    jid=`sbatch <<- HIC30 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $debugdir/hic30-%j.out
	#SBATCH -e $debugdir/hic30-%j.err
	#SBATCH -t $long_queue_time
	#SBATCH -c $threadsHic
	#SBATCH --ntasks=1
	#SBATCH --mem=150G
	#SBATCH -J "${groupname}_hic30"
	${sbatch_waitstats30}
        $userstring	

	${load_java}
	export IBM_JAVA_OPTIONS="-Xmx150000m -Xgcthreads1"
	export _JAVA_OPTIONS="-Xmx150000m -Xms150000m"
	date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
	mkdir ${outputdir}"/HIC30_tmp"
	# multithreaded and index doesn't exist yet
	if [[ $threadsHic -gt 1 ]] && [[ ! -s ${outputdir}/merged30_index.txt ]] 
	then
	  time ${juiceDir}/scripts/index_by_chr.awk ${outputdir}/merged30.txt 500000 > ${outputdir}/merged30_index.txt
	fi

	time ${juiceDir}/scripts/juicer_tools pre -n -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $resstr $fragstr $threadHic30String $outputdir/merged30.txt $outputdir/inter_30.hic $genomePath
	time ${juiceDir}/scripts/juicer_tools addNorm $threadNormString ${outputdir}/inter_30.hic
	rm -Rf ${outputdir}"/HIC30_tmp"
	date
HIC30`

    dependhic30="${dependhic}:$jid"
    sbatch_wait="#SBATCH -d $dependhic30"
else
    sbatch_wait=""
fi

if [[ "$isNots" -eq 1 ]] || [[ "$isVoltron" -eq 1 ]]
then
    if [[  "$isNots" -eq 1 ]] 
    then
	sbatch_req="#SBATCH --gres=gpu:kepler:1"
    fi
    if [ "$qc" != 1 ]
    then
	jid=`sbatch <<- HICCUPS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=4G
	${sbatch_req}
	#SBATCH -o $debugdir/hiccups_wrap-%j.out
	#SBATCH -e $debugdir/hiccups_wrap-%j.err
	#SBATCH -t $queue_time
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_hiccups_wrap"
	${sbatch_wait}
        $userstring

	${load_gpu}
	echo "load: $load_gpu"
	${load_java}
	date
	nvcc -V
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
	${juiceDir}/scripts/juicer_hiccups.sh -j ${juiceDir}/scripts/juicer_tools -i $outputdir/inter_30.hic -m ${juiceDir}/references/motif -g $genomeID
	date
HICCUPS`
	dependhiccups="afterok:$jid"
    fi
else
    dependhiccups="afterok"
fi

if [ "$qc" != 1 ]
then
    jid=`sbatch <<- ARROWS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=8G
	#SBATCH -o $debugdir/arrowhead_wrap-%j.out
	#SBATCH -e $debugdir/arrowhead_wrap-%j.err
	#SBATCH -t $queue_time
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_arrowhead_wrap"
	${sbatch_wait}
        $userstring			

	${load_java}
	date
        if [ -f "${errorfile}" ]
        then 
            echo "***! Found errorfile. Exiting." 
            exit 1 
        fi 
	${juiceDir}/scripts/juicer_arrowhead.sh -j ${juiceDir}/scripts/juicer_tools -i $outputdir/inter_30.hic
	date;
ARROWS`
    dependarrows="#SBATCH -d ${dependhiccups}:$jid"
else
    dependarrows=${sbatch_wait}
fi

if [ "$qc_apa" = 1 ]
then
    jid=`sbatch <<- QC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=4G
	#SBATCH -o $debugdir/qc-%j.out
	#SBATCH -e $debugdir/qc-%j.err
	#SBATCH -t $queue_time
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_qc_apa"
	${sbatch_wait}
	$userstring
	${load_java}
	date
	export IBM_JAVA_OPTIONS="-Xmx4000m -Xgcthreads1"                                                                    
        export _JAVA_OPTIONS="-Xmx4000m -Xms4000m" 
	java -jar ${juiceDir}/scripts/juicer_3.25.21_aggNormAPA.jar apa --ag-norm -k NONE -n 300 -w 100 -r 1000 -q 60 --threads 1 $outputdir/inter_30.hic ${juiceDir}/scripts/GSE63525_GM12878_primary_replicate_HiCCUPS_looplist_with_motifs_unique_localized.txt $outputdir/qc_apa
	date
QC`
fi

jid=`sbatch <<- FINCLN1 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $debugdir/fincln-%j.out
	#SBATCH -e $debugdir/fincln-%j.err
	#SBATCH -t 1200
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_prep_done"
	$dependarrows
        $userstring			

	date
	export splitdir=${splitdir}
	export outputdir=${outputdir}


	if ${juiceDir}/scripts/check.sh 
        then
           if [ "$cleanup" = 1 ] 
           then
	      ${juiceDir}/scripts/cleanup.sh
           fi
	fi
	date
FINCLN1`

echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id $jid"
