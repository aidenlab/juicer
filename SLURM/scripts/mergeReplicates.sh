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
# arguments (default human, DpnII). Optional arguments are the queue for the 
# alignment (default short), key for menu entry, description for stats file, 
# using the short read aligner, read end (to align one read end using short 
# read aligner), early exit (to exit before the final creation of .hic files),
# merge (to start jobs at merge step), and the top-level directory (default 
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
#             are so many jobs that the cluster won't run them
# read1str  - portion of fastq filename that indicates this is the "read 1"
#             file; used to loop over only the read 1 and within that loop,
#             also align read 2 and merge.  If this is not set correctly,
#             script will not work. The error will often manifest itself
#             through a "*" in the name because the wildcard was not able to
#             match any files with the read1str.   
shopt -s extglob

## Set the following variables to work with your system
# path additionals, make sure paths are correct for your system
#myPath=/bin:$PATH
# set global tmpdir so no problems with /var/tmp
## use cluster load commands:
#usePath=""
load_bwa="module load biobuilds"
#load_java="module load Java/7.1.2.0"
#load_cluster=""
#load_coreutils=""
load_gpu="module load cuda/7.5"
# Juicer directory, contains scripts/. references/, and restriction_sites/ moved to referencegenomes
juiceDir="/home/keagen/programs/juicer/SLURM"
# Reference genome direction, contains fasta files and restriction site files somewhere in here
refDir="/scratch/PI/kornberg/keagen/referencegenomes/fly"
# default queue, can also be set in options
queue="kornberg,owners,bigmem,normal"
# default long queue, can also be set in options
long_queue="kornberg,owners,normal"

# top level directory, can also be set in options
topDir=$(pwd)
# unique name for jobs in this run
groupname=$(basename $topDir)
# parent directory
parentdir=$(dirname $topDir)
# experiment name, needed to identify replicates, should be related to groupname.  replicates should be lettered (e.g. a, b, c, etc)
experiment="${groupname//merge}"
# Find the replicates (i.e. any directories with same experiment name but without "merge" in it:
replicates=$(ls $parentdir | grep $experiment | grep -v "merge")

## Default options, overridden by command line arguments

#output messages
outDir="$topDir/debug"
# restriction enzyme, can also be set in options
site="DpnII"
# genome ID, default to human, can also be set in options
genomeID="dm3"
# normally both read ends are aligned with long read aligner; 
# if one end is short, this is set                 
shortreadend=0
# description, default empty
about="Drosophila melanogaster Kc167 cultured cells; merged 2 biological replicates"

## Read arguments
usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-q queue] [-l long queue] [-s site] [-a about] [-R end] [-S stage] [-p chrom.sizes path]\n\t\t[-y restriction site file path] [-z reference sequence path] [-r] [-h]"
genomeHelp="   genomeID must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \"$genomeID\")\n   alternatively, it can be defined using the -z command"
dirHelp="   [topDir] is the top level directory (default \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
queueHelp="   [queue] is the queue for running alignments (default \"$queue\")"
longQueueHelp="   [long queue] is the queue for running longer jobs such as the hic file creation (default \"$long_queue\")"
siteHelp="   [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" (default \"$site\")"
shortHelp2="   [end]: use the short read aligner on read end, must be one of 1 or 2 "
stageHelp="   [stage]: must be one of \"merge\", \"dedup\", \"final\", or \"early\".\n\t\tUse \"merge\" when alignment has finished but the merged_sort file has not yet been created.\n\t\tUse \"dedup\" when the files have been merged into merged_sort but merged_nodups has not yet been created.\n\t\tUse \"final\" when the reads have been deduped into merged_nodups and the final stats have been calculated but the hic files have not yet been created.\n\t\tUse \"early\" for an early exit, before the final creation of the stats and hic files"
shortHelp="   -r: use the short read version of the aligner, bwa aln (default: long read, bwa mem)"
aboutHelp="   -a: enter description of experiment, enclosed in single quotes"
pathHelp="   -p: enter path for chrom.sizes file"
refSeqHelp="   -z: enter path for reference sequence file, BWA index file must be in same directory"
siteFileHelp="   -y: enter path for restriction site file (locations of restriction sites in genome)"
helpHelp="   -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$genomeHelp"
    echo -e "$dirHelp"
    echo "$queueHelp"
    echo "$longQueueHelp"
    echo "$siteHelp"
    echo "$shortHelp2"
    echo "$shortHelp"
    echo -e "$stageHelp"
    echo "$aboutHelp"
    echo "$pathHelp"
    echo "$siteFileHelp"
    echo "$refSeqHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:R:k:a:hrmefq:s:p:l:y:z:S:" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	l) long_queue=$OPTARG ;;
	q) queue=$OPTARG ;;
	s) site=$OPTARG ;;
	R) shortreadend=$OPTARG ;;
	r) shortread=1 ;;  #use short read aligner
	a) about=$OPTARG ;;
	p) genomePath=$OPTARG ;;  
	y) site_file=$OPTARG ;;
	z) refSeq=$OPTARG ;;
	S) stage=$OPTARG ;;
	[?]) printHelpAndExit 1;;
    esac
done

if [ ! -z "$stage" ]
then
    case $stage in
        merge) merge=1 ;;
        dedup) dedup=1 ;;
        early) earlyexit=1 ;;
        final) final=1 ;; 
        *)  echo "$usageHelp"
			echo "$stageHelp"
			exit 1
	esac
fi

## Set reference sequence based on genome ID
if [ -z "$refSeq" ]
then 
    case $genomeID in
		mm9) refSeq="${juiceDir}/references/Mus_musculus_assembly9_norandom.fasta";;
		mm10) refSeq="${juiceDir}/references/Mus_musculus_assembly10.fasta";;
		hg38) refSeq="${juiceDir}/references/hg38.fa";;
		hg19) refSeq="${juiceDir}/references/Homo_sapiens_assembly19.fasta";;
		dm3) refSeq="${refDir}/${genomeID}/BWA/dm3.fa";;
		dm6) refSeq="${refDir}/${genomeID}/BWA/dm6.fa";;


    *)  echo "$usageHelp"
        echo "$genomeHelp"
        exit 1
    esac
    
   	## Also set genomePath since genomeID is set (added by K. Eagen).
   	genomePath="/scratch/PI/kornberg/keagen/referencegenomes/fly/${genomeID}/${genomeID}.chrom.sizes"

else
    ## Reference sequence passed in, so genomePath must be set for the .hic file
    ## to be properly created
    if [ -z "$genomePath" ]
        then
        echo "***! You must define a chrom.sizes file via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq";
        exit 100;
    fi
fi

## Check that refSeq exists 
if [ ! -e "$refSeq" ]; then
    echo "***! Reference sequence $refSeq does not exist";
    exit 100;
fi

## Check that index for refSeq exists
if [ ! -e "${refSeq}.bwt" ]; then
    echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run bwa index on this file before running juicer.";
    exit 100;
fi

## Set ligation junction based on restriction enzyme
case $site in
	HindIII) ligation="AAGCTAGCTT";;
	DpnII) ligation="GATCGATC";;
	MboI) ligation="GATCGATC";;
	NcoI) ligation="CCATGCATGG";;
	*)  ligation="XXXX"
		echo "$site not listed as recognized enzyme. Using $site_file as site file"
		echo "Ligation junction is undefined"
		exit 100
esac

## If short read end is set, make sure it is 1 or 2
case $shortreadend in
	0) ;;
	1) ;;
	2) ;;
	*)	echo "$usageHelp"
		echo "$shortHelp2"
		exit 100
esac

if [ -z "$site_file" ]
then
    site_file="${refDir}/${genomeID}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [ ! -e "$site_file" ]; then
    echo "***! $site_file does not exist. It must be created before running this script."
    exit 100
fi

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/juicer_splits"
outputdir=${topDir}"/juicer_output"
tmpdir=${topDir}"/juicer_HIC_tmp"
dumpdir=${outputdir}"/dumped_data"


## Create output directory, only if not in dedup or final stages
if [[ -d "$outputdir" && -z "$final" && -z "$dedup" ]] 
then
	echo "***! Move or remove directory \"$outputdir\" before proceeding."
	echo "***! Type \"juicer.sh -h \" for help"
	exit 100			
else
    if [[ -z "$final" && -z "$dedup" ]]; then
        mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 100; } 
    fi
fi

## Create split directory
if [ -d "$splitdir" ]; then
	splitdirexists=1
else
	mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 100; }
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$dedup" ]; then
	mkdir "$tmpdir"
	chmod 777 "$tmpdir"
fi

## Create output directory, used for reporting commands output
if [ ! -d "$outDir" ]; then
        mkdir "$outDir"
        chmod 777 "$outDir"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

srun --ntasks=1 -c 1 -p "$queue" -t 1 -o ${outDir}/head-%j.out -e ${outDir}/head-%j.err -J "${groupname}_cmd" echo "$0 $@"

if [ -z $final ] && [ -z $postproc ]
then

jid=`sbatch <<- MERGE_REPLICATES | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $outDir/merge_replicates-%j.out
	#SBATCH -e $outDir/merge_replicates-%j.err
	#SBATCH -t 1440
	#SBATCH -c 8
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=15500M
	#SBATCH -J "${groupname}_merge_replicates"
	
	module load coreutils

	if ! sort --parallel=8 -S 120G -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $parentdir/$experiment*/juicer_output/merged_nodups.txt > $outputdir/merged_nodups.txt
	#if ! sort -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $parentdir/$experiment*/juicer_output/merged_nodups.txt > $outputdir/merged_nodups.txt
	then 
		echo "***! Some problems occurred somewhere in creating sorted merged_nodups.txt file."
		exit 22
	else
		echo "(-: Finished sorting all merged_nodups.txt files into a single, sorted merged_nodups.txt file."	
	fi

	if ! sort --parallel=8 -S 120G -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $parentdir/$experiment*/juicer_output/dups.txt > $outputdir/dups.txt
	#if ! sort -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $parentdir/$experiment*/juicer_output/dups.txt > $outputdir/dups.txt
	then 
		echo "***! Some problems occurred somewhere in creating sorted dups.txt file."
		exit 22
	else
		echo "(-: Finished sorting all dups.txt files into a single, sorted dups.txt file."	
	fi

	if ! sort --parallel=8 -S 120G -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $parentdir/$experiment*/juicer_output/opt_dups.txt > $outputdir/opt_dups.txt
	#if ! sort -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $parentdir/$experiment*/juicer_output/opt_dups.txt > $outputdir/opt_dups.txt
	then 
		echo "***! Some problems occurred somewhere in creating sorted opt_dups.txt file."
		exit 22
	else
		echo "(-: Finished sorting all opt_dups.txt files into a single, sorted opts_dups.txt file."	
	fi
	
	cat $parentdir/$experiment*/juicer_splits/*.res.txt > $outputdir/$experiment.res.txt
	
MERGE_REPLICATES`

dependmerge_replicates="afterok:$jid"

# Set the number of chunks.  Can change this, but may also need to change memory requirements for randomization jobs below.
chunks=256

jid=`sbatch <<- SPLIT_MERGED_NODUPS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $outDir/split_merged_nodups-%j.out
	#SBATCH -e $outDir/split_merged_nodups-%j.err
	#SBATCH -t 16:00:00
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=16G
	#SBATCH -J "${groupname}_split_merged_nodups"
	#SBATCH -d ${dependmerge_replicates}

	module load coreutils

	# Split into chunks (number of chunks set above), and need to make sure lines don't get split apart.
	srun -N 1 -n 1 split -n l/$chunks -a 3 -d --additional-suffix=.txt $outputdir/merged_nodups.txt $splitdir/merged_nodups_

SPLIT_MERGED_NODUPS`

dependsplitmerged_nodups="afterok:$jid"

dependrandomizepos="afterok"

for i in $(seq -w 0 $(expr $chunks - 1)); do
	jid=`sbatch <<- RANDOMIZE_POSITIONS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $outDir/randomize_positions-$i-%j.out
	#SBATCH -e $outDir/randomize_positions-$i-%j.err
	#SBATCH -t 16:00:00
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=15G
	#SBATCH -J "${groupname}_randomize_positions-$i"
	#SBATCH -d ${dependsplitmerged_nodups}

	srun python ${juiceDir}/scripts/randomizePositions.py $i $site_file $splitdir
	
	RANDOMIZE_POSITIONS`
	dependrandomizepos="${dependrandomizepos}:$jid"
done

jid=`sbatch <<- MERGE_RANDS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $outDir/merge_randoms-%j.out
	#SBATCH -e $outDir/merge_randoms-%j.err
	#SBATCH -t 1440
	#SBATCH -c 8
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=15500M
	#SBATCH -J "${groupname}_merge_randoms"
	#SBATCH -d ${dependrandomizepos}

	module load coreutils
	
	if ! srun -N 1 -n 1 -c 8 sort --parallel=8 -S 120G -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/merged_nodups_randomized_*.txt  > $outputdir/merged_nodups_randomized.txt
	#if ! sort -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/merged_nodups_randomized_*.txt  > $outputdir/merged_nodups_randomized.txt

	then 
		echo "***! Some problems occurred somewhere in creating merged, sorted merged_nodups_randomized.txt file."
		exit 22
	else
		echo "(-: Finished sorting all merged_nodups_randomized.txt files into a single merge."
			rm -Rf $tmpdir;
	fi
	
	# Also finish intrafragment statistics	
	cat $splitdir/intrafragment_statistics_*.txt > $outputdir/intrafragment_statistics.tmp
	srun python ${juiceDir}/scripts/intrafragStats.py $outputdir
	rm $outputdir/intrafragment_statistics.tmp
	
MERGE_RANDS`

dependmergerands="afterok:$jid"

# Run preseq for complexity analysis/plots.  Must have preseq installed for this.

if [ ! -d "$topDir/preseq_output" ]; then
	mkdir $topDir/preseq_output
fi

jid=`sbatch <<- PRESEQ | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $outDir/preseq-%j.out
	#SBATCH -e $outDir/preseq-%j.err
	#SBATCH -t 16:00:00
	#SBATCH -n 1
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=125G
	#SBATCH -J "${groupname}_preseq"
	#SBATCH -d $dependmerge_replicates
	
	srun python ${juiceDir}/scripts/preseqReadCounts.py
		
	module load preseq
	preseq c_curve -v -V -s 100000 -o $topDir/preseq_output/${groupname}_c_curve.txt $topDir/preseq_output/obsReadCounts.txt

	preseq lc_extrap -v -e 1000000000000 -V -o $topDir/preseq_output/${groupname}_lc_extrap.txt $topDir/preseq_output/obsReadCounts.txt

	srun xvfb-run python ${juiceDir}/scripts/preseqPlots.py
	
PRESEQ`

jid=`sbatch <<- STATS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $outDir/stats-%j.out
	#SBATCH -e $outDir/stats-%j.err
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=125G
	#SBATCH -d $dependmerge_replicates
	#SBATCH -J "${groupname}_stats"
	$load_java
	export JAVA_OPTIONS="-Xmx16384m -XX:ParallelGCThreads=1"

	echo 'Experiment description: $about' > $outputdir/inter.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/stats_dups.txt $outputdir/dups.txt
	cat $outputdir/$experiment.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt
	java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt
	
	echo 'Experiment description: $about' > $outputdir/inter_30.txt
	cat $outputdir/$experiment.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter_30.txt
	java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter_30.txt >> $outputdir/inter_30.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt
	
	rm $outputdir/$experiment.res.txt
	
STATS`

dependstats="afterok:$jid"
dependmergerands="${dependmergerands}:$jid"

fi

if [ $final ] && [ $postproc ]
then
	sbatch_wait=""
	sbatch_wait_random=""

else
	sbatch_wait="#SBATCH -d $dependstats"
	sbatch_wait_random="#SBATCH -d $dependmergerands"

fi

if [ -z "$genomePath" ]
then
        #If no path to genome is give, use genome ID as default.
        genomePath=$genomeID
fi

jid=`sbatch <<- HIC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $outDir/hic-%j.out
    #SBATCH -e $outDir/hic-%j.err	
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=125G
	#SBATCH -J "${groupname}_hic"
	${sbatch_wait}
	${load_java}
	export JAVA_OPTIONS="-Xmx48192m -XX:ParallelGCThreads=1"

	${juiceDir}/scripts/juicebox48g pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/${groupname}_inter.hic $genomePath
HIC`

	dependhic="afterok:$jid"

jid=`sbatch <<- HIC30 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $outDir/hic30-%j.out
	#SBATCH -e $outDir/hic30-%j.err
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=125G
	#SBATCH -J "${groupname}_hic30"
	${sbatch_wait}
	${load_java}
	export JAVA_OPTIONS="-Xmx48192m -XX:ParallelGCThreads=1"

	${juiceDir}/scripts/juicebox48g pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/${groupname}_inter_30.hic $genomePath
HIC30`

dependhic30="${dependhic}:$jid"

jid=`sbatch <<- HIC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $outDir/hic_randomized-%j.out
    #SBATCH -e $outDir/hic_randomized-%j.err	
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=125G
	#SBATCH -J "${groupname}_hic_randomized"
	${sbatch_wait_random}
	${load_java}
	export JAVA_OPTIONS="-Xmx48192m -XX:ParallelGCThreads=1"

	${juiceDir}/scripts/juicebox48g pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups_randomized.txt $outputdir/${groupname}_randomizedpositions_inter.hic $genomePath
HIC`

dependhicrandomized="afterok:$jid"

jid=`sbatch <<- HIC30 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $outDir/hic30_randomized-%j.out
	#SBATCH -e $outDir/hic30_randomized-%j.err
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=125G
	#SBATCH -J "${groupname}_hic30_randomized"
	${sbatch_wait_random}
	${load_java}
	export JAVA_OPTIONS="-Xmx48192m -XX:ParallelGCThreads=1"

	${juiceDir}/scripts/juicebox48g pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups_randomized.txt $outputdir/${groupname}_randomizedpositions_inter_30.hic $genomePath
HIC30`

dependhic30randomized="${dependhicrandomized}:$jid"

if [ ! -d $dumpdir ]; then
	mkdir $dumpdir
fi

dependdump="afterok"

for res in 5000 2000 1000 500; do
	for chr in chr2L chr2R chr3L chr3R chrX chr4; do
		jid=`sbatch <<- DUMPDATA | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -o $outDir/dump_data_res${res}_${chr}-%j.out
		#SBATCH -e $outDir/dump_data_res${res}_${chr}-%j.err
		#SBATCH -t 60
		#SBATCH -c 1
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=125G
		#SBATCH -J "${groupname}_dump_data_res${res}_${chr}"
		#SBATCH -d $dependhic30randomized

		${juiceDir}/scripts/juicebox48g dump observed KR $outputdir/${groupname}_randomizedpositions_inter_30.hic $chr $chr BP $res $dumpdir/${chr}_observed_randomizedpositions_${res}.txt
		${juiceDir}/scripts/juicebox48g dump oe KR $outputdir/${groupname}_randomizedpositions_inter_30.hic $chr $chr BP $res | tail -n +2 | head -n -1 > $dumpdir/${chr}_oe_randomizedpositions_${res}.txt

		DUMPDATA`
		dependdump="${dependdump}:$jid"
	done
done

jid=`sbatch <<- COVERAGE | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH -o $outDir/coverage-%j.out
	#SBATCH -e $outDir/coverage-%j.err
	#SBATCH -t 60
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=16G
	#SBATCH -d $dependdump
	#SBATCH -J "${groupname}_coverage"

	python ${juiceDir}/scripts/coverageCalc.py $dumpdir $genomePath $outputdir
	
	COVERAGE`

dependarrowhead="afterok"
for res in 500 1000 2000 5000; do
	jid=`sbatch <<- ARROWHEAD | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=16G
	#SBATCH -o $outDir/arrowhead_${res}-%j.out
	#SBATCH -e $outDir/arrowhead_${res}-%j.err
	#SBATCH -t 60
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_arrowhead_${res}"
	#SBATCH -d $dependhic30randomized

	echo "Running Arrowhead for resolution "$res
	srun ${juiceDir}/scripts/juicebox arrowhead -m 2000 -k KR -r $res $outputdir/${groupname}_randomizedpositions_inter_30.hic $outputdir/${groupname}_randomizedpositions_inter_30_TADs_res${res}.txt
	
	ARROWHEAD`
	dependarrowhead="${dependarrowhead}:$jid"
done

jid=`sbatch <<- ARROWHEAD_POST | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=4G
	#SBATCH -o $outDir/arrowhead_post-%j.out
	#SBATCH -e $outDir/arrowhead_post-%j.err
	#SBATCH -t 60
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_arrowhead_post"
	#SBATCH -d $dependarrowhead

	srun python ${juiceDir}/scripts/arrowheadPost.py

	ARROWHEAD_POST`

jid=`sbatch <<- HICCUPS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p owners,gpu
	#SBATCH --mem-per-cpu=16G
	#SBATCH --gres=gpu:1
	#SBATCH -o $outDir/hiccups-%j.out
	#SBATCH -e $outDir/hiccups-%j.err
	#SBATCH -t 60
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_hiccups"
	#SBATCH -d $dependhic30randomized

	${load_java}
	${load_gpu}
	mkdir /local-scratch/$USER/juicer
	export TMPDIR=/local-scratch/keagen/juicer

	srun ${juiceDir}/scripts/juicebox hiccups $outputdir/${groupname}_randomizedpositions_inter_30.hic -m 512 -k KR \
	-r 2000 \
	-f 0.1 \
	-p 8 \
	-i 14 \
	-t 0.02,1.5,1.75,2 \
	-d 20000 \
	$outputdir/${groupname}_randomizedpositions_inter_30_loops.txt

HICCUPS`

echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee.."
