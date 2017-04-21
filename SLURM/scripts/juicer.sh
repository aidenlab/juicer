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
# arguments (default human, MboI). Optional arguments are the queue for the 
# alignment (default short), description for stats file, 
# using the short read aligner, read end (to align one read end using short 
# read aligner), stage to relaunch at, paths to various files if needed,
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
shopt -s extglob
juicer_version="1.5" 
## Set the following variables to work with your system
# path additionals, make sure paths are correct for your system
#myPath=/bin:$PATH
# set global tmpdir so no problems with /var/tmp
## use cluster load commands:
#usePath=""
load_bwa="module load BioBuilds/2015.04;"
load_java="module load Java/7.1.2.0;"
#load_cluster=""
#load_coreutils=""
load_gpu="module load gcccuda/2.6.10;export PATH=$PATH:/usr/local/cuda/bin/;"

# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir="/poscratch/aidenlab/juicer"
# default queue, can also be set in options via -q
queue="commons"
# default queue time, can also be set in options via -Q
queue_time="1200"
# default long queue, can also be set in options via -l
long_queue="commons"
# default long queue time, can also be set in options via -L
long_queue_time="3600"
# size to split fastqs. adjust to match your needs. 4000000=1M reads per split
# can also be changed via the -C flag
splitsize=90000000
# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1" 
read2str="_R2" 

# unique name for jobs in this run
groupname="a$(date +%s)"

## Default options, overridden by command line arguments

# top level directory, can also be set in options
topDir=$(pwd)
#output messages
outDir="$topDir/debug"
# restriction enzyme, can also be set in options
site="MboI"
# genome ID, default to human, can also be set in options
genomeID="hg19"
# normally both read ends are aligned with long read aligner; 
# if one end is short, this is set                 
shortreadend=0
# description, default empty
about=""
nofrag=0

## Read arguments                                                     
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]\n                 [-a about] [-R end] [-S stage] [-p chrom.sizes path]\n                 [-y restriction site file] [-z reference genome file]\n                 [-C chunk size] [-D Juicer scripts directory]\n                 [-Q queue time limit] [-L long queue time limit] [-r] [-h] [-x]"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \n  \"$genomeID\"); alternatively, it can be defined using the -z command"
dirHelp="* [topDir] is the top level directory (default\n  \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
queueHelp="* [queue] is the queue for running alignments (default \"$queue\")"
longQueueHelp="* [long queue] is the queue for running longer jobs such as the hic file\n  creation (default \"$long_queue\")"
siteHelp="* [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" \n  (default \"$site\")"
aboutHelp="* [about]: enter description of experiment, enclosed in single quotes"
shortHelp="* -r: use the short read version of the aligner, bwa aln\n  (default: long read, bwa mem)"
shortHelp2="* [end]: use the short read aligner on read end, must be one of 1 or 2 "
stageHelp="* [stage]: must be one of \"merge\", \"dedup\", \"final\", \"postproc\", or \"early\".\n    -Use \"merge\" when alignment has finished but the merged_sort file has not\n     yet been created.\n    -Use \"dedup\" when the files have been merged into merged_sort but\n     merged_nodups has not yet been created.\n    -Use \"final\" when the reads have been deduped into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the stats and\n     hic files"
pathHelp="* [chrom.sizes path]: enter path for chrom.sizes file"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
chunkHelp="* [chunk size]: number of lines in split files, must be multiple of 4\n  (default ${splitsize}, which equals $(awk -v ss=${splitsize} 'BEGIN{print ss/4000000}') million reads)"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
refSeqHelp="* [reference genome file]: enter path for reference sequence file, BWA index\n  files must be in same directory"
queueTimeHelp="* [queue time limit]: time limit for queue, i.e. -W 12:00 is 12 hours\n  (default ${queue_time})"
longQueueTimeHelp="* [long queue time limit]: time limit for long queue, i.e. -W 168:00 is one week\n  (default ${long_queue_time})"
excludeHelp="* -x: exclude fragment-delimited maps from hic file creation"
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$genomeHelp"
    echo -e "$dirHelp"
    echo -e "$queueHelp"
    echo -e "$longQueueHelp"
    echo -e "$siteHelp"
    echo -e "$aboutHelp"
    echo -e "$shortHelp"
    echo -e "$shortHelp2"
    echo -e "$stageHelp"
    echo -e "$pathHelp"
    echo -e "$siteFileHelp"
    echo -e "$refSeqHelp"
    echo -e "$chunkHelp"
    echo -e "$scriptDirHelp"
    echo -e "$queueTimeHelp"
    echo -e "$longQueueTimeHelp"
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}
while getopts "d:g:R:k:a:hrq:s:p:l:y:z:S:C:D:Q:L:x" opt; do
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
	C) splitsize=$OPTARG ;;
	D) juiceDir=$OPTARG ;;
	Q) queue_time=$OPTARG ;;
	L) long_queue_time=$OPTARG ;;
	x) nofrag=1 ;;
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
	dMel)      refSeq="${juiceDir}/references/Dmel_release5_first6.fasta";;
        mm9) refSeq="${juiceDir}/references/Mus_musculus_assembly9_norandom.fasta";;
	mm10)      refSeq="${juiceDir}/references/Mus_musculus_assembly10/v0/Mus_musculus_assembly10.fasta";;
	hg38)      refSeq="${juiceDir}/references/hg38/hg38.fa";;
        hg19) refSeq="${juiceDir}/references/Homo_sapiens_assembly19.fasta";;
	*)  echo "$usageHelp"
            echo "$genomeHelp"
            exit 1
    esac
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
    none) ligation="XXXX";;
    *)  ligation="XXXX"
	echo "$site not listed as recognized enzyme. Using $site_file as site file"
	echo "Ligation junction is undefined"
esac

## If DNAse-type experiment, no fragment maps
if [ "$site" == "none" ]
then
    nofrag=1;
fi

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
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [ ! -e "$site_file" ] && [ "$nofrag" -ne 1 ]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    exit 100
fi

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/HIC_tmp"

## Check that fastq directory exists and has proper fastq files
if [ ! -d "$topDir/fastq" ]; then
	echo "Directory \"$topDir/fastq\" does not exist."
	echo "Create \"$topDir/fastq\" and put fastq files to be aligned there."
	echo "Type \"juicer.sh -h\" for help"
	exit 100
else 
	if stat -t ${fastqdir} >/dev/null 2>&1
		then
		echo "(-: Looking for fastq files...fastq files exist"
	else
		if [ ! -d "$splitdir" ]; then 
			echo "***! Failed to find any files matching ${fastqdir}"
			echo "***! Type \"juicer.sh -h \" for help"
			exit 100		
		fi
	fi
fi

## Create output directory, only if not in merge, dedup, final, or postproc stages
if [[ -d "$outputdir" && -z "$final" && -z "$merge" && -z "$dedup" && -z "$postproc" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 100			
else
    if [[ -z "$final" && -z "$dedup" && -z "$merge" && -z "$postproc" ]]; then
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
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$dedup" ] && [ -z "$postproc" ]; then
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

#source $usePath
#$load_cluster


fastqsize=$(ls -lL  ${fastqdir} | awk '{sum+=$5}END{print sum}')
if [ "$fastqsize" -gt "2592410750" ]
then
    threads="-t 24"
fi
testname=$(ls -l ${fastqdir} | awk 'NR==1{print $9}')
if [ "${testname: -3}" == ".gz" ]
then
    skipsplit=1
    read1=${splitdir}"/*${read1str}*.fastq.gz"
else
    read1=${splitdir}"/*${read1str}*.fastq"
fi

# Add header containing command executed and timestamp:
jid=`sbatch <<- HEADER | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -p $queue
#SBATCH -t 2
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH -o $outDir/head-%j.out
#SBATCH -e $outDir/head-%j.err
#SBATCH -J "${groupname}_cmd"
date
echo "$0 $@"
echo "Juicer version:$juicer_version" 
HEADER`

## Not in merge, dedup, final, or postproc stage, i.e. need to split and align files.
if [ -z $merge ] && [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
	echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID with site file $site_file"

    ## Split fastq files into smaller portions for parallelizing alignment 
    ## Do this by creating a text script file for the job on STDIN and then 
    ## sending it to the cluster	
dependsplit="afterok"
	if [ ! $splitdirexists ]
	then
		echo "(-: Created $splitdir and $outputdir."
 
		if [ -n "$threads" ] && [ -z "$skipsplit" ] 
		then
      echo " Splitting files"
    for i in ${fastqdir}
    do
	filename=$(basename "$i")
        filename=${filename%.*}      
        jid=`sbatch <<- SPLITEND | egrep -o -e "\b[0-9]+$"
			#!/bin/bash -l
			#SBATCH -p $queue
			#SBATCH -t 1440
			#SBATCH -c 1
			#SBATCH --ntasks=1
			#SBATCH -o $outDir/split-%j.out
			#SBATCH -e $outDir/split-%j.err
			#SBATCH -J "${groupname}_split_${i}"
			date
			echo "Split file: $filename"
			split -a 3 -l $splitsize -d --additional-suffix=.fastq $i $splitdir/$filename
			date
SPLITEND`
        dependsplit="$dependsplit:$jid"
    done

	srun --ntasks=1 -c 1 -p "$queue" -t 1 -o $outDir/wait-%j.out -e $outDir/wait-%j.err -d $dependsplit -J "${groupname}_wait" sleep 1
        else
            cp -rs ${fastqdir} ${splitdir}
            wait
        fi
	else
		## No need to re-split fastqs if they already exist
		echo -e "---  Using already created files in $splitdir\n"
	fi

    ## Launch job. Once split/move is done, set the parameters for the launch. 
    
	echo "(-: Starting job to launch other jobs once splitting is complete"

    ## Loop over all read1 fastq files and create jobs for aligning read1,
    ## aligning read2, and merging the two. Keep track of merge names for final 
    ## merge. When merge jobs successfully finish, can launch final merge job.  


dependmerge="afterok"

for i in ${read1}
do
	ext=${i#*$read1str}
	name=${i%$read1str*} 
	# these names have to be right or it'll break
	name1=${name}${read1str}
	name2=${name}${read2str}	
	jname=$(basename "$name")${ext}
        usegzip=0
        if [ "${ext: -3}" == ".gz" ]
        then
            usegzip=1
	fi

	# count ligations
	sbatch <<- CNTLIG
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -t 1440
		#SBATCH -c 1
		#SBATCH --ntasks=1
		#SBATCH -o $outDir/count_ligation-%j.out
		#SBATCH -e $outDir/count_ligation-%j.err
		#SBATCH -J "${groupname}${jname}_Count_Ligation"
		date
		export ARG1=${usegzip}; export name=${name}; export name1=${name1}; export name2=${name2}; export ext=${ext}; export ligation=${ligation}; ${juiceDir}/scripts/countligations.sh
		date
CNTLIG

	# align read1 fastq
	if [ -n "$shortread" ] || [ "$shortreadend" -eq 1 ]
	then
		alloc_mem=16384
	else
		alloc_mem=12288
	fi	
	
	jid=`sbatch <<- ALGNR1 | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -o $outDir/align1-%j.out
		#SBATCH -e $outDir/align1-%j.err
		#SBATCH -t 1440
		#SBATCH -n 1
		#SBATCH -c 2
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=$alloc_mem
		#SBATCH -J "${groupname}_align1_${jname}"
		$load_bwa
		# Align read1
		date
		if [ -n "$shortread" ] || [ "$shortreadend" -eq 1 ]
		then		
            echo 'Running command bwa aln -q 15 $refSeq $name1$ext > $name1$ext.sai && bwa samse $refSeq $name1$ext.sai $name1$ext > $name1$ext.sam'
			srun --ntasks=1 bwa aln -q 15 $refSeq $name1$ext > $name1$ext.sai && srun --ntasks=1 bwa samse $refSeq $name1$ext.sai $name1$ext > $name1$ext.sam
			if [ \$? -ne 0 ]
			then 
				exit 100
			else
				echo "(-: Short align of $name1$ext.sam done successfully"
			fi		
		else
            echo 'Running command bwa mem $threads $refSeq $name1$ext > $name1$ext.sam '
			srun --ntasks=1 bwa mem -t 24 $refSeq $name1$ext > $name1$ext.sam
			if [ \$? -ne 0 ]
			then 
				exit 100
			else
				echo "(-: Mem align of $name1$ext.sam done successfully"
			fi		
		fi
		date
ALGNR1`

	dependalign="afterok:$jid"

	# align read2 fastq
	if [ -n "$shortread" ] || [ "$shortreadend" -eq 2 ]
	then
		alloc_mem=16384
	else
		alloc_mem=12288
	fi	
	
	jid=`sbatch <<- ALGNR2 | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -o $outDir/align2-%j.out
		#SBATCH -e $outDir/align2-%j.err
		#SBATCH -t 1440
		#SBATCH -n 1
		#SBATCH -c 2
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=$alloc_mem
		#SBATCH -J "${groupname}_align2_${jname}"
		$load_bwa
		date
		# Align read2
		if [ -n "$shortread" ] || [ "$shortreadend" -eq 2 ]
		then		
            echo 'Running command bwa aln -q 15 $refSeq $name2$ext > $name2$ext.sai && bwa samse $refSeq $name2$ext.sai $name2$ext > $name2$ext.sam '
			srun --ntasks=1 bwa aln -q 15 $refSeq $name2$ext > $name2$ext.sai && srun --ntasks=1 bwa samse $refSeq $name2$ext.sai $name2$ext > $name2$ext.sam
			if [ \$? -ne 0 ]
			then 
				exit 100
			else
				echo "(-: Short align of $name2$ext.sam done successfully"
			fi		
	
		else	
            echo 'Running command bwa mem $threads $refSeq $name2$ext > $name2$ext.sam'
			srun --ntasks=1 bwa mem -t 24 $refSeq $name2$ext > $name2$ext.sam
			if [ \$? -ne 0 ]
			then 
				exit 100
			else
				echo "(-: Mem align of $name2$ext.sam done successfully"
			fi		
		fi
		date		
ALGNR2`

	dependalign="$dependalign:$jid"

	# wait for top two, merge
	jid=`sbatch <<- MRGALL | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $long_queue
		#SBATCH -o $outDir/merge-%j.out
		#SBATCH -e $outDir/merge-%j.err
		#SBATCH --mem-per-cpu=2G
		#SBATCH -t 1440
		#SBATCH -c 8 
		#SBATCH --ntasks=1
		#SBATCH -d $dependalign
		#SBATCH -J "${groupname}_merge_${jname}"
		export LC_COLLATE=C
		date
		# sort read 1 aligned file by readname
		sort --parallel=8 -S 14G -T $tmpdir -k1,1 $name1$ext.sam > $name1${ext}_sort.sam
		#sort -T $tmpdir -k1,1 ${name1}${ext}.sam > ${name1}${ext}_sort.sam
		if [ \$? -ne 0 ]
		then 
			echo "***! Error while sorting $name1$ext.sam"
			exit 100
		else
			echo "(-: Sort read 1 aligned file by readname completed."
		fi
		
		# sort read 2 aligned file by readname 
		sort --parallel=8 -S 14G -T $tmpdir -k1,1 $name2$ext.sam > $name2${ext}_sort.sam
		#sort -T $tmpdir -k1,1 $name2$ext.sam > $name2${ext}_sort.sam
		if [ \$? -ne 0 ]
		then
		echo "***! Error while sorting $name2$ext.sam"
		exit 100
		else
			echo "(-: Sort read 2 aligned file by readname completed."
		fi
		
		# remove header, add read end indicator toreadname
		awk 'NF >= 11{\\\$1 = \\\$1"/1";print}' ${name1}${ext}_sort.sam > ${name1}${ext}_sort1.sam
		awk 'NF >= 11{\\\$1 = \\\$1"/2";print}' ${name2}${ext}_sort.sam > ${name2}${ext}_sort1.sam

		# merge the two sorted read end files
		sort --parallel=8 -S 14G -T $tmpdir -k1,1 -m $name1${ext}_sort1.sam $name2${ext}_sort1.sam > $name$ext.sam
		#sort -T $tmpdir -k1,1 -m $name1${ext}_sort1.sam $name2${ext}_sort1.sam > $name$ext.sam
		if [ \$? -ne 0 ]
		then
			echo "***! Failure during merge of read files"
			exit 100
		else
			rm $name1$ext.sa* $name2$ext.sa* $name1${ext}_sort*.sam $name2${ext}_sort*.sam
			echo "$name$ext.sam created successfully."
		fi

		# call chimeric_blacklist.awk to deal with chimeric reads; sorted file is sorted by read name at this point
		awk -v "fname1"=${name}${ext}_norm.txt -v "fname2"=${name}${ext}_abnorm.sam -v "fname3"=${name}${ext}_unmapped.sam -f $juiceDir/scripts/chimeric_blacklist.awk ${name}${ext}.sam

		# if any normal reads were written, find what fragment they correspond to and store that
		if [ -e "${name}${ext}_norm.txt" ] && [ "$site" != "none" ] 
		then
			$juiceDir/scripts/fragment.pl ${name}${ext}_norm.txt ${name}${ext}.frag.txt $site_file
			# sort by chromosome, fragment, strand, and position
			sort -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n ${name}${ext}.frag.txt > ${name}${ext}.sort.txt
			rm ${name}${ext}_norm.txt ${name}${ext}.frag.txt
		elif [ "$site" == "none" ] 
		then
		    awk '{printf("%s %s %s %d %s %s %s %d", \\\$1, \\\$2, \\\$3, 0, \\\$4, \\\$5, \\\$6, 1); for (i=7; i<=NF; i++) {printf(" %s",\\\$i);}printf("\n");}' $name${ext}_norm.txt > $name${ext}.frag.txt
		    # sort by chromosome, fragment, strand, and position
		    sort -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n ${name}${ext}.frag.txt > ${name}${ext}.sort.txt
		    rm ${name}${ext}.frag.txt
		else     
		    echo "***! No $name${ext}_norm.txt file created"    
		    exit 100
		fi
		date
		MRGALL`

	dependmerge="$dependmerge:$jid"
	# done looping over all fastq split files
done
fi

if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    if [ -z $merge ]
    then
        sbatch_wait="#SBATCH -d $dependmerge"
    else
        sbatch_wait=""
    fi

    # merge the sorted files into one giant file that is also sorted.
    jid=`sbatch <<- MRGSRT | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH --mem-per-cpu=16G
	#SBATCH -o $outDir/fragmerge-%j.out
	#SBATCH -e $outDir/fragmerge-%j.err
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_fragmerge"
	$sbatch_wait
	export LC_COLLATE=C
	date
	#Set to 3600 minutes, but should be set to ~4 minutes per split file.
	if ! sort --parallel=8 -S 120G -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt  > $outputdir/merged_sort.txt
	#if ! sort -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt  > $outputdir/merged_sort.txt
	then 
		echo "***! Some problems occurred somewhere in creating  sorted align files."
		exit 22
	else
		echo "(-: Finished sorting all sorted files into a single merge."
	fi
	date
MRGSRT`
dependmrgsrt="afterok:$jid"

fi

# Remove the duplicates from the big sorted file
if [ -z $final ] && [ -z $postproc ]
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
	#SBATCH -o $outDir/dedupguard-%j.out
	#SBATCH -e $outDir/dedupguard-%j.err
	#SBATCH -t 10
	#SBATCH -c 1
	#SBATCH -H
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_dedup_guard"
	${sbatch_wait}
	date
DEDUPGUARD`

	dependguard="afterok:$guardjid"

	# if jobs succeeded, kill the cleanup job, remove the duplicates from the big sorted file
    	jid=`sbatch <<- DEDUP | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $outDir/dedup-%j.out
	#SBATCH -e $outDir/dedup-%j.err
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_dedup"
	${sbatch_wait}
	date
	awk -v queue=$long_queue -v groupname=$groupname -v outDir=$outDir -v dir=$outputdir -v topDir=$topDir -v juicedir=$juiceDir -v site=$site -v genomeID=$genomeID -v genomePath=$genomePath -v user=$USER -v guardjid=$guardjid -f $juiceDir/scripts/split_rmdups.awk $outputdir/merged_sort.txt
	##Schedule new job to run after last dedup part:
	##Push guard to run after last dedup is completed:
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
	#SBATCH -o $outDir/post_dedup-%j.out
	#SBATCH -e $outDir/post_dedup-%j.err
	#SBATCH -t 100
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_post_dedup"
	#SBATCH -d ${dependguard}
	date
	rm -Rf $tmpdir;
	find $outDir -type f -size 0 | xargs rm
	date
MSPLITWAIT`

	dependmsplit="afterok:$jid"
fi

if [ -z "$genomePath" ]
then
        #If no path to genome is give, use genome ID as default.
        genomePath=$genomeID
fi

	#Skip if post-processing only is required
if [ -z $postproc ]
then

    if [ -z $final ]
    then
        sbatch_wait="#SBATCH -d $dependmsplit"
    else
        sbatch_wait=""
    fi
    
#Statistics of mapq>0:
jid=`sbatch <<- STATS0 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $outDir/stats0-%j.out
	#SBATCH -e $outDir/stats0-%j.err
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=6G
	#SBATCH -J "${groupname}_stats_0"
	${sbatch_wait}	
	$load_java
	date
	export IBM_JAVA_OPTIONS="-Xmx16384m -Xgcthreads1"

	echo -e 'Experiment description: $about' > $outputdir/inter.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/stats_dups.txt $outputdir/dups.txt
	cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt
	java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt

	date
STATS0`

	dependstats="afterok:$jid"
	
#Statistics for mapq>30:
jid=`sbatch <<- STATS30 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $outDir/stats30-%j.out
	#SBATCH -e $outDir/stats30-%j.err
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=6G
	#SBATCH -J "${groupname}_stats_30"
	${sbatch_wait}  
	$load_java
	date
	export IBM_JAVA_OPTIONS="-Xmx16384m -Xgcthreads1"

	echo -e 'Experiment description: $about' > $outputdir/inter_30.txt
	cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter_30.txt
	java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter_30.txt >> $outputdir/inter_30.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt
	
	date
STATS30`

	dependstats="$dependstats:$jid"	

#Generate abnormal files:
jid=`sbatch <<- ABNRML | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $outDir/abnormals-%j.out
	#SBATCH -e $outDir/abnormals-%j.err
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=6G
	#SBATCH -J "${groupname}_abnormals"
	${sbatch_wait}
	date

	cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam
	cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam
	awk -f ${juiceDir}/scripts/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt
	date
ABNRML`

	dependabnormals="afterok:$jid"	
	returnCode=$jid

# if early exit, we stop here, once the merged_nodups.txt file is created, and statistics were calculated.
if [ -z "$earlyexit" ]
then

	jid=`sbatch <<- HIC | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $long_queue
	#SBATCH -o $outDir/hic-%j.out
        #SBATCH -e $outDir/hic-%j.err	
	#SBATCH -t 1440
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH --mem-per-cpu=32G
	#SBATCH -J "${groupname}_hic"
	#SBATCH -d $dependstats
	${load_java}
	export IBM_JAVA_OPTIONS="-Xmx48192m -Xgcthreads1"
	date
	if [ -n "$nofrag" ]
	then 
	    ${juiceDir}/scripts/juicer_tools48g pre -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath 
	else 
	    ${juiceDir}/scripts/juicer_tools48g pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath
	fi
	date
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
	#SBATCH --mem-per-cpu=32G
	#SBATCH -J "${groupname}_hic30"
	#SBATCH -d ${dependstats}
	${load_java}
	export IBM_JAVA_OPTIONS="-Xmx48192m -Xgcthreads1"
	date
	if [ -n "$nofrag" ]
	then 
	    ${juiceDir}/scripts/juicer_tools48g pre -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath 
	else 
	    ${juiceDir}/scripts/juicer_tools48g pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath
	fi
	date
HIC30`

	dependhic30="${dependhic}:$jid"
fi

	if [ -z $postproc ]
	then
		sbatch_wait="#SBATCH -d $dependhic30"
	else
		sbatch_wait=""
	fi
        
	jid=`sbatch <<- HICCUPS | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=2G
	#SBATCH --gres=gpu:kepler:1
	#SBATCH -o $outDir/hiccups_wrap-%j.out
	#SBATCH -e $outDir/hiccups_wrap-%j.err
	#SBATCH -t 1440
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_hiccups_wrap"
	${sbatch_wait}
	${load_java}
	${load_gpu}
	date
	${juiceDir}/scripts/juicer_hiccups.sh -j ${juiceDir}/scripts/juicer_tools -i $outputdir/inter_30.hic -m ${juiceDir}/references/motif -g $genomeID
	date
HICCUPS`
        dependhiccups="afterok:$jid"

	jid=`sbatch <<- ARROWHEAD | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $outDir/arrowhead_wrap-%j.out
	#SBATCH -e $outDir/arrowhead_wrap-%j.err
	#SBATCH -t 1440
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_arrowhead_wrap"
	${sbatch_wait}
	${load_java}
	date
	${juiceDir}/scripts/juicer_arrowhead.sh -j ${juiceDir}/scripts/juicer_tools -i $outputdir/inter_30.hic
	date
ARROWHEAD`
        dependarrows="${dependhiccups}:$jid"

	jid=`sbatch <<- FINCLN1 | egrep -o -e "\b[0-9]+$"
	#!/bin/bash
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $outDir/fincln-%j.out
	#SBATCH -e $outDir/fincln-%j.err
	#SBATCH -t 1200
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_prep_done"
	#SBATCH -d $dependarrows
	date
	export splitdir=${splitdir}; export outputdir=${outputdir}; ${juiceDir}/scripts/check.sh
	${postjuicer}
	date
FINCLN1`
else

	jid=`sbatch <<- FINCLN1 | egrep -o -e "\b[0-9]+$" 
	#!/bin/bash
	#SBATCH -p $queue
	#SBATCH --mem-per-cpu=2G
	#SBATCH -o $outDir/fincln1-%j.out
	#SBATCH -e $outDir/fincln1-%j.err
	#SBATCH -t 1200
	#SBATCH -c 1
	#SBATCH --ntasks=1
	#SBATCH -J "${groupname}_prep_done"     
	#SBATCH -d $dependarrows
	date
	export splitdir=${splitdir}; export outputdir=${outputdir}; ${juiceDir}/scripts/check.sh
	${postjuicer}
	date
FINCLN1`
	dependprepdone="afterok:$jid"

fi
echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee.."
exit $returnCode
