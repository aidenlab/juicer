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
# Juicer version 1.5
shopt -s extglob
juicer_version="1.5.6" 
## Set the following variables to work with your system

# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir="/opt/juicer"
# default queue, can also be set in options via -q
queue="normal"
# default queue time, can also be set in options via -Q
queue_time="1200"
# default long queue, can also be set in options via -l
long_queue="long"
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
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]\n                 [-a about] [-R end] [-S stage] [-p chrom.sizes path]\n                 [-y restriction site file] [-z reference genome file]\n                 [-C chunk size] [-D Juicer scripts directory]\n                 [-Q queue time limit] [-L long queue time limit] [-b ligation] [-t threads]\n                 [-r] [-h] [-x]"
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
ligationHelp="* [ligation junction]: use this string when counting ligation junctions"
threadsHelp="* [threads]: number of threads when running BWA alignment"
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
    echo -e "$ligationHelp"
    echo -e "$threadsHelp"
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:R:a:hrq:s:p:l:y:z:S:C:D:Q:L:b:t:x" opt; do
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
	C) splitsize=$OPTARG; splitme=1 ;;
	D) juiceDir=$OPTARG ;;
	Q) queue_time=$OPTARG ;;
	L) long_queue_time=$OPTARG ;;
	x) nofrag=1 ;;
	b) ligation=$OPTARG ;;
	t) threads=$OPTARG ;;
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
	mm9) refSeq="${juiceDir}/references/Mus_musculus_assembly9_norandom.fasta";;
	mm10) refSeq="${juiceDir}/references/Mus_musculus_assembly10.fasta";;
	hg38) refSeq="${juiceDir}/references/hg38.fa";;
	hg19) refSeq="${juiceDir}/references/Homo_sapiens_assembly19.fasta";;
	*)  echo "$usageHelp"
            echo "$genomeHelp"
            exit 1
    esac
else
    # Reference sequence passed in, so genomePath must be set for the .hic file
    # to be properly created
    if [ -z "$genomePath" ]
        then
        echo "***! You must define a chrom.sizes file via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq";
        exit 1;
    fi
fi

## Check that refSeq exists 
if [ ! -e "$refSeq" ]; then
    echo "***! Reference sequence $refSeq does not exist";
    exit 1;
fi

## Check that index for refSeq exists
if [ ! -e "${refSeq}.bwt" ]; then
    echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run bwa index on this file before running juicer.";
    exit 1;
fi

## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]
then
    case $site in
        HindIII) ligation="AAGCTAGCTT";;
        DpnII) ligation="GATCGATC";;
        MboI) ligation="GATCGATC";;
        NcoI) ligation="CCATGCATGG";;
        none) ligation="XXXX";;
        *)  ligation="XXXX"
	    echo "$site not listed as recognized enzyme. Using $site_file as site file"
	    echo "Ligation junction is undefined";;
    esac
fi

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
	exit 1
esac

if [ -z "$site_file" ]
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [ ! -e "$site_file" ] && [ "$nofrag" -ne 1 ]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    exit 1
fi

## Set threads for sending appropriate parameters to cluster and string for BWA call
if [ -z "$threads" ]
then
    threads="$(getconf _NPROCESSORS_ONLN)"
fi

threadstring="-t $threads"
#alloc_mem=$(($threads * 5000))
alloc_mem=12000

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/HIC_tmp"
debugdir=${topDir}"/debug"

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

## Create output directory, only if not in dedup, final, or postproc stages
if [[ -d "$outputdir" && -z "$final" && -z "$dedup" && -z "$postproc" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 1			
else
    if [[ -z "$final" && -z "$dedup" && -z "$postproc" ]]; then
        mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 1; } 
    fi
fi

## Create split directory
if [ -d "$splitdir" ]; then
    splitdirexists=1
else
    mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 1; }
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] 
then
    mkdir $tmpdir
    chmod 777 $tmpdir
fi

## Create debug directory, used for reporting commands output
if [ ! -d "$debugdir" ]
then
    mkdir $debugdir
    chmod 777 $debugdir
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

# If chunk size sent in, split. Otherwise check size before splitting
if [ -z $splitme ]
then
    fastqsize=`ls -lL ${fastqdir} | awk '{sum+=$5}END{print sum}'`
    if [ "$fastqsize" -gt "2592410750" ]
    then
	splitme=1
    fi
fi

testname=`ls -l ${fastqdir} | awk 'NR==1{print $9}'`
if [ ${testname: -3} == ".gz" ]
then
    read1=${splitdir}"/*${read1str}*.fastq.gz"
    gzipped=1
else
    read1=${splitdir}"/*${read1str}*.fastq"
fi

bsub -o ${debugdir}/head-${groupname}.out -q "$queue" -W "$queue_time" -J "${groupname}_cmd" <<EOF
	#!/bin/bash
	date
	# Experiment description
	if [ -n "${about}" ]
	then
		printf 'Experiment description: ${about}; '
	else
		printf 'Experiment description: '
	fi

	# Get version numbers of all software
	printf "Juicer version $juicer_version;"
	bwa 2>&1 | awk '\$1=="Version:"{printf(" BWA %s; ", \$2)}' 
	printf "$threads threads; "
	if [ -n "$splitme" ]
	then
		printf "splitsize $splitsize; "
	fi  
	java -version 2>&1 | awk 'NR==1{printf("%s; ", \$0);}' 
	${juiceDir}/scripts/juicer_tools -V 2>&1 | awk '\$1=="Juicer" && \$2=="Tools"{printf("%s; ", \$0);}' 

  	echo "$0 $@"
EOF
headfile="${debugdir}/head-${groupname}.out"

## Record if we failed while aligning, so we don't waste time on other jobs
## Remove file if we're relaunching Juicer
errorfile=${debugdir}/${groupname}_alignfail
if [ -f $errorfile ]
then
    rm $errorfile
fi

## Not in merge, dedup, final, or postproc stage, i.e. need to split and align files.
if [ -z $merge ] && [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    if [ "$nofrag" -eq 0 ]
    then
        echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID with site file $site_file"
    else
        echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID with no fragment delimited maps."
    fi
    
    ## Split fastq files into smaller portions for parallelizing alignment 
    ## Do this by creating a text script file for the job on STDIN and then 
    ## sending it to the cluster
    if [ ! $splitdirexists ]; then
        echo "(-: Created $splitdir and $outputdir.  Splitting files"
	
        if [ -n "$splitme" ]
        then
	    countjobs=0
	    declare -a ARRAY
            for i in ${fastqdir}
            do
		filename=$(basename $i)
		filename=${filename%.*}      
		if [ -z "$gzipped" ]
                then
		    bsub << SPLITEND

			#!/bin/bash
			#BSUB -q Split
			#BSUB -W $queue_time
			#BSUB -K
			#BSUB -o ${debugdir}/split-${groupname}.out 
			#BSUB -e ${debugdir}/split-${groupname}.err
			#BSUB -J "${groupname}_split_${i}"
	        	#echo "Split file: $filename"
                	#Following line is used in coreutils >= 8.16, if not, simpler version is used below.
	        	#split -a 3 -l $splitsize -d --additional-suffix=.fastq ${i} $splitdir/$filename
                	split -a 3 -l $splitsize -d ${i} $splitdir/$filename
SPLITEND
		    wait
				
		    ARRAY[countjobs]="${groupname}_split_${i}"
		    countjobs=$(( countjobs + 1 ))
		else
		    bsub << SPLITEND
			#!/bin/bash
			#BSUB -q Split
			#BSUB -W $queue_time
			#BSUB -K
			#BSUB -o ${debugdir}/split-${groupname}.out 
			#BSUB -e ${debugdir}/split-${groupname}.err
			#BSUB -J "${groupname}_split_${i}"
	        	#echo "Split file: $filename"
                	#Following line is used in coreutils >= 8.16, if not, simpler version is used below.
	        	#split -a 3 -l $splitsize -d --additional-suffix=.fastq ${i} $splitdir/$filename
			zcat $i | split -a 3 -l $splitsize -d - $splitdir/$filename
SPLITEND
		    wait
                    # if we split files, the splits are named .fastq
                    read1=${splitdir}"/*${read1str}*.fastq"
                fi
	    done


            # Once split succeeds, rename the splits as fastq files
	    ## LSF users change queue below to $queue 
	    bsub << SPLITMV
            #!/bin/bash
            #BSUB -q Split
            #BSUB -W $queue_time
            #BSUB -K
            #BSUB -o ${debugdir}/splitmv-${groupname}.out 
	    #BSUB -e ${debugdir}/splitmv-${groupname}.err
            #BSUB -J "${groupname}_move"
            for j in $splitdir/*;do mv \${j} \${j}.fastq;done
SPLITMV

	    ARRAY[countjobs]="${groupname}_move"
	    countjobs=$(( $countjobs + 1 ))

            ## List of split jobs to wait for

	    for (( i=0; i < countjobs; i++ ))
	    do
		if [ $i -eq 0 ]; then
		    exitjobs="exit(${ARRAY[i]}) "
		else
		    exitjobs="$exitjobs || exit(${ARRAY[i]})"
		fi
	    done
	else # we're not splitting but just copying
	    cp -rs ${fastqdir} ${splitdir}
            wait
        fi
    else
        ## No need to re-split fastqs if they already exist
        echo -e "---  Using already created files in $splitdir\n"
        # unzipped files will have .fastq extension, softlinked gz 
        testname=$(ls -l ${splitdir} | awk '$9~/fastq$/||$9~/gz$/{print $9; exit}')
        if [ ${testname: -3} == ".gz" ]
        then
            read1=${splitdir}"/*${read1str}*.fastq.gz"
        else
	    read1=${splitdir}"/*${read1str}*.fastq"
        fi
    fi

    ## Launch job. Once split/move is done, set the parameters for the launch. 
    
    echo "(-: Starting job to launch other jobs once splitting is complete"
    
    ## Loop over all read1 fastq files and create jobs for aligning read1,
    ## aligning read2, and merging the two. Keep track of merge names for final
    ## merge. When merge jobs successfully finish, can launch final merge job.
    ## ARRAY holds the names of the jobs as they are submitted
    countjobs=0
    declare -a ARRAY

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
        ## LSF users change queue below to $queue
	bsub <<-CNTLIG
	#!/bin/bash
	#BSUB -q Align
        #BSUB -W $queue_time
	#BSUB -o ${debugdir}/count_ligations-${groupname}-${jname}.out
	#BSUB -e ${debugdir}/count_ligations-${groupname}-${jname}.err
	#BSUB -J "${groupname}${jname}_Count_Ligation"
	export usegzip=${usegzip}; export name=${name}; export name1=${name1}; export name2=${name2}; export ext=${ext}; export ligation=${ligation}; ${juiceDir}/scripts/countligations.sh
CNTLIG

	## LSF users change queue below to $queue
	bsub <<- ALGNR1
	#!/bin/bash
	#BSUB -q Align
	#BSUB -W $queue_time
	#BSUB -o ${debugdir}/alignR1-${groupname}-${jname}.out
	#BSUB -e ${debugdir}/alignR1-${groupname}-${jname}.err
	#BSUB -J "${groupname}_align1${jname}"
	#BSUB -R "rusage[mem=${alloc_mem}]"

	# Align read1 
	if [ -n "$shortread" ] || [ "$shortreadend" -eq 1 ]
	then		
           echo 'Running command bwa aln -q 15 $refSeq $name1$ext > $name1$ext.sai && bwa samse $refSeq $name1$ext.sai $name1$ext > $name1$ext.sam'
	   bwa aln -q 15 $refSeq $name1$ext > $name1$ext.sai && bwa samse $refSeq $name1$ext.sai $name1$ext > $name1$ext.sam
	   if [ \$? -ne 0 ]
	   then 
              echo "***! Error, failed to align $name1$ext" 
              touch $errorfile
  	      exit 1
	   else
	      echo "(-: Short align of $name1$ext.sam done successfully"
	   fi		
			
	else	
           echo 'Running command bwa mem $threadstring $refSeq $name1$ext > $name1$ext.sam '
	   bwa mem $threadstring $refSeq $name1$ext > $name1$ext.sam
	   if [ \$? -ne 0 ]
	   then
	      echo "***! Error, failed to align $name1$ext" 
              touch $errorfile
  	      exit 1 
	   else
	      echo "(-: Mem align of $name1$ext.sam done successfully"
	   fi		
	fi
ALGNR1

	bsub <<- ALGNR2
	#!/bin/bash
	#BSUB -q Align
	#BSUB -o ${debugdir}/alignR2-${groupname}-${jname}.out
	#BSUB -e ${debugdir}/alignR2-${groupname}-${jname}.err
        #BSUB -W $queue_time	
	#BSUB -R "rusage[mem=${alloc_mem}]"
	#BSUB -J "${groupname}_align2${jname}"
	# Align read2
	if [ -n "$shortread" ] || [ "$shortreadend" -eq 2 ]
	then		
           echo 'Running command bwa aln -q 15 $refSeq $name2$ext > $name2$ext.sai && bwa samse $refSeq $name2$ext.sai $name2$ext > $name2$ext.sam '
	   bwa aln -q 15 $refSeq $name2$ext > $name2$ext.sai && bwa samse $refSeq $name2$ext.sai $name2$ext > $name2$ext.sam
	   if [ \$? -ne 0 ]
	   then 
              echo "***! Error, failed to align $name2$ext" 
              touch $errorfile
 	      exit 1
	   else
	      echo "(-: Short align of $name2$ext.sam done successfully"
	   fi		
			
	else	
          echo 'Running command bwa mem $threadstring $refSeq $name2$ext > $name2$ext.sam'
	  bwa mem $threadstring $refSeq $name2$ext > $name2$ext.sam
	  if [ \$? -ne 0 ]
	  then 
	      echo "***! Error, failed to align $name2$ext" 
              touch $errorfile
 	      exit 1
	  else
	     echo "(-: Mem align of $name2$ext.sam done successfully"
	  fi		
       fi		
ALGNR2

	# wait for top two, merge
	bsub <<- MRGALL
	#!/bin/bash
	#BSUB -q Align
        #BSUB -W $long_queue_time
	#BSUB -o ${debugdir}/merge-${groupname}-${jname}.out 
	#BSUB -e ${debugdir}/merge-${groupname}-${jname}.err
	#BSUB -R "rusage[mem=8000]"
	#BSUB -w " done(${groupname}_align1${jname}) && done(${groupname}_align2${jname}) "
	#BSUB -J "${groupname}_merge${jname}"
	export LC_ALL=C
		
	# sort read 1 aligned file by readname 
	sort -T $tmpdir -k1,1 $name1$ext.sam > $name1${ext}_sort.sam
	if [ \$? -ne 0 ]
	then 
	   echo "***! Error while sorting $name1$ext.sam"
	   echo "Sort of $name1$ext.sam failed. Check ${debugdir} for results"
	   exit 1
	else
	   echo "(-: Sort read 1 aligned file by readname completed."
	fi
		
	# sort read 2 aligned file by readname 
	sort -T $tmpdir -k1,1 $name2$ext.sam > $name2${ext}_sort.sam
	if [ \$? -ne 0 ]
	then 
	   echo "***! Error while sorting $name2$ext.sam"  
	   echo "Sort of $name2$ext.sam failed. Check ${debugdir} for results"
    	   exit 1
	else
	   echo "(-: Sort read 2 aligned file by readname completed."
	fi
		
	# remove header, add read end indicator toreadname
	awk 'NF >= 11{\$1 = \$1"/1";print}' $name1${ext}_sort.sam > $name1${ext}_sort1.sam
	awk 'NF >= 11{\$1 = \$1"/2";print}' $name2${ext}_sort.sam > $name2${ext}_sort1.sam	
	# merge the two sorted read end files
	sort -T $tmpdir -k1,1 -m $name1${ext}_sort1.sam $name2${ext}_sort1.sam > $name$ext.sam
	if [ \$? -ne 0 ]
	then
	   echo "***! Failure during merge of read files"
           echo "Merge of $name$ext.sam failed. Check ${debugdir} for results"
	   exit 1
	else
	   rm $name1$ext.sa* $name2$ext.sa* $name1${ext}_sort*.sam $name2${ext}_sort*.sam
	   echo "$name$ext.sam created successfully."
	fi 
MRGALL
       		
	bsub <<- CHIMERIC
	#!/bin/bash
	#BSUB -q Align
        #BSUB -W $queue_time
	#BSUB -o ${debugdir}/chimeric-${groupname}-${jname}.out 
	#BSUB -e ${debugdir}/chimeric-${groupname}-${jname}.err 
	#BSUB -R "rusage[mem=8000]"
	#BSUB -w " done(${groupname}_merge${jname})"
	#BSUB -J "${groupname}_chimeric${jname}"
	export LC_ALL=C	
	
	# so there are no complaints later if empty
        touch $name${ext}_abnorm.sam 
	# call chimeric_blacklist.awk to deal with chimeric reads; sorted file is sorted by read name at this point
	awk -v "fname1"=$name${ext}_norm.txt -v "fname2"=$name${ext}_abnorm.sam -v "fname3"=$name${ext}_unmapped.sam -f ${juiceDir}/scripts/chimeric_blacklist.awk $name$ext.sam
		
        if [ \$? -ne 0 ]                                              
        then                                       
	    echo "***! Failure during chimera handling of $name${ext}"
            touch $errorfile
            exit 1  
        fi
  
	# if any normal reads were written, find what fragment they correspond to and store that
        # check if site file exists and if so write the fragment number even if nofrag set
        # one is not obligated to provide a site file if nofrag set; but if one does, frag
        # numbers will be calculated correctly
        if [ -e "$name${ext}_norm.txt" ] && [ "$site" != "none" ] && [ -e "$site_file" ]
        then
            ${juiceDir}/scripts/fragment.pl ${name}${ext}_norm.txt ${name}${ext}.frag.txt $site_file
        elif [ "$site" == "none" ] || [ "$nofrag" -eq 1 ]
        then
            awk '{printf("%s %s %s %d %s %s %s %d", \$1, \$2, \$3, 0, \$4, \$5, \$6, 1); for (i=7; i<=NF; i++) {printf(" %s",\$i);}printf("\n");}' $name${ext}_norm.txt > $name${ext}.frag.txt 
        else 
           echo "***! No $name${ext}_norm.txt file created"
           touch $errorfile
           exit 1
        fi
	if [ \$? -ne 0 ] 
        then 
           echo "***! Failure during fragment assignment of $name${ext}"
           touch $errorfile
           exit 1 
        fi
	# sort by chromosome, fragment, strand, and position
	sort -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $name${ext}.frag.txt > $name${ext}.sort.txt
        if [ \$? -ne 0 ]
        then  
          echo "***! Failure during sort of $name${ext}"
          touch $errorfile
          exit 1                                                        
        else
	  rm $name${ext}_norm.txt $name${ext}.frag.txt
	fi	
CHIMERIC

	# list of all jobs. if any fail, i.e., exit(jobid) != 0, we will kill
	# the remaining jobs
	ARRAY[countjobs]="exit(${groupname}_align1${jname}) || exit(${groupname}_align2${jname}) || exit(${groupname}_merge${jname}) || exit(${groupname}_chimeric${jname}) "
	countjobs=$(( countjobs + 1 ))
	# done looping over all fastq split files		
    done
    
    for (( i=0; i < countjobs; i++ ))
    do
	# clean up jobs if any fail
        ## LSF users change queue below to $queue
        ## LSF users change bkill to bkill -g ${groupname} 0
	bsub <<- CLNFAIL
	#!/bin/bash
	#BSUB -q Clean2
        #BSUB -W $queue_time
	#BSUB -o ${debugdir}/clean${i}-${groupname}.out
	#BSUB -e ${debugdir}/clean${i}-${groupname}.err
	#BSUB -w " ${ARRAY[i]} "
	#BSUB -J "cleanup_${groupname}_${i}"
	bkill -q Align 0
	echo 
CLNFAIL
    done

    # kill the kill jobs if everything went well
    ## LSF users change queue below to $queue
    ## LSF users change bkill to kill cleanup jobs
    bsub <<- FNLKILL
    #!/bin/bash
    #BSUB -q Merge
    #BSUB -W $queue_time
    #BSUB -o ${debugdir}/cleankill-${groupname}.out
    #BSUB -e ${debugdir}/cleankill-${groupname}.err
    #BSUB -w " done(${groupname}_chimeric*) " 
    #BSUB -J "${groupname}_fragmerge1"
    bkill -q Clean2 0
    echo
FNLKILL
fi

if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    # merge the sorted files into one giant file that is also sorted.
    ## LSF users change queue below to $long_queue
    bsub <<- MRGSRT
    #!/bin/bash
    #BSUB -q Merge
    #BSUB -W $long_queue_time
    #BSUB -o ${debugdir}/fragmerge-${groupname}.out 
    #BSUB -e ${debugdir}/fragmerge-${groupname}.err
    #BSUB -w " done(${groupname}_chimeric*) "
    #BSUB -J "${groupname}_fragmerge"
    export LC_ALL=C
    if [ -d $donesplitdir ]
    then
       mv $donesplitdir/* $splitdir/.
    fi
    if ! sort -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt  > $outputdir/merged_sort.txt
    then 
	echo "***! Some problems occurred somewhere in creating  sorted align files."
    else
	echo "(-: Finished sorting all sorted files into a single merge."
        rm -r ${tmpdir}		
    fi
MRGSRT
    ## LSF users change queue below to $queue
    ## LSF users change bkill to kill groupname
    bsub <<- CLNDIE
    #!/bin/bash
    #BSUB -q Clean3
    #BSUB -W $queue_time
    #BSUB -o ${debugdir}/cleanfragmerge-${groupname}.out 
    #BSUB -e ${debugdir}/cleanfragmerge-${groupname}.err
    #BSUB -J "${groupname}_clean1"
    #BSUB -w " exit(${groupname}_fragmerge) "
    bkill -q Merge 0
    echo

CLNDIE
fi
		
# if jobs succeeded, kill the cleanup job, 
# Remove the duplicates from the big sorted file
if [ -z $dedup ] && [ -z $final ] && [ -z $postproc ]
  then
    waitstring="#BSUB -w \" done(${groupname}_fragmerge*) \""
fi

if [ -z $final ] && [ -z $postproc ]
then
    ## LSF users change queue below to $queue
    ## LSF users change bkill to kill cleanup job
    bsub <<- KILLCLNUP
    #!/bin/bash
    #BSUB -q Clean3
    #BSUB -W $queue_time
    #BSUB -o ${debugdir}/dedup-${groupname}.out 
    #BSUB -e ${debugdir}/dedup-${groupname}.err  
    ${waitstring}
    #BSUB -J "${groupname}_osplit"
    bkill -J ${groupname}_clean1
    awk -v queue=$long_queue -v outfile=${debugdir}/dedup-${groupname}.out -v juicedir=${juiceDir}  -v dir=$outputdir -v queuetime=$long_queue_time -v groupname=$groupname -f ${juiceDir}/scripts/split_rmdups.awk $outputdir/merged_sort.txt
KILLCLNUP

    # if it dies, cleanup and write to relaunch script
    ## LSF users change queue below to $queue
    ## LSF users change bkill to kill groupname
    bsub <<- DIECLNRLNCH
    #!/bin/bash
    #BSUB -q Stat
    #BSUB -W $queue_time
    #BSUB -o ${debugdir}/cleandedup-${groupname}.out 
    #BSUB -e ${debugdir}/cleandedup-${groupname}.err
    #BSUB -J "${groupname}_clean2"
    #BSUB -w " exit(${groupname}_osplit) "
    bkill -q Clean3 0
    echo
DIECLNRLNCH
fi

if [ -z "$genomePath" ]
then
    #If no path to genome is give, use genome ID as default.
    genomePath=$genomeID
fi

# if early exit, we stop here, once the merged_nodups.txt file is created.
if [ -z "$earlyexit" ]
then
    diefinal=" -w \"exit(${groupname}_postproc)\""
    waitstring5="-w \"done(${groupname}_postproc)\""
    #Skip if post-processing only is required
    if [ -z $postproc ]
    then
	
	if [ -z $final ]
	then
	    waitstring2="#BSUB -w \" done(${groupname}_osplit) \""
	    execstring="bkill -J ${groupname}_clean2 "
	    waitstring22="-w \"done(${groupname}_rmsplit) && done(${groupname}_osplit)\""
	fi
        ## LSF users change queue below to $queue
	bsub <<- DOSTAT
        #!/bin/bash
        #BSUB -q Stat
        #BSUB -W $queue_time
        #BSUB -o ${debugdir}/stats-${groupname}.out 
	#BSUB -e ${debugdir}/stats-${groupname}.err
        ${waitstring2}
        #BSUB -J "${groupname}_launch"
        echo "(-: Alignment and merge done, launching other jobs."
        $execstring
        bsub -o ${debugdir}/stats-${groupname}.out -e ${debugdir}/stats-${groupname}.err -q $long_queue -W $long_queue_time -R "rusage[mem=16000]" $waitstring22 -J "${groupname}_stats" "df -h;_JAVA_OPTIONS=-Xmx16384m; export LC_ALL=en_US.UTF-8; tail -n1 ${headfile} | awk '{printf\"%-1000s\n\", \\\$0}' > $outputdir/inter.txt; ${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/stats_dups.txt $outputdir/dups.txt; cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt; java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt; ${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt; cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam; cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam; awk -f ${juiceDir}/scripts/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt"
        bsub -o ${debugdir}/hic-${groupname}.out -e ${debugdir}/hic-${groupname}.err -q $long_queue -W $long_queue_time -R "rusage[mem=16000]" -w "done(${groupname}_stats)" -J "${groupname}_hic" "df -h;export _JAVA_OPTIONS=-Xmx16384m; if [ \"$nofrag\" -eq 1 ]; then ${juiceDir}/scripts/juicer_tools pre -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath; else ${juiceDir}/scripts/juicer_tools pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath; fi ;"
        bsub -o ${debugdir}/stats30hic30-${groupname}.out -e ${debugdir}/stats30hic30-${groupname}.err -q $long_queue -W $long_queue_time -R "rusage[mem=16000]" $waitstring22 -J "${groupname}_hic30" "df -h;export _JAVA_OPTIONS=-Xmx16384m; export LC_ALL=en_US.UTF-8; tail -n1 $headfile | awk '{printf\"%-1000s\n\", \\\$0}' > $outputdir/inter_30.txt; cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter_30.txt; java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter_30.txt >> $outputdir/inter_30.txt; ${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt; export _JAVA_OPTIONS=-Xmx8192m; if [ \"$nofrag\" -eq 1 ]; then ${juiceDir}/scripts/juicer_tools pre -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath; else ${juiceDir}/scripts/juicer_tools pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath; fi"
	echo
DOSTAT
	waitstring3="#BSUB -w \" done(${groupname}_launch) \""
	waitstring4="-w \"done(${groupname}_hic30)\""
	waitstring5="-w \"done(${groupname}_stats) && done(${groupname}_postproc) && done(${groupname}_hic30)  && done(${groupname}_hic)\""
	diefinal=" -w \"exit(${groupname}_postproc) || exit(${groupname}_stats) || exit(${groupname}_hic) || exit(${groupname}_hic30)\""
    fi
    ## LSF users change queue below to $queue
    bsub <<- POSTPROC
    #!/bin/bash
    #BSUB -q Stat
    #BSUB -W $queue_time
    #BSUB -o ${debugdir}/postproc-${groupname}.out 
    #BSUB -e ${debugdir}/postproc-${groupname}.err 
    ${waitstring3}
    #BSUB -J "${groupname}_postproc_wrap"
    bsub -o ${debugdir}/postproc-${groupname}.out -e ${debugdir}/postproc-${groupname}.err  -q $long_queue -W $long_queue_time -R "rusage[mem=12000]" $waitstring4 -J "${groupname}_postproc" "${juiceDir}/scripts/juicer_postprocessing.sh -j ${juiceDir}/scripts/juicer_tools -i ${outputdir}/inter_30.hic -m ${juiceDir}/references/motif -g ${genomeID}"
    echo
POSTPROC
    ## LSF users change queue below to $queue
    bsub <<- FINCLN1
    #!/bin/bash
    #BSUB -q CleanFnl
    #BSUB -W $queue_time
    #BSUB -o ${debugdir}/finalcheck-${groupname}.out 
    #BSUB -e ${debugdir}/finalcheck-${groupname}.err
    #BSUB -J "${groupname}_prep_done"
    #BSUB -w " done(${groupname}_postproc_wrap) "
    bsub -o ${debugdir}/finalcheck-${groupname}.out -e ${debugdir}/finalcheck-${groupname}.err -q $queue -W $queue_time $waitstring5 -J "${groupname}_done" "bkill -J ${groupname}_clean3; export splitdir=${splitdir}; export outputdir=${outputdir}; ${juiceDir}/scripts/check.sh;"
    echo
FINCLN1
    # if it dies, cleanup 
    ## LSF users change queue below to $queue
    ## LSF users change bkill to kill groupname
    bsub <<- DIEFINAL
    #!/bin/bash
    #BSUB -q Stat
    #BSUB -W $queue_time
    #BSUB -o ${debugdir}/finalclean-${groupname}.out 
    #BSUB -e ${debugdir}/finalclean-${groupname}.err
    $waitstring3
    bsub -o ${debugdir}/finalclean-${groupname}.out -e ${debugdir}/finalclean-${groupname}.err -q $queue -W $queue_time $diefinal -J "${groupname}_clean3" "bkill -q $long_queue 0; bkill -q CleanFnl 0; bkill -q $queue 0;"
    echo
DIEFINAL
else
    ## LSF users change queue below to $queue
    bsub <<- FINCLN1
    #!/bin/bash
    #BSUB -q CleanFnl
    #BSUB -o ${debugdir}/finalclean-${groupname}.out 
    #BSUB -e ${debugdir}/finalclean-${groupname}.err
    #BSUB -J "${groupname}_prep_done"
    #BSUB -w " done(${groupname}_launch) "
    bsub -o ${debugdir}/finalclean-${groupname}.out -e ${debugdir}/finalclean-${groupname}.err -q $queue -W $queue_time -w "done(${groupname}_clean2)" -J "${groupname}_done" "export splitdir=${splitdir}; export outputdir=${outputdir}; ${juiceDir}/scripts/check.sh"
    echo
FINCLN1
fi
echo "(-: Finished adding all jobs... please wait while processing."
