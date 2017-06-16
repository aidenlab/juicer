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
juicer_version="1.5" 
## Set the following variables to work with your system

# set global tmpdir so no problems with /var/tmp
export TMPDIR=/broad/hptmp

## use cluster load commands:
usePath=/broad/software/scripts/useuse
load_bwa="use .bwa-0.7.12"
load_java="use Java-1.7"
load_cluster="use UGER"
load_coreutils="use Coreutils"
# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir="/broad/aidenlab"
# default queue, can also be set in options via -q
queue="short"
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
groupname="a"`date +%s`

## Default options, overridden by command line arguments

# top level directory, can also be set in options
topDir=`pwd`
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
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]\n                 [-a about] [-R end] [-S stage] [-p chrom.sizes path]\n                 [-y restriction site file] [-z reference genome file]\n                 [-C chunk size] [-D Juicer scripts directory]\n                 [-Q queue time limit] [-L long queue time limit] [-b ligation]\n                 [-r] [-h] [-x]"
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
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:R:k:a:hrq:s:p:l:y:z:S:C:D:Q:L:b:x" opt; do
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
	b) ligation=$OPTARG ;;
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
        mm9) refSeq="/seq/references/Mus_musculus_assembly9/v1/Mus_musculus_assembly9.fasta" ;;
        mm10) refSeq="/seq/references/Mus_musculus_assembly10/v0/Mus_musculus_assembly10.fasta" ;;
        hg38) refSeq="/seq/references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta" ;;
        hg19) refSeq="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta" ;;
	*)  echo "$usageHelp"
            echo "$genomeHelp"
            exit 1
    esac
else
    ## Reference sequence passed in, so genomePath must be set for the .hic file
    ## to be properly created
    if [ -z "$genomePath" ]
    then
        echo "***! WARNING: You must define a chrom.sizes file via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq";
        echo "***! .hic file may not be created properly";
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
    *) echo "$usageHelp"
	echo "$shortHelp2"
	exit 1
esac

if [ -z $site_file ] 
then
    site_file="${juiceDir}/restriction_sites/${genomeID}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [ ! -e "$site_file" ] && [ "$nofrag" -ne 1 ]
then
    echo "***! $site_file does not exist. It must be created before running this script."
    exit 1
fi

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
    echo "Create \"$topDir/$fastq\" and put fastq files to be aligned there."
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

## Create output directory, only if not in merge, dedup, final, or postproc stages
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
    mkdir $splitdir || { echo '***! Unable to create ${splitdir}, check permissions.' ; exit 1; }
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] 
then
    mkdir $tmpdir
    chmod 777 $tmpdir
fi

if [ ! -d "$debugdir" ]
then
    mkdir $debugdir
    chmod 777 $debugdir
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

source $usePath
$load_cluster
threads="4"
fastqsize=`ls -lL ${fastqdir} | awk '{sum+=$5}END{print sum}'`
if [ "$fastqsize" -gt "2592410750" ]
then
    splitme=1
fi

testname=`ls -l ${fastqdir} | awk 'NR==1{print $9}'`
if [ ${testname: -3} == ".gz" ]
then
    read1=${splitdir}"/*${read1str}*.fastq.gz"
    gzipped=1
else
    read1=${splitdir}"/*${read1str}*.fastq"
fi

myjid=$(qsub -terse -o ${debugdir}/version-${groupname}.out -j y -q $queue -r y -N ${groupname}cmd <<EOF
  date
  echo "Juicer version:$juicer_version" 
  echo "$0 $@"
EOF
)
echo "Version job $myjid" > ${debugdir}/jobs-${groupname}.out

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
        echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID. Ignoring fragment delimitation."
    fi
    ## Split fastq files into smaller portions for parallelizing alignment 
    ## Do this by creating a text script file for the job on STDIN and then 
    ## sending it to the cluster
    if [ ! $splitdirexists ]; then
        echo "(-: Created $splitdir and $outputdir.  Splitting files"
	
        if [ -n "$splitme" ]
        then
            for i in ${fastqdir}
            do
		filename=$(basename $i)
		filename=${filename%.*}      
                if [ -z "$gzipped" ]
                then
                    qsub -o ${debugdir}/split-${groupname}.out -e ${debugdir}/split-${groupname}.err -q $queue -r y -N ${groupname}split${countjobs} -sync y <<- SPLITEND
                    #!/bin/bash
                    source $usePath
                    $load_coreutils
                    echo "Split file: $filename"
                    split -a 3 -l $splitsize -d --additional-suffix=.fastq $i $splitdir/$filename
SPLITEND
		    wait
                else
                    qsub -o ${debugdir}/split-${groupname}.out -e ${debugdir}/split-${groupname}.err -q $queue -r y -N ${groupname}split${countjobs} -sync y <<- SPLITEND
                    #!/bin/bash
                    source $usePath
                    $load_coreutils
                    echo "Split file: $filename"
                    zcat $i | split -a 3 -l $splitsize -d --additional-suffix=.fastq - $splitdir/$filename
SPLITEND
		    wait
                    # if we split files, the splits are named .fastq
                    read1=${splitdir}"/*${read1str}*.fastq"
                fi
            done 
        else
            cp -rs ${fastqdir} ${splitdir}
            wait
        fi
    else
        ## No need to re-split fastqs if they already exist
        echo -e "---  Using already created files in $splitdir\n"
	# unzipped files have .fastq extension
	read1=${splitdir}"/*${read1str}*.fastq"
    fi

    ## Launch job. Once split/move is done, set the parameters for the launch. 
    
    echo "(-: Starting job to launch other jobs once splitting is complete..."

    ## Loop over all read1 fastq files and create jobs for aligning read1,
    ## aligning read2, and merging the two. Keep track of merge names for final 
    ## merge. When merge jobs successfully finish, can launch final merge job.  

    ## ARRAY holds the names of the jobs as they are submitted
    countjobs=0
    declare -a ARRAY
    declare -a JIDS
    declare -a TOUCH
    
    for i in ${read1}
    do
	ext=${i#*$read1str}
	name=${i%$read1str*} 
        # these names have to be right or it'll break
	name1=${name}${read1str}
	name2=${name}${read2str}	
        jname=$(basename $name)${ext}
        usegzip=0
        if [ ${ext: -3} == ".gz" ]
        then
            usegzip=1
        fi
        myjid=$(qsub -terse -o ${debugdir}/count_ligations-${groupname}-${jname}.out -e ${debugdir}/count_ligations-${groupname}-${jname}.err -q ${queue} -r y -N ${groupname}${jname}countligations -v usegzip=${usegzip},name=${name},name1=${name1},name2=${name2},ext=${ext},ligation=${ligation} ${juiceDir}/scripts/countligations.sh)
        echo "Count ligations $jname job $myjid" >> ${debugdir}/jobs-${groupname}.out
        # align read1 fastq
	if [ -n "$shortread" ] || [ "$shortreadend" -eq 1 ]
        then
 	    alloc_mem=16
        else
 	    alloc_mem=6
        fi	
#-pe smp 12
        touchfile1=${tmpdir}/${jname}1
        jid1=$(qsub -terse -o ${debugdir}/alignR1-${groupname}-${jname}.out -e ${debugdir}/alignR1-${groupname}-${jname}.err -q ${queue} -r y -N ${groupname}_align1${jname} -l h_vmem=${alloc_mem}g -pe smp ${threads} -binding linear:${threads}  <<- ALGNR1
        #!/bin/bash 
        source $usePath
        $load_bwa
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
                touch $touchfile1
  		echo "(-: Short align of $name1$ext.sam done successfully"
  	    fi		
			
        else
            echo 'Running command bwa mem -t $threads $refSeq $name1$ext > $name1$ext.sam '
            bwa mem -t $threads $refSeq $name1$ext > $name1$ext.sam
  	    if [ \$? -ne 0 ]
  	    then 
                echo "***! Error, failed to align $name1$ext" 
                touch $errorfile
  	        exit 1
  	    else
                touch $touchfile1
  		echo "(-: Mem align of $name1$ext.sam done successfully"
  	    fi		
        fi
ALGNR1
)
        echo "Align read1 $jname job $jid1" >> ${debugdir}/jobs-${groupname}.out
        touchfile2=${tmpdir}/${jname}2
 	# align read2 fastq
	if [ -n "$shortread" ] || [ "$shortreadend" -eq 2 ] 	
 	then
 	    alloc_mem=16
 	else
 	    alloc_mem=6
 	fi	
    
	jid2=$(qsub -terse -o ${debugdir}/alignR2-${groupname}-${jname}.out -e ${debugdir}/alignR2-${groupname}-${jname}.err -q ${queue} -r y -N ${groupname}_align2${jname} -l h_vmem=${alloc_mem}g -pe smp ${threads} -binding linear:${threads}  <<- ALGNR2
 	#!/bin/bash 
        source $usePath
        $load_bwa
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
                touch $touchfile2
 		echo "(-: Short align of $name2$ext.sam done successfully"
 	    fi				
 	else	
            echo 'Running command bwa mem -t $threads $refSeq $name2$ext > $name2$ext.sam'
 	    bwa mem -t $threads $refSeq $name2$ext > $name2$ext.sam
 	    if [ \$? -ne 0 ]
 	    then 
                echo "***! Error, failed to align $name2$ext" 
                touch $errorfile
 		exit 1
 	    else
                touch $touchfile2
 		echo "(-: Mem align of $name2$ext.sam done successfully"
 	    fi		
 	fi		
ALGNR2
);
        echo "Align read2 $jname job $jid2" >> ${debugdir}/jobs-${groupname}.out
        touchfile3=${tmpdir}/${jname}3
 	# wait for top two, merge
	jid3=$(qsub -terse  -o ${debugdir}/merge-${groupname}-${jname}.out -e ${debugdir}/merge-${groupname}-${jname}.err -q ${queue} -r y -N ${groupname}_merge${jname} -l h_vmem=8g -hold_jid ${groupname}_align1${jname},${groupname}_align2${jname}  <<- MRGALL
 	#!/bin/bash 
 	export LC_ALL=C
        source $usePath
	$load_coreutils
        if [ ! -f "${touchfile1}" ] || [ ! -f "${touchfile2}" ]
        then
            echo "***! Error, cluster did not finish aligning ${jname}"
            echo "Type the below to see what happened:"
            echo "qacct -j $jid1"
            echo "qacct -j $jid2"
            touch $errorfile 
            exit 1
        fi
 	# sort read 1 aligned file by readname 
 	sort -S 2G -T $tmpdir -k1,1f $name1$ext.sam > $name1${ext}_sort.sam
 	if [ \$? -ne 0 ]
 	then 
 	    echo "***! Error while sorting $name1$ext.sam"
            touch $errorfile
 	    exit 1
 	else
 	    echo "(-: Sort read 1 aligned file by readname completed."
 	fi
		
 	# sort read 2 aligned file by readname 
 	sort -S 2G -T $tmpdir -k1,1f $name2$ext.sam > $name2${ext}_sort.sam
 	if [ \$? -ne 0 ]
 	then
            echo "***! Error while sorting $name2$ext.sam"			 
            touch $errorfile
 	    exit 1
 	else
 	    echo "(-: Sort read 2 aligned file by readname completed."
 	fi
 	# remove header, add read end indicator to readname
 	awk 'NF >= 11{\$1 = \$1"/1";print}' $name1${ext}_sort.sam > $name1${ext}_sort1.sam
 	awk 'NF >= 11{\$1 = \$1"/2";print}' $name2${ext}_sort.sam > $name2${ext}_sort1.sam
		
 	# merge the two sorted read end files
 	sort -S 2G -T $tmpdir -k1,1f -m $name1${ext}_sort1.sam $name2${ext}_sort1.sam > $name$ext.sam
 	if [ \$? -ne 0 ]
 	then
            echo "***! Failure during merge of read files"			
            touch $errorfile
 	    exit 1
 	else
 	    rm $name1$ext.sa* $name2$ext.sa* $name1${ext}_sort*.sam $name2${ext}_sort*.sam
 	    echo "(-: $name$ext.sam created successfully."
            touch $touchfile3
 	fi 
MRGALL
)
        echo "Merge $jname job $jid3" >> ${debugdir}/jobs-${groupname}.out
        touchfile4=${tmpdir}/${jname}4
	jid4=$(qsub -terse  -o ${debugdir}/chimeric-${groupname}-${jname}.out -e ${debugdir}/chimeric-${groupname}-${jname}.err -q ${queue} -r y -hold_jid ${groupname}_merge${jname} -l h_vmem=4g -N ${groupname}_chimeric${jname}  <<- CHIMERIC
        #!/bin/bash                   
        export LC_ALL=C
        source $usePath
        $load_coreutils 
	if [ ! -f "${touchfile3}" ] 
        then
            echo "***! Error, cluster did not finish merging ${jname}"
            echo "Type qacct -j $jid3 to see what happened"
            touch $errorfile 
            exit 1
        fi
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
            awk '{printf("%s %s %s %d %s %s %s %d", \$1, \$2, \$3, 0, \$4, \$5, \$6, 1); for (i=7; i<=NF; i++) {printf(" %\s",\$i);}printf("\n");}' $name${ext}_norm.txt > $name${ext}.frag.txt
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
        sort -S 2G -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $name${ext}.frag.txt > $name${ext}.sort.txt
        if [ \$? -ne 0 ]        
        then
            echo "***! Failure during sort of $name${ext}"
            touch $errorfile
            exit 1
        else
            rm $name${ext}_norm.txt $name${ext}.frag.txt
            touch $touchfile4
        fi

CHIMERIC
)
        echo "Chimeric $jname job $jid4" >> ${debugdir}/jobs-${groupname}.out
        ARRAY[countjobs]="${groupname}_chimeric${jname}"
        JIDS[countjobs]="$jid4"
        TOUCH[countjobs]="$touchfile4"
        countjobs=$(( $countjobs + 1 ))
    done # done looping over all fastq split files

   # list of all jobs.  hold next steps until they have finished
    for (( i=0; i < countjobs; i++ ))
    do
        if [ $i -eq 0 ]; then
            holdjobs="-hold_jid ${ARRAY[i]}"
        else
            holdjobs="${holdjobs},${ARRAY[i]}"
        fi
    done    
    myjid=$(qsub -terse -o ${debugdir}/check-${groupname}.out -e ${debugdir}/check-${groupname}.err -q $queue -r y -N ${groupname}_check $holdjobs <<-CHECK
    # list of all jobs. print errors if failed
    for (( i=0; i < countjobs; i++ ))
    do
        echo "Checking ${TOUCH[i]}"
        if [ ! -f "${TOUCH[i]}" ]
        then
            echo "***! Error in ${ARRAY[i]}"
            echo "Type qacct -j ${JIDS[i]} to see what happened"        
        fi
    done
    for (( i=0; i < countjobs; i++ ))
    do
        if [ ! -f "${TOUCH[i]}" ]
        then
            touch $errorfile
        fi
    done

CHECK
)
    echo "Check alignment job $myjid" >> ${debugdir}/jobs-${groupname}.out
fi

if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    # merge the sorted files into one giant file that is also sorted.
    # note that if the $merge flag is set, this will run with $holdjobs empty
    # if a small file, do not bother parallelizing this sort
    if [ -n "$splitme" ] 
    then
        fragflags="-pe smp 8 -binding linear:8"
        pflags="--parallel=8"
    else
        fragflags=""
        pflags=""
    fi

    myjid=$(qsub -terse -o ${debugdir}/fragmerge-${groupname}.out -e ${debugdir}/fragmerge-${groupname}.err -q ${long_queue} -r y -N ${groupname}_fragmerge -l h_vmem=2g $fragflags $holdjobs <<- MRGSRT
    if [ -f "${errorfile}" ]
    then
        echo "Some alignment failed. Pipeline halting. Check ${debugdir} for results" | mail -s "Pipeline error" ${USER}@broadinstitute.org
        exit 1
    fi
    export LC_ALL=C
    source $usePath                  
    $load_coreutils 
    if [ -d $donesplitdir ]
    then
        mv $donesplitdir/* $splitdir/.
    fi

    if ! sort $pflags -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt  > $outputdir/merged_sort.txt
    then             
 	echo "***! Some problems occurred somewhere in creating sorted align files."
        echo "Failed to create merged_sort. Check ${debugdir} for results"| mail -s "Pipeline error" ${USER}@broadinstitute.org 
        exit 1
    else
    	echo "(-: Finished sorting all sorted files into a single merge."
    fi
MRGSRT
)
    echo "Fragmerge job $myjid" >> ${debugdir}/jobs-${groupname}.out
    if [ ! -z $merge ]
    then
        holdjobs="-hold_jid ${groupname}_fragmerge";
    else
        holdjobs="$holdjobs,${groupname}_fragmerge";
    fi
fi

touchfile5=${tmpdir}/${groupname}5
# Remove the duplicates from the big sorted file
if [ -z $final ] && [ -z $postproc ]
then
    myjid=$(qsub -terse -o ${debugdir}/dedup-${groupname}.out -e ${debugdir}/dedup-${groupname}.err -q long -r y -N ${groupname}_osplit $holdjobs <<- DEDUP
    if [ -f "${errorfile}" ]
    then
        exit 1
    fi
    source $usePath
    $load_cluster
    awk -v queue=$long_queue -v outfile=${debugdir}/dedup-${groupname}.out -v errfile=${debugdir}/dedup-${groupname}.err -v juicedir=${juiceDir} -v dir=$outputdir -v groupname=$groupname -f ${juiceDir}/scripts/split_rmdups.awk $outputdir/merged_sort.txt
    echo "(-: Finished launching dedup jobs"
    touch $touchfile5
DEDUP
)
    echo "Dedup job $myjid" >> ${debugdir}/jobs-${groupname}.out
    
    if [ ! -z $dedup ]
    then
        holdjobs="-hold_jid ${groupname}_osplit,${groupname}_rmsplit"; 
    else
        holdjobs="$holdjobs,${groupname}_osplit,${groupname}_rmsplit";
    fi
else
    touch $touchfile5
    touch ${debugdir}/dedup-${groupname}.err
fi

if [ -z "$genomePath" ]
then
    #If no path to genome is give, use genome ID as default.
    genomePath=$genomeID
fi
# if early exit, we stop here, once the merged_nodups.txt file is created.
if [ -z "$earlyexit" ]
then
    holdjobs3="-hold_jid ${groupname}_postproc"
    #Skip if post-processing only is required
    if [ -z $postproc ]
    then
	holdjobs1="-hold_jid ${groupname}_finallaunch"
	holdjobs2="-hold_jid ${groupname}_hic30"
	holdjobs3="-hold_jid ${groupname}_hic30,${groupname}_hic,${groupname}_stats,${groupname}_postproc"
	myjid=$(qsub -terse -o ${debugdir}/finallaunch-${groupname}.out -e ${debugdir}/finallaunch-${groupname}.err -q ${queue} -r y -N ${groupname}_finallaunch $holdjobs <<- DOSTAT
        #!/bin/bash 
        if [ -f "${errorfile}" ] 
        then
            exit 1
        fi
 	echo "(-: Alignment and merge done, launching other jobs."
        source $usePath
        $load_cluster
        qsub -o ${debugdir}/stats-${groupname}.out -e ${debugdir}/stats-${groupname}.err -q $long_queue -l h_vmem=16g -l m_mem_free=16g -N ${groupname}_stats -r y $holdjobs <<- EOF  
        if [ ! -f $touchfile5 ]; then touch $errorfile; exit 1; fi; if [ -s ${debugdir}/dedup-${groupname}.err ]; then touch $errorfile; exit 1; fi; source $usePath; $load_java; export _JAVA_OPTIONS=-Xmx16384m; export LC_ALL=en_US.UTF-8; echo -e 'Experiment description: $about' > $outputdir/inter.txt; ${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/stats_dups.txt $outputdir/dups.txt; cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt; java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt; ${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt; cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam; cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam; awk -f ${juiceDir}/scripts/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt
EOF
       qsub -o ${debugdir}/hic-${groupname}.out -e ${debugdir}/hic-${groupname}.err -q $long_queue -l m_mem_free=16g -l h_vmem=16g -N ${groupname}_hic -hold_jid ${groupname}_stats -r y <<- EOF 
        if [ ! -f $touchfile5 ]; then touch $errorfile; exit 1; fi; if [ -s ${debugdir}/dedup-${groupname}.err ]; then touch $errorfile; exit 1; fi; source $usePath; $load_java; export _JAVA_OPTIONS=-Xmx16384m; if [ "$nofrag" -eq 1 ]; then ${juiceDir}/scripts/juicebox pre -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath; else ${juiceDir}/scripts/juicebox pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath; fi ;
EOF
 	qsub  -o ${debugdir}/hic30-${groupname}.out -e ${debugdir}/hic30-${groupname}.err -q $long_queue -l h_vmem=16g -l m_mem_free=16g -N ${groupname}_hic30 -r y $holdjobs <<- EOF 
if [ ! -f $touchfile5 ]; then touch $errorfile; exit 1; fi; if [ -s ${debugdir}/dedup-${groupname}.err ]; then touch $errorfile; exit 1; fi; source $usePath; $load_java; export _JAVA_OPTIONS=-Xmx16384m; export LC_ALL=en_US.UTF-8; echo -e 'Experiment description: $about' > $outputdir/inter_30.txt; cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter_30.txt; java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter_30.txt >> $outputdir/inter_30.txt; ${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt; if [ "$nofrag" -eq 1 ]; then ${juiceDir}/scripts/juicebox pre -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath; else ${juiceDir}/scripts/juicebox pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath; fi
EOF
DOSTAT
)
	echo "Statistics job $myjid" >> ${debugdir}/jobs-${groupname}.out
    fi
    myjid=$(qsub -terse  -o ${debugdir}/postproc-${groupname}.out -e ${debugdir}/postproc-${groupname}.err -q ${queue} -r y -N ${groupname}_postprocessing $holdjobs1 -cwd <<- POSTPROC
    if [ -f "${errorfile}" ]
    then
        exit 1
    fi
    source $usePath
    $load_cluster 
    qsub -o ${debugdir}/postproc-${groupname}.out -e ${debugdir}/postproc-${groupname}.err -q $long_queue -l h_vmem=16g -l m_mem_free=16g -N ${groupname}_postproc $holdjobs2 -cwd <<-EOF
    if [ -f "${errorfile}" ]
    then
       exit 1
    fi
    source $usePath; 
    $load_java; 
    export _JAVA_OPTIONS=-Xmx16384m; 
    export LC_ALL=en_US.UTF-8; 
    ${juiceDir}/scripts/juicer_postprocessing.sh -j ${juiceDir}/scripts/juicebox -i $outputdir/inter_30.hic -m ${juiceDir}/references/motif -g $genomeID
    rm -r ${tmpdir}
EOF
POSTPROC
)
    echo "Postprocessing job $myjid" >> ${debugdir}/jobs-${groupname}.out
    myjid=$(qsub -terse -o ${debugdir}/finalcheck-${groupname}.out -e ${debugdir}/finalcheck-${groupname}.err -q ${queue} -r y -N ${groupname}_done $holdjobs1 <<- FINCLN2
    source $usePath
    $load_cluster
    qsub  -o ${debugdir}/finalcheck-${groupname}.out -e ${debugdir}/finalcheck-${groupname}.err -q $queue -N ${groupname}_prep_finalize -v splitdir=${splitdir},outputdir=${outputdir} $holdjobs3 ${juiceDir}/scripts/check.sh
FINCLN2
)
    echo "Final check job $myjid" >> ${debugdir}/jobs-${groupname}.out
else
qsub -o ${debugdir}/finalcheck-${groupname}.out -e ${debugdir}/finalcheck-${groupname}.err -q ${queue} -r y -N ${groupname}_finallaunch $holdjobs <<- CHECKEARLY
    source $usePath
    $load_cluster
    qsub  -o ${debugdir}/finalcheck-${groupname}.out -e ${debugdir}/finalcheck-${groupname}.err -q $queue -N ${groupname}_done -v splitdir=${splitdir},outputdir=${outputdir},early=1 $holdjobs ${juiceDir}/scripts/check.sh
CHECKEARLY
fi
echo "(-: Finished adding all jobs... please wait while processing."
