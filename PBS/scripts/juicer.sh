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
juicer_version="1.5.6" 

# timestamp the pipeline and report the juicer app version, and shell
date +"%Y-%m-%d %H:%M:%S"
echo "Juicer version:$juicer_version"
echo "$0 $@"
    
## Set the following variables to work with your system
# path additionals, make sure paths are correct for your system
#myPath=/bin:$PATH
# set global tmpdir so no problems with /var/tmp
## use cluster load commands:
#usePath=""
load_bwa="module load bwa/0.7.15"
load_java='module load java/jdk1.8.0_131'
#load_cluster=""
#load_coreutils=""
load_cuda='module load cuda/7.5.18/gcc/4.4.7'

# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir=""
# default queue, can also be set in options via -q
queue="batch"
# default queue time, can also be set in options via -Q
walltime="walltime=24:00:00"
# default long queue, can also be set in options via -l
long_queue="batch"
# default long queue time, can also be set in options via -L
long_walltime="walltime=120:00:00"
# size to split fastqs. adjust to match your needs. 4000000=1M reads per split
# can also be changed via the -C flag
# give your email address to be used in #PBS -M to receive notifications when job error occurs.
# Must be either set with an email address or skipped
# This email is not included in the launch stat and postprocessing steps, add manually if needed
EMAIL='#PBS -M xxx@gmail.com'
splitsize=90000000
# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1"
read2str="_R2"


# unique groupname for jobs submitted in each run. lab initial with an timestamp
# Length of $groupname in this PBS version needs to be no longer than 8 characters.
# Groupname Length limitatioin is critical, because jobID will be obtained through qstat output using jobName.
# The string being searched from qstat output is "jobstr_${groupname}", where "jobstr" is 
# the job specific name string in the qsub -N option. max length of "jobstr" in this PBS version script is 7.
# our PBS cluster has an maximum jobName column width of 16.
# Change the $groupname max length accordingly based on your PBS cluster qstat output setup.
groupname="C$(date "+%s"|cut -c 6-11)"
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
usageHelp="Usage: ${0##*/} [-W group_list=genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]\n                 [-a about] [-R end] [-S stage] [-p chrom.sizes path]\n                 [-y restriction site file] [-z reference genome file]\n                 [-C chunk size] [-D Juicer scripts directory]\n                 [-Q queue time limit] [-L long queue time limit] [-r] [-h] [-x]"
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
queueTimeHelp="* [queue time limit]: time limit for queue, i.e. -l 12:00 is 12 hours\n  (default ${walltime})"
longQueueTimeHelp="* [long queue time limit]: time limit for long queue, i.e. -l 168:00 is one week\n  (default ${long_walltime})"
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
    Q) walltime=$OPTARG ;;
    L) long_walltime=$OPTARG ;;
    x) nofrag=1 ;; #no fragment maps
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
case $site in
    HindIII) ligation="AAGCTAGCTT";;
    DpnII) ligation="GATCGATC";;
    MboI) ligation="GATCGATC";;
    NcoI) ligation="CCATGCATGG";;
    none) ligation="XXXX";;
    *)  ligation="XXXX"
    echo "$site not listed as recognized enzyme. Using $site_file as site file"
    echo "Ligation junction is undefined"
    exit 1
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
    *)    echo "$usageHelp"
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

## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq"
outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/HIC_tmp"
logdir=${topDir}"/logs"
## Check that fastq directory exists and has proper fastq files
if [ ! -d "${topDir}/fastq" ]; then
    echo "Directory \"${topDir}/fastq\" does not exist."
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

## Create output directory, only if not in merge, dedup, final, or postproc stages
if [[ -d "$outputdir" && -z "$final" && -z "$merge" && -z "$dedup" && -z "$postproc" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 1            
else
    if [[ -z "$final" && -z "$dedup" && -z "$merge" && -z "$postproc" ]]; then
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
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$dedup" ] && [ -z "$postproc" ]; then
    mkdir "$tmpdir"
    chmod 777 "$tmpdir"
fi
 
## Create a log directory, used to store log files (standout and standerr of each submitted jobs)
if [ ! -d "$logdir" ] && [ -z "$final" ] && [ -z "$dedup" ] && [ -z "$postproc" ]; then
    mkdir "$logdir"
    chmod 777 "$tmpdir"
fi

## Arguments have been checked and directories created. Now begins the real work of the pipeline
#source $usePath
#$load_cluster

## jobIDstring holds the jobIDs as they are submitted
countjobs=0
jobIDstring=""
#check the total fastq input files size, determine the threads number
#threads value could be customized according to user's job and cluster resources
fastqsize=$(ls -lL ${fastqdir} | awk '{sum+=$5}END{print sum}')
if [ "$fastqsize" -gt "2592410750" ]
then
    threads=16
fi

testname=$(ls -l ${fastqdir} | awk 'NR==1;{print $9}')
if [ "${testname: -3}" == ".gz" ]
then
    skipsplit=1
    read1=${splitdir}"/*${read1str}*.fastq.gz"
else
    read1=${splitdir}"/*${read1str}*.fastq"
fi



## Not in merge, dedup, final, or postproc stage, i.e. need to split and align files.
if [[ -z "$final" && -z "$merge" && -z "$dedup" && -z "$postproc" ]] 
then
    echo -e "(-: Aligning files matching $fastqdir\n in queue $queue to genome $genomeID with site file $site_file"

    ## Split fastq files into smaller portions for parallelizing alignment 
    ## Do this by creating a text script file for the job on STDIN and then 
    ## sending it to the cluster
    if [ ! $splitdirexists ]
    then
        echo "(-: Created $splitdir and $outputdir."
        if [ -n "$threads" ] && [ -z "$skipsplit" ] 
        then
            echo " Splitting files" 
            jIDs_spit=""
            for i in ${fastqdir}
            do
                filename=$(basename "$i")
                filename=${filename%.*}

                #wait till the previous qsub job has finished sucessfully
                #submitted job might get delayed due to time in the queue.
                timestamp=$(date +"%s" | cut -c 4-10)
                jID_split=$(qsub <<SPLITEND
                #PBS -S /bin/bash
                #PBS -q $queue
                #PBS -l $walltime
                #PBS -l nodes=1:ppn=1:AMD
                #PBS -l mem=20gb
                ${EMAIL}
                #PBS -m a
                #PBS -o ${logdir}/${timestamp}_split_${filename}_${groupname}.log
                #PBS -j oe
                #PBS -N split_${filename}_${groupname}

                date +"%Y-%m-%d %H:%M:%S"
                echo "Split file: $filename"
                #Following line is used in coreutils >= 8.16, if not, simpler version is used below.
                #split -a 3 -l $splitsize -d --additional-suffix=.fastq ${i} $splitdir/$filename
                split -a 3 -l $splitsize -d ${i} ${splitdir}/${filename}
SPLITEND
)
                jIDs_split="${jIDs_split}:${jID_split}"
           done

            ## wait till the fastq spliting job has sucessfully finished
            ## Once split succeeds, rename the splits as fastq files
            ## PBS users change queue below to $queue 
            timestamp=$(date +"%s" | cut -c 4-10)
            jID_splitmv=$(qsub << SPLITMV
            #PBS -S /bin/bash
            #PBS -q $queue
            #PBS -l $walltime
            #PBS -l nodes=1:ppn=1:AMD
            #PBS -l mem=20gb
            ${EMAIL}
            #PBS -m a
            #PBS -o ${logdir}/${timestamp}_move_${groupname}.log
            #PBS -j oe
            #PBS -N move_${groupname}
            #PBS -W depend=afterok${jIDs_split}
            date +"%Y-%m-%d %H:%M:%S"
            for j in ${splitdir}/*;do mv \${j} \${j}.fastq;done
SPLITMV
)
        
        else # we're not splitting but just copying
            cp -rs ${fastqdir} ${splitdir}
            wait
        fi
    else # split dir already exists, no need to re-split fastqs 
        echo -e "---  Using already created files in $splitdir\n"
    fi

    ## Launch job. Once split/move is done, set the parameters for the launch. 
    
    echo "(-: Starting job to launch other jobs once splitting is complete"
    
    # Loop over all read1 fastq files and create jobs for aligning read1,
    # aligning read2, and merging the two. Keep track of merge names for final 
    # merge. When merge jobs successfully finish, can launch final merge job.  

    if [ -n "$threads" ] && [ -z "$skipsplit" ]
    then
        waitstring_alnwrp="#PBS -W depend=afterok:${jID_splitmv}"
    fi
    echo "waitstring_alnwrp is ${waitstring_alnwrp}"

    timestamp=$(date +"%s" | cut -c 4-10)
    qsub <<ALIGNWRAP
    #PBS -S /bin/bash
    #PBS -q $queue
    #PBS -l $walltime
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=6gb
    #PBS -o ${logdir}/${timestamp}_alnwrap_${groupname}.log
    #PBS -j oe
    #PBS -N AlnWrp${groupname}
    $waitstring_alnwrp
    ${EMAIL}
    #PBS -m a
    for i in ${read1}
    do

        countjobs=\$(( countjobs + 1 ))
        ext=\${i#*$read1str}
        name=\${i%$read1str*} 
        # these names have to be right or it'll break
        name1=\${name}${read1str}
        name2=\${name}${read2str}    
        jname=\$(basename "\$name")\${ext}
        usegzip=0
        if [ "\${ext}" == ".gz" ]
        then
            usegzip=1
        fi

        ## count ligations
        timestamp=\$(date +"%s" | cut -c 4-10)
        qsub <<-CNTLIG
        #PBS -S /bin/bash
        #PBS -q $queue
        #PBS -l $walltime
        #PBS -l nodes=1:ppn=1:AMD
        #PBS -l mem=4gb
        ${EMAIL}
        #PBS -m a
        #PBS -o ${logdir}/\${timestamp}_\${jname}_CntLig_\${countjobs}_${groupname}.log
        #PBS -j oe
        #PBS -N CtLig\${countjobs}${groupname}
        #PBS -v name=\${name}
        #PBS -v name1=\${name1}
        #PBS -v name2=\${name2}
        #PBS -v ext=\${ext}
        #PBS -v ligation=${ligation}
        #PBS -v usezip=\${usezip}

        date +"%Y-%m-%d %H:%M:%S"
        export usegzip=\${usegzip}
        export name=\${name}
        export name1=\${name1}
        export name2=\${name2}
        export ext=\${ext}
        export ligation=${ligation}
        ${juiceDir}/scripts/countligations.sh
CNTLIG
        jID_cntlig=\$(qstat | grep "CtLig\${countjobs}${groupname}" | cut -d ' ' -f 1 )
        echo "jID_cntlig \${countjobs} id is \${jID_cntlig}"
        ## Align read1        
        # align read1 fastq
        # allocate memory for alignment according to threads number
        # set-up the max memory according to the cluster set-up
        # this local cluster uses shared mem per node when in resource quest
        alloc_mem_value=\$(( $threads * 5 ))
        alloc_mem=\${alloc_mem_value}"gb"
        #set the max mem cap accordingly to ensure the job resource request meet the cluster-setup
        if [ \$alloc_mem_value -gt 240 ]
        then
            alloc_mem=240gb
        fi
 
        echo "starting read1 alignemnt"
        
        timestamp=\$(date +"%s" | cut -c 4-10)
        qsub <<ALGNR1
        #PBS -S /bin/bash
        #PBS -q $queue
        #PBS -l $walltime
        #PBS -l nodes=1:ppn=${threads}:AMD
        #PBS -l mem=\${alloc_mem}
        ${EMAIL}
        #PBS -m a
        #PBS -o ${logdir}/\${timestamp}_\${jname}_align1_\${countjobs}_${groupname}.log
        #PBS -j oe
        #PBS -N ALN1\${countjobs}${groupname}
        #PBS -W depend=afterok:\$jID_cntlig
        #PBS -v name=\${name}
        #PBS -v name1=\${name1}
        #PBS -v ext=\${ext}
     
        date +"%Y-%m-%d %H:%M:%S"
        $load_bwa
        if [ -n "$shortread" ] || [ "$shortreadend" -eq 1 ]
        then        
            echo 'Running command bwa aln -q 15 $refSeq \${name1}\${ext} > \${name1}\${ext}.sai && bwa samse $refSeq \${name1}\${ext}.sai \${name1}\${ext} > \${name1}\${ext}.sam'
            bwa aln -q 15 $refSeq \${name1}\${ext} > \${name1}\${ext}.sai && bwa samse $refSeq \${name1}\${ext}.sai \${name1}\${ext} > \${name1}\${ext}.sam
            if [ \$? -ne 0 ]
            then 
                echo "Alignment of \${name1}\${ext} failed. Check ${topDir}/pbs.out for results"
                exit 1
            else
                echo " Short align of \${name1}\${ext}.sam done successfully"
            fi
        else    
            echo 'Running command bwa mem $threads $refSeq \${name1}\${ext} > \${name1}\${ext}.sam '
            bwa mem -t $threads $refSeq \${name1}\${ext} > \${name1}\${ext}.sam
            if [ \$? -ne 0 ]
            then
                echo "***!Error, failed to align \${name1}\${ext}"
                exit 1
            else
                echo "Mem align of \${name1}\${ext}.sam done successfully"
            fi
        fi
        echo "below is the number of lines in align1 .sam file"
        cat \${name1}\${ext}.sam |wc -l
ALGNR1
        wait
        # Get the jobID from qstat ouput by searching job specific string,"read1\${countjobs}" , in job Name.
        jID_1=\$( qstat | grep "ALN1\${countjobs}${groupname}" |cut -d ' ' -f 1 )
        echo "align1 \$countjobs id is: \${jID_1}"
        # Align read2
        echo "starting read2 alignment" 
        # align read2 fastq    
        timestamp=\$(date +"%s" | cut -c 4-10)
        qsub <<ALGNR2
        #PBS -S /bin/bash
        #PBS -q $queue
        #PBS -l $walltime
        #PBS -l nodes=1:ppn=${threads}:AMD
        #PBS -l mem=\$alloc_mem
        ${EMAIL}
        #PBS -m a
        #PBS -o ${logdir}/\${timestamp}_\${jname}_align2_\${countjobs}_${groupname}.log
        #PBS -j oe
        #PBS -N ALN2\${countjobs}${groupname}
        #PBS -W depend=afterok:\$jID_cntlig
        #PBS -v name=\${name}
        #PBS -v name2=\${name2}
        #PBS -v ext=\${ext}
    
        date +"%Y-%m-%d %H:%M:%S"
        $load_bwa
        if [ -n "$shortread" ] || [ "$shortreadend" -eq 2 ]
        then        
            echo 'Running command bwa aln -q 15 $refSeq \${name2}\${ext} > \${name2}\${ext}.sai && bwa samse $refSeq \${name2}\${ext}.sai \${name2}\${ext} > \${name2}\${ext}.sam '
            bwa aln -q 15 $refSeq \${name2}\${ext} > \${name2}\${ext}.sai && bwa samse $refSeq \${name2}\${ext}.sai \${name2}\${ext} > \${name2}\${ext}.sam
            if [ \$? -ne 0 ]
            then 
                echo "Alignment of \${name2}\${ext} failed. Check ${topDir}/pbs.out for results"
                exit 1
            else
                echo "Short align of \${name2}\${ext}.sam done successfully"
            fi
            
        else    
            echo 'Running command bwa mem $threads $refSeq \${name2}\${ext} > \${name2}\${ext}.sam'
            bwa mem -t $threads $refSeq \${name2}\${ext} > \${name2}\${ext}.sam
            if [ \$? -ne 0 ]
            then 
                exit 1
            else
                echo "Mem align of \${name2}\${ext}.sam done successfully"
            fi        
        fi
        echo "below is the number of lines in aligned .sam file before sorting"
        cat \${name2}\${ext}.sam |wc -l        
ALGNR2
        wait

        # Get the jobID of job above from qstat ouput using job specific string,"read2\${countjobs}" , in jobName.
        jID_2=\$(qstat | grep "ALN2\${countjobs}${groupname}" |cut -d ' ' -f 1 )
        echo "align2 \${countjobs} id is: \${jID_2}"
        echo "starting merging from read1 and read2"
        # wait for align1 and align2 jobs finish,then merge
        timestamp=\$(date +"%s" | cut -c 4-10)
        qsub <<- MRGALL
        #PBS -S /bin/bash
        #PBS -q $queue
        #PBS -l $long_walltime
        #PBS -l nodes=1:ppn=1:AMD
        #PBS -l mem=24gb
        ${EMAIL}
        #PBS -m a
        #PBS -o ${logdir}/\${timestamp}_\${jname}_merge_\${countjobs}_${groupname}.log
        #PBS -j oe
        #PBS -N Mrg\${countjobs}${groupname}
        #PBS -W depend=afterok:\${jID_1}:\${jID_2}
        #PBS -v name=\${name}
        #PBS -v name1=\${name1}
        #PBS -v name2=\${name2}
        #PBS -v ext=\${ext}
        #PBS -v countjobs=\${countjobs}

        date +"%Y-%m-%d %H:%M:%S"
        export LC_ALL=C

        # sort read 1 aligned file by readname 
        sort -T $tmpdir -k1,1f \${name1}\${ext}.sam > \${name1}\${ext}_sort.sam
        if [ \$? -ne 0 ]
        then 
            echo "***! Error while sorting \${name1}\${ext}.sam"
            echo "Sort of \${name1}\${ext}.sam failed."
            exit 1
        else
            echo "Sort read 1 aligned file by readname completed."
        fi
        echo "below is the number of lines in sorted read1 .sam files"
        cat \${name1}\${ext}_sort.sam | wc -l

        # sort read 2 aligned file by readname 
        sort -T $tmpdir -k1,1f \${name2}\${ext}.sam > \${name2}\${ext}_sort.sam
        if [ \$? -ne 0 ]
        then 
            echo "***! Error while sorting \${name2}\${ext}.sam"  
            echo "Sort of \${name2}\${ext}.sam failed."
            exit 1
        else
            echo "Sort read 2 aligned file by readname completed."
        fi
        echo "below is the number of lines in sorted read2 .sam files"
        cat \${name2}\${ext}_sort.sam | wc -l

        # remove header, add read end indicator to read name
        awk -f ${juiceDir}/scripts/read1_sortproc.awk \${name1}\${ext}_sort.sam > \${name1}\${ext}_sort1.sam
        awk -f ${juiceDir}/scripts/read2_sortproc.awk \${name2}\${ext}_sort.sam > \${name2}\${ext}_sort1.sam    
        echo "below is the number of lines in \${name1}\${ext}_sort1.sam"
        cat \${name1}\${ext}_sort1.sam | wc -l
        echo "below is the number of lines in \${name2}\${ext}_sort1.sam"
        cat \${name2}\${ext}_sort1.sam | wc -l

        # merge the two sorted read end files
        sort -T $tmpdir -k1,1f -m \${name1}\${ext}_sort1.sam \${name2}\${ext}_sort1.sam > \${name}\${ext}.sam
        if [ $? -ne 0 ]
        then
            echo "***! Failure during merge of read files"
            echo "Merge of \${name}\${ext}.sam failed"
            exit 1
        else
            rm \${name1}\${ext}.sa* \${name2}\${ext}.sa* \${name1}\${ext}_sort*.sam \${name2}\${ext}_sort*.sam
            echo "\${name}\$next.sam created successfully."
        fi 
        echo "below is the number of lines in \${name}\${ext}_sort.sam after merging sorted aln1 and aln2"
        cat \${name}\${ext}_sort.sam | wc -l

MRGALL
    wait

    jID_3=\$(qstat | grep "Mrg\${countjobs}${groupname}" | cut -d ' ' -f 1 )
    echo "merging align1 and align2 \${coutjobs} id is \${jID_3}"
    echo "starting chimeric step after alignment"
    timestamp=\$(date +"%s" | cut -c 4-10)
    qsub <<- CHIMERIC
    #PBS -S /bin/bash
    #PBS -q $queue
    #PBS -l $walltime
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=24gb
    ${EMAIL}
    #PBS -m a
    #PBS -o ${logdir}/\${timestamp}_\${jname}_chimeric_\${countjobs}_${groupname}.log
    #PBS -j oe
    #PBS -N Chmr\${countjobs}${groupname}
    #PBS -W depend=afterok:\${jID_3}
    #PBS -v name=\${name}
    #PBS -v ext=\${ext}
 
    date +"%Y-%m-%d %H:%M:%S"
    export LC_ALL=C    
    # call chimeric_blacklist.awk to deal with chimeric reads; sorted file is sorted by read name at this point
    awk -v fname1="\${name}\${ext}_norm.txt" -v fname2="\${name}\${ext}_abnorm.sam" -v fname3="\${name}\${ext}_unmapped.sam" -f ${juiceDir}/scripts/chimeric_blacklist.awk \${name}\${ext}.sam

    if [ \$? -ne 0 ]                                              
    then                                       
        echo "***! Failure during chimera handling of \${name}\${ext}"
        echo "Chimera handling of \${name}\${ext}.sam failed."
        exit 1                                                     
    fi  
    # if any normal reads were written, find what fragment they correspond to and store that
    if [ -e "\${name}\${ext}_norm.txt" ] && [ "$site" != "none" ] 
    then  
        ${juiceDir}/scripts/fragment.pl \${name}\${ext}_norm.txt \${name}\${ext}.frag.txt $site_file
    elif [ "$site" == "none" ]
    then
        awk -f ${juiceDir}/scripts/chimeric_nonsites.awk \${name}\${ext}_norm.txt > \${name}\${ext}.frag.txt 
    else 
        echo "***! No \${name}\${ext}_norm.txt file created"
        echo "Creation of \${name}\${ext}_norm.txt failed."
        exit 1
    fi
    if [ \$? -ne 0 ] 
    then 
        echo "***! Failure during fragment assignment of \${name}\${ext}"
        echo "Fragment assignment of \${name}\${ext}.sam failed."
        exit 1 
    fi
    # sort by chromosome, fragment, strand, and position
    sort -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n \${name}\${ext}.frag.txt > \${name}\${ext}.sort.txt
    if [ \$? -ne 0 ]
    then  
        echo "***! Failure during sort of \${name}\${ext}"
        echo "Sort of \${name}\${ext}.frag.txt failed."
        exit 1                                                       
    else
        echo "removing temperary files \${name}\${ext}_norm.txt \${name}\${ext}.frag.txt"
        rm \${name}\${ext}_norm.txt \${name}\${ext}.frag.txt
    fi
CHIMERIC
    wait

    jID_4=\$(qstat | grep "Chmr\${countjobs}${groupname}" | cut -d ' ' -f 1)
    echo "chimeric \$countjobs id is \$jID_4"
    exitstatus=\$(qstat -f \${jID_4} |grep "exit_status" )
    echo "the exit status of \{jID_4} is \${exitstatus}"
    jobIDstring="\${jobIDstring}:\${jID_4}"
    echo "jobIDstring \$countjobs is \${jobIDstring}"

    # done looping over all fastq split files
    done
 
   
    # if error occored, we will kill the remaining jobs
    # output an error message of error detection and killing the remaining jobs
    timestamp=\$(date +"%s" | cut -c 4-10)
    qsub <<- CKALIGNFAIL
    #PBS -S /bin/bash
    #PBS -q $queue  
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=2gb
    #PBS -l $walltime
    ${EMAIL}
    #PBS -m a
    #PBS -o ${logdir}/\${timestamp}_check_alnOK_${groupname}.log
    #PBS -j oe
    #PBS -W depend=afterok\${jobIDstring}
    #PBS -N AlnOK_${groupname}
    
    date +"%Y-%m-%d %H:%M:%S"
    echo "Sucess: All alignment jobs were successfully finished without failure!"
CKALIGNFAIL

    timestamp=\$(date +"%s" | cut -c 4-10)
    qsub <<- CKALIGNFAILCLN
    #PBS -S /bin/bash
    #PBS -q $queue  
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=4gb
    #PBS -l $walltime
    ${EMAIL}
    #PBS -m a
    #PBS -o ${logdir}/\${timestamp}_alignfailclean_${groupname}.log
    #PBS -j oe
    #PBS -W depend=afternotok\${jobIDstring}
    #PBS -N Alncln${groupname}

    date +"%Y-%m-%d %H:%M:%S"
    echo "Error with alignment jobs, deleting all remaining jobs of this pipeline."
    RemJob=\$(qstat |grep ${groupname} |grep " Q \| H \| R " | awk 'BEGIN{FS=" "}{print \$1}')
    echo \${RemJob}
    qdel \${RemJob}
    
CKALIGNFAILCLN

ALIGNWRAP
    #done fastq alignment && alignment jobs failure checking.
    
fi
wait    
jID_alignwrap=$( qstat | grep AlnWrp${groupname} | cut -d ' ' -f 1 )
if [ -z $merge ]
then
    waitstring_mrgsrtwrp="#PBS -W depend=afterok:${jID_alignwrap}"
fi
echo "below is the jID_alignwrap jobid"
echo ${jID_alignwrap}
echo ${waitstring_mrgsrtwrp}
if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    ## merge the sorted files into one giant file that is also sorted.
    ## change queue below to $long_queue
    timestamp=$(date +"%s" | cut -c 4-10)
    qsub <<MRGSRTWRAP
    #PBS -S /bin/bash
    #PBS -q $queue  
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=24gb
    #PBS -l $walltime
    ${EMAIL}
    #PBS -m a
    #PBS -o ${logdir}/${timestamp}_mergesortwrap_${groupname}.log
    #PBS -j oe
    #PBS -N MStWrp${groupname}
    ${waitstring_mrgsrtwrp}
    date +"%Y-%m-%d %H:%M:%S"
    echo "all alignment done, all aplitting and alignment jobs succeeded!" 
    jID_alnOK=\$( qstat | grep AlnOK_${groupname} | cut -d ' ' -f 1 )
    echo "jID_aln-OK job id is \$jID_alnOK "
    timestamp=\$(date +"%s" | cut -c 4-10)
    if [ -z $merge ]
    then
        waitstring_alnOK="#PBS -W depend=afterok:\${jID_alnOK}" 
    fi
    echo "waitstring_anlOK is \${waitstring_alnOK}"  
    echo "below without backslash"
    echo ${waitstring_alnOK}
    echo "below with backslash"
    echo \${waitstring_alnOK}
    qsub <<MRGSRT
        #PBS -S /bin/bash
        #PBS -q $queue  
        #PBS -l nodes=1:ppn=1:AMD
        #PBS -l mem=24gb
        #PBS -l $walltime
        ${EMAIL}
        #PBS -m a
        #PBS -o ${logdir}/\${timestamp}_fragmerge_${groupname}.log
        #PBS -j oe
        #PBS -N frgmrg${groupname}
        \${waitstring_alnOK}
        date +"%Y-%m-%d %H:%M:%S"
        export LC_ALL=C
        if [ -d $donesplitdir ]
        then
            mv $donesplitdir/* $splitdir/.
        fi
        if ! sort -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt  > $outputdir/merged_sort.txt
        then 
            echo "***! Some problems occurred somewhere in creating  sorted align files."
        else
            echo "Finished sorting all sorted files into a single merge."
            rm -r ${tmpdir}        
        fi
MRGSRT
        jID_mrgsrt=\$( qstat | grep frgmrg${groupname} | cut -d ' ' -f 1 )

        ##kill all remaining jobs if previous mergesort step exited with error
        timestamp=\$(date +"%s" | cut -c 4-10)
        qsub <<MRGSRTFAILCK
        #PBS -S /bin/bash
        #PBS -q $queue  
        #PBS -l nodes=1:ppn=1:AMD
        #PBS -l mem=2gb
        #PBS -l $walltime
        ${EMAIL}
        #PBS -m a
        #PBS -o ${logdir}/\${timestamp}_clean1_${groupname}.log
        #PBS -j oe
        #PBS -N clean1${groupname}
        #PBS -W depend=afternotok:\${jID_mrgsrt}

        date +"%Y-%m-%d %H:%M:%S"
        echo "Error with merging sorted files job, ${jID_mrgsort}, deleting all remaining jobs of this pipeline."
        RemJob1=\$(qstat |grep "$groupname" |grep " Q \| H \| R " | awk 'BEGIN{FS=" "}{print $1}'| cut -d ' ' -f 1)
        qdel \${RemJob1}
MRGSRTFAILCK
MRGSRTWRAP
fi
wait
jID_mrgsrtwrap=$( qstat| grep MStWrp${groupname} | cut -d ' ' -f 1 )
if [ -z $dedup ]
then
    waitstring_RDpWrp="#PBS -W depend=afterok:${jID_mrgsrtwrap}"
fi
echo "waitstring_RDpWrp is below:"
echo ${waitstring_RDpWrp}
wait
if [ -z $final ] && [ -z $postproc ]
then
    ##remove duplicates from the big sorted file if merge sorted job exited successfully
    timestamp=$(date +"%s" | cut -c 4-10)
    qsub <<RMDUPWRAP
    #PBS -S /bin/bash
    #PBS -q $queue  
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=4gb
    #PBS -l $walltime
    ${EMAIL}
    #PBS -m a
    #PBS -o ${logdir}/${timestamp}_rmdupwrap_${groupname}.log
    #PBS -j oe
    #PBS -N RDpWrp${groupname}
    ${waitstring_RDpWrp}
    date +"%Y-%m-%d %H:%M:%S"
    echo ${waitstring_RDpWrp}
    jID_mrgsrt=\$( qstat | grep frgmrg${groupname} | cut -d ' ' -f 1 )
    if [ -z $dedup ]
    then
        waitstring_osplit="#PBS -W depend=afterok:\${jID_mrgsrt}"
    fi
    echo "jID_mrgsrt jobid is \$jID_mrgsrt "
    echo "waitstring_osplit is:\${waitstring_osplit}"
    timestamp=\$(date +"%s" | cut -c 4-10)
    qsub <<RMDUPLICATE
        #PBS -S /bin/bash
        #PBS -q $queue  
        #PBS -l nodes=1:ppn=1:AMD
        #PBS -l mem=4gb
        #PBS -l $walltime
        ${EMAIL}
        #PBS -m a   
        #PBS -o ${logdir}/\${timestamp}_osplit_${groupname}.log
        #PBS -j oe
        #PBS -N osplit${groupname}
        #PBS -v timestamp=\${timestamp}
        \${waitstring_osplit}
        date +"%Y-%m-%d %H:%M:%S"
        echo "Sucess: All mergefragments jobs were successfully finished!"
        echo "now starts to remove duplicates from the big sorted file"
        awk -v queue=${long_queue} -v outfile=${logdir}/\${timestamp}_awksplit_rmdunps -v juicedir=${juiceDir} -v dir=$outputdir -v groupname=$groupname -v walltime=$long_walltime -f ${juiceDir}/scripts/split_rmdups.awk $outputdir/merged_sort.txt    
RMDUPLICATE
        
RMDUPWRAP
fi

jID_rmdupwrap=$( qstat | grep RDpWrp${groupname} | cut -d ' ' -f 1 )
echo "jID_rmdupwrap ID: $jID_rmdupwrap"
wait
if [ -z "$genomePath" ]
then
    #If no path to genome is given, use genome ID as default.
    genomePath=$genomeID
fi

# if early exit, we stop here, once the merged_nodups.txt file is created.
if [ -z "$earlyexit" ]
then
    waitstring0="#PBS -W depend=afterok:${jID_rmdupwrap}"
   
    #Skip if post-processing only is required
    if [ -z $postproc ]
    then
        if [ -z $final ]
        then
            echo "final not set, superwrap1 job depend=afterok:jID_rmdupwrap"
        else
            waitstring0=""
        fi
        echo "waitstring0 is: $waitstring0"        
        timestamp=$(date +"%s" | cut -c 4-10)
		qsub <<SUPERWRAP1
        #PBS -S /bin/bash
        #PBS -q $queue
        #PBS -l nodes=1:ppn=1:AMD
        #PBS -l mem=1gb
        #PBS -l $walltime
        ${EMAIL}
        #PBS -m a
        #PBS -o ${logdir}/${timestamp}_superwrap1_${groupname}.log
        #PBS -j oe
        #PBS -N SpWrp1${groupname}
        ${waitstring0}
        timestamp=$(date +"%s" | cut -c 4-10)		
        echo "start submitting the lauch job!"
        export groupname=$groupname
        export juiceDir=$juiceDir
        export load_java="${load_java}"
        export about=$about
        export site_file=$site_file
        export ligation=$ligation
        export logdir=${logdir}
        export outputdir=$outputdir
        export splitdir=$splitdir
        export nofrag=$nofrag
        export genomePath=$genomePath
        export final=$final
        export queue=$queue
        export walltime=$walltime
        export long_walltime=$long_walltime
        export nofrag=$nofrag
        ${juiceDir}/scripts/launch_stats.sh
SUPERWRAP1
    fi
    jID_superwrap1="$( qstat | grep SpWrp1${groupname} | cut -d ' ' -f 1 )"
    wait
    if [ -z $postproc ]
    then
        waitstring6="#PBS -W depend=afterany:${jID_superwrap1}"
    fi
    wait
    echo "waitstring6 is : ${waitstring6}"
    
    timestamp=$(date +"%s" | cut -c 4-10)
    qsub <<SUPERWRAP2
    #PBS -S /bin/bash
    #PBS -q $queue
    #PBS -l $walltime
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=4gb
    #PBS -o ${logdir}/${timestamp}_super_wrap2_${groupname}.log
    #PBS -j oe
    #PBS -N SpWrp2${groupname}
    ${EMAIL}
    #PBS -m a
    ${waitstring6}    
    wait
    export groupname=$groupname
    export load_java="${load_java}"
    export load_cuda="${load_cuda}"
    export juiceDir=$juiceDir
    export genomeID=$genomeID
    export outputdir=$outputdir
    export logdir=${logdir}
    export splitdir=$splitdir
    export queue=$queue
    export walltime=$walltime
    export long_walltime=$long_walltime
    export postproc=$postproc
    jID_launch=\$(qstat | grep Lnch_${groupname} | cut -d ' ' -f 1)
    echo \$jID_launch
    echo "waitstring3 is \${waitstring3}"
    ${juiceDir}/scripts/postprocessing.sh
SUPERWRAP2



## After removing duplicates,if early exit is set,we directly go to the final check step.
else
    echo "earlyexit is set, stat,hic, and postprocess were not done."
    timestamp=$(date +"%s" | cut -c 4-10)
    qsub <<FINCK2
    #PBS -S /bin/bash
    #PBS -q $queue
    #PBS -l $walltime
    #PBS -l nodes=1:ppn=1:AMD 
    #PBS -l mem=1gb
    #PBS -o ${logdir}/${timestamp}_prep_done_${groupname}.out
    #PBS -j oe
    ${EMAIL}
    #PBS -m a
    #PBS -N prepd_${groupname}
    #PBS -W depend=afterok:${jID_rmdupwrap}
    date +"%Y-%m-%d %H:%M:%S"

    jID_osplit=\$( qstat | grep osplit${groupname} | cut -d ' ' -f 1 )
    jID_rmsplit=\$( qstat | grep RmSplt${groupname} | cut -d ' ' -f 1)        
    wait
    timestamp=\$(date +"%s" | cut -c 4-10)
    qsub <<PREPDONE
        #PBS -S /bin/bash
        #PBS -q $queue
        #PBS -l $walltime
        #PBS -l nodes=1:ppn=1:AMD 
        #PBS -l mem=1gb
        #PBS -o ${logdir}/\${timestamp}_done_${groupname}.log
        #PBS -j oe
        #PBS -N ${groupname}_done
        ${EMAIL}
        #PBS -m a
        #PBS -W depend=afterok:\${jID_osplit}:\${jID_rmsplit}

        date +"%Y-%m-%d %H:%M:%S"        
        export splitdir=${splitdir}
        export outputdir=${outputdir}
        ${juiceDir}/scripts/check.sh
PREPDONE
FINCK2

fi
echo "Finished adding all jobs... please wait while processing."

