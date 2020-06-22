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
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##########
#
# Single CPU version of Juicer.
#
# Alignment script. Sets the reference genome and genome ID based on the input
# arguments (default human, none). Optional arguments are description for 
# stats file, stage to relaunch at, paths to various files if 
# needed, path to scripts directory, and the top-level directory (default 
# current directory). In lieu of setting the genome ID, you can instead set the
# reference sequence and the chrom.sizes file path, but the directory 
# containing the reference sequence must also contain the BWA index files.
#
# Aligns the fastq files, sorts, and merges.
#
# If all is successful, takes the final merged file, removes name duplicates,
# removes PCR duplicates, and creates the hic job and stats job.  Final
# product will be hic file and stats file in the aligned directory.
#                                                                       
# [topDir]/fastq  - Should contain the fastq files. This code assumes that
#                   there is an "_R" in the appropriate files, i.e. *_R*.fastq
# From the top-level directory, the following two directories are created:
#                                                                              
# [topDir]/splits  - Where to write the scratch split files (fastq files and
#                    intermediate SAM files). This can be deleted after 
#                    execution.
# [topDir]/aligned - Where to write the final output files.
#
# The following globals should be set correctly before proceeding:
#
# read1str  - portion of fastq filename that indicates this is the "read 1"
#             file; used to loop over only the read 1 and within that loop,
#             also align read 2 and merge.  If this is not set correctly,
#             script will not work. The error will often manifest itself
#             through a "*" in the name because the wildcard was not able to
#             match any files with the read1str.   
#set -e ## This is causing problems; need better error detection
# Juicer version 1.6
shopt -s extglob
export LC_ALL=C

juicer_version="1.6" 
### LOAD BWA AND SAMTOOLS


# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1" 
read2str="_R2" 

## Default options, overridden by command line arguments

# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir="/opt/juicer"
# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
site="none"
# genome ID, default to human, can also be set in options
genomeID="hg19"
# description, default empty
about=""
# do not include fragment delimited maps by default
nofrag=1
# use wobble for dedupping by default (not just exact matches)
justexact=0

## Read arguments                                                     
usageHelp="Usage: ${0##*/} [-g genomeID] [-d topDir] [-s site] [-a about] \n                 [-S stage] [-p chrom.sizes path] [-y restriction site file]\n                 [-z reference genome file] [-D Juicer scripts directory]\n                 [-b ligation] [-t threads] [-h] [-f] [-j]"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \n  \"$genomeID\"); alternatively, it can be defined using the -z command"
dirHelp="* [topDir] is the top level directory (default\n  \"$topDir\")\n     [topDir]/fastq must contain the fastq files\n     [topDir]/splits will be created to contain the temporary split files\n     [topDir]/aligned will be created for the final alignment"
siteHelp="* [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" \n  (default \"$site\")"
aboutHelp="* [about]: enter description of experiment, enclosed in single quotes"
stageHelp="* [stage]: must be one of \"merge\", \"dedup\", \"final\", \"postproc\", or \"early\".\n    -Use \"merge\" when alignment has finished but the merged_sort file has not\n     yet been created.\n    -Use \"dedup\" when the files have been merged into merged_sort but\n     merged_nodups has not yet been created.\n    -Use \"final\" when the reads have been deduped into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the stats and\n     hic files"
pathHelp="* [chrom.sizes path]: enter path for chrom.sizes file"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
scriptDirHelp="* [Juicer scripts directory]: set the Juicer directory,\n  which should have scripts/ references/ and restriction_sites/ underneath it\n  (default ${juiceDir})"
refSeqHelp="* [reference genome file]: enter path for reference sequence file, BWA index\n  files must be in same directory"
ligationHelp="* [ligation junction]: use this string when counting ligation junctions"
threadsHelp="* [threads]: number of threads when running BWA alignment"
excludeHelp="* -f: include fragment-delimited maps in hic file creation"
justHelp="* -j: just exact duplicates excluded at dedupping step"
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
    echo "$excludeHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:g:a:hs:p:y:z:S:D:fjet:b:" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	d) topDir=$OPTARG ;;
	s) site=$OPTARG ;;
	R) shortreadend=$OPTARG ;;
	r) shortread=1 ;;  #use short read aligner
	a) about=$OPTARG ;;
	p) genomePath=$OPTARG ;;  
	y) site_file=$OPTARG ;;
	z) refSeq=$OPTARG ;;
	S) stage=$OPTARG ;;
	D) juiceDir=$OPTARG ;;
	f) nofrag=0 ;; #use fragment maps
	b) ligation=$OPTARG ;;
        t) threads=$OPTARG ;;
	j) justexact=1 ;;
	e) earlyexit=1 ;;
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
	chimeric) chimeric=1 ;;
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
	hg38) refSeq="${juiceDir}/references/Homo_sapiens_assembly38.fasta";;
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
        echo "***! You must define a chrom.sizes file via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refSeq; you may use \"-p hg19\" or other standard genomes";
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
        Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
	none) ligation="XXXX";;
	*)  ligation="XXXX"
	    echo "$site not listed as recognized enzyme. Using $site_file as site file"
	    echo "Ligation junction is undefined"
    esac
fi

## If DNAse-type experiment, no fragment maps
if [ "$site" == "none" ]
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
if [ ! -z "$threads" ]
then
    threadstring="-t $threads"
else
    threads="$(getconf _NPROCESSORS_ONLN)"
    threadstring="-t $threads"
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

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

testname=$(ls -l ${fastqdir} | awk 'NR==1{print $9}')
if [ "${testname: -3}" == ".gz" ]
then
    read1=${splitdir}"/*${read1str}*.fastq.gz"
    gzipped=1
else
    read1=${splitdir}"/*${read1str}*.fastq"
fi

headfile=${outputdir}/header
date > $headfile
# Experiment description
if [ -n "${about}" ]
then
    echo -ne 'Experiment description: ${about}; ' >> $headfile
else
    echo -ne 'Experiment description: ' >> $headfile
fi

# Get version numbers of all software
echo -ne "Juicer version $juicer_version;" >> $headfile
bwa 2>&1 | awk '$1=="Version:"{printf(" BWA %s; ", $2)}' >> $headfile 
echo -ne "$threads threads; " >> $headfile
java -version 2>&1 | awk 'NR==1{printf("%s; ", $0);}' >> $headfile 
${juiceDir}/scripts/juicer_tools -V 2>&1 | awk '$1=="Juicer" && $2=="Tools"{printf("%s; ", $0);}' >> $headfile
echo "$0 $@" >> $headfile

## ALIGN FASTQ AS SINGLE END, SORT BY READNAME, HANDLE CHIMERIC READS
 
## Not in merge, dedup, final, or postproc stage, i.e. need to align files.
if [ -z $merge ] && [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    echo -e "(-: Aligning files matching $fastqdir\n to genome $genomeID with site file $site_file"
    if [ ! $splitdirexists ]
    then
        echo "(-: Created $splitdir and $outputdir."
        filename=$(basename "$i")
        filename=${filename%.*}
	ln -s ${fastqdir} ${splitdir}/.
    else
        echo -e "---  Using already created files in $splitdir\n"
    fi

    ## Loop over all read1 fastq files and create jobs for aligning read1,
    ## aligning read2, and merging the two. Keep track of merge names for final
    ## merge. When merge jobs successfully finish, can launch final merge job.
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

	source ${juiceDir}/scripts/common/countligations.sh
	if [ -z "$chimeric" ]
        then
        # Align fastq pair
            echo "bwa mem -SP5M $threadstring $refSeq $name1$ext $name2$ext > $name$ext.sam"

            bwa mem -SP5M $threadstring $refSeq $name1$ext $name2$ext > $name$ext.sam
            if [ $? -ne 0 ]
            then
                echo "***! Alignment of $name1$ext $name2$ext failed."
                exit 1
            else                                                            
		echo "(-:  Align of $name$ext.sam done successfully"
            fi                                    
        fi                                                              
    
        # call chimeric_blacklist.awk to deal with chimeric reads; 
        # sorted file is sorted by read name at this point
	touch $name${ext}_abnorm.sam $name${ext}_unmapped.sam  
	awk -v "fname1"=$name${ext}_norm.txt -v "fname2"=$name${ext}_abnorm.sam -v "fname3"=$name${ext}_unmapped.sam -f ${juiceDir}/scripts/common/chimeric_blacklist.awk $name$ext.sam
	if [ $? -ne 0 ]
	then
            echo "***! Failure during chimera handling of $name${ext}"
            exit 1
	fi
        # if any normal reads were written, find what fragment they correspond to 
        # and store that
	if [ -e "$name${ext}_norm.txt" ] && [ "$site" != "none" ] && [ -e "$site_file" ]
	then
            ${juiceDir}/scripts/common/fragment.pl $name${ext}_norm.txt $name${ext}.frag.txt $site_file                                                                
	elif [ "$site" == "none" ] || [ "$nofrag" -eq 1 ] 
	then
            awk '{printf("%s %s %s %d %s %s %s %d", $1, $2, $3, 0, $4, $5, $6, 1); for (i=7; i<=NF; i++) {printf(" %s",$i);}printf("\n");}' $name${ext}_norm.txt > $name${ext}.frag.txt
	else                                                                    
            echo "***! No $name${ext}_norm.txt file created"
            exit 1
	fi                                                                      
	if [ $? -ne 0 ]
	then
            echo "***! Failure during fragment assignment of $name${ext}"
            exit 1
	fi                              
        # sort by chromosome, fragment, strand, and position                    
	sort -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $name${ext}.frag.txt > $name${ext}.sort.txt
	if [ $? -ne 0 ]
	then
            echo "***! Failure during sort of $name${ext}"
            exit 1
	else
            rm $name${ext}_norm.txt $name${ext}.frag.txt
	fi
    done
fi

#MERGE SORTED AND ALIGNED FILES
if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    if [ -d $donesplitdir ]
    then
        mv $donesplitdir/* $splitdir/.
    fi
    if ! sort -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt  > $outputdir/merged_sort.txt
    then 
        echo "***! Some problems occurred somewhere in creating sorted align files."
        exit 1
    else
        echo "(-: Finished sorting all sorted files into a single merge."
        rm -r ${tmpdir}
    fi
fi
#REMOVE DUPLICATES
if [ -z $final ] && [ -z $postproc ]
then
    touch ${outputdir}/dups.txt
    touch ${outputdir}/optdups.txt
    touch ${outputdir}/merged_nodups.txt
    
    if [ "$justexact" -eq 1 ]
    then
	awk -f ${juiceDir}/scripts/common/dups.awk -v name=${outputdir}/ -v nowobble=1 ${outputdir}/merged_sort.txt
    else
	awk -f ${juiceDir}/scripts/common/dups.awk -v name=${outputdir}/ ${outputdir}/merged_sort.txt
    fi
    # for consistency with cluster naming in split_rmdups
    mv ${outputdir}/optdups.txt ${outputdir}/opt_dups.txt 
fi

#STATISTICS
#Skip if post-processing only is required
if [ -z $postproc ]
then        
    export _JAVA_OPTIONS=-Xmx16384m
    export LC_ALL=en_US.UTF-8 
    tail -n1 $headfile | awk '{printf"%-1000s\n", $0}' > $outputdir/inter.txt;
    cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/common/stats_sub.awk >> $outputdir/inter.txt
    ${juiceDir}/scripts/common/juicer_tools LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt
    cp $outputdir/inter.txt $outputdir/inter_30.txt 
    ${juiceDir}/scripts/common/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt 
    cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam
    cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam
    awk -f ${juiceDir}/scripts/common/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt
    # Collisions dedupping: two pass algorithm, ideally would make one pass
    gawk -v fname=$outputdir/collisions.txt -f ${juiceDir}/scripts/common/collisions_dedup_rearrange_cols.awk $outputdir/collisions.txt | sort -T $tmpdir -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n | awk -v name=$outputdir/ -f ${juiceDir}/scripts/common/collisions_dups.awk
fi

if [ -z "$genomePath" ]
then
    #If no path to genome is give, use genome ID as default.
    genomePath=$genomeID
fi

#CREATE HIC FILES
# if early exit, we stop here, once the statistics are calculated
if [ -z "$earlyexit" ]
then
    #Skip if post-processing only is required
    if [ -z $postproc ]
    then        
        if [ "$nofrag" -eq 1 ]
        then 
            ${juiceDir}/scripts/common/juicer_tools pre -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath
        else 
            ${juiceDir}/scripts/common/juicer_tools pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath 
        fi 
        ${juiceDir}/scripts/common/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt
        if [ "$nofrag" -eq 1 ]
        then 
            ${juiceDir}/scripts/common/juicer_tools pre -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath 
        else 
            ${juiceDir}/scripts/common/juicer_tools pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath
        fi
    fi
    # POSTPROCESSING
    ${juiceDir}/scripts/common/juicer_postprocessing.sh -j ${juiceDir}/scripts/common/juicer_tools -i ${outputdir}/inter_30.hic -m ${juiceDir}/references/motif -g ${genomeID}
fi
#CHECK THAT PIPELINE WAS SUCCESSFUL
export early=$earlyexit
export splitdir=$splitdir
source ${juiceDir}/scripts/common/check.sh
