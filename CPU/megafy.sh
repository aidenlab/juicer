#!/bin/bash
{
#### Description: Wrapper script to calculate a "mega" hic file and accessibility track from a set of bams.
#### Usage: bash ./megafy.sh -c|--chrom-sizes <path_to_chrom_sizes_file> <path_to_merged_dedupped_bam_1> ... <path_to_merged_dedup_bam_N>.
#### Input: list of merged_dedup.bam files corresponding to individual Hi-C experiments.
#### Output: "mega" hic map, "mega" chromatin accessibility bw file, "mega" annotations [latter set to be expanded].
#### Dependencies: Java, Samtools, GNU Parallel, KentUtils, Juicer.
#### Written by: OD

echo "*****************************************************" >&1
echo "cmd log: "$0" "$* >&1
echo "*****************************************************" >&1

USAGE="
*****************************************************
Optimized mega script for ENCODE DCC hic-pipeline.

USAGE: ./megafy.sh -c|--chrom-sizes <path_to_chrom_sizes_file> --juicer-dir <path_to_juicer_dir> -i|--infer-stats <path_to_inter_30.hic_1,...path_to_inter_30.hic_N> [-g|--genome-id <genome_id>] [-q|--mapq <mapq>] [-r|--resolutions <resolutions_string>] [-t|--threads <thread_count>] [-T|--threads-hic <hic_thread_count>] [--from-stage <stage>] [--to-stage <stage>] <path_to_merged_dedup_bam_1> ... <path_to_merged_dedup_bam_N>

DESCRIPTION:
This is an optimized mega.sh script to produce aggregate Hi-C maps and chromatic accessibility tracks from multiple experiments.

ARGUMENTS:
path_to_merged_dedup_bam
						Path to bam file containing deduplicated alignments of Hi-C reads in bam format (output by Juicer2). Multiple bam files are expected to be passed as arguments.

OPTIONS:
-h|--help
						Shows this help.

ACTUAL WORK:

-q|--mapq   [mapq]
                        Mapping quality threshold to be used. Deafult: 30. 

-c|--chrom-sizes [path_to_chrom_sizes_file]                         
                        Path to chrom.sizes file for the reference used when processing the individual Hi-C experiments.

-g|--genome-id [genome_id]
                        Genome id, e.g. hg38, for some of the common references used by Juicer. Used to run the motif finder.

-r|--resolutions    [string]
                        Comma-separated resolutions at which to build the hic files. Default: 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100,50,20,10.

-i|--infer-stats [inter_30.hic_1,...inter_30.hic_N]
                        Comma-separated list of hic files to use to quickly calculate the mega stats file. Make sure the files used were generated at the same mapping quality as expected of the mega-map (default: -q 30).

WORKFLOW CONTROL:

-t|--threads [num]
        				Indicate how many threads to use. Default: half of available cores as calculated by parallel --number-of-cores.

-T|--threads-hic [num]
						Indicate how many threads to use when generating the Hi-C file. Default: 24.

--juicer-dir [path_to_juicer_dir]
                        Path to Juicer directory, contains scripts/, references/, and restriction_sites/

--from-stage [pipeline_stage]
						Fast-forward to a particular stage of the pipeline. The pipeline_stage argument can be \"prep\", \"hic\", \"annotations\", \"dhs\", \"cleanup\".

--to-stage [pipeline_stage]
						Exit after a particular stage of the pipeline. The argument can be \"prep\", \"hic\", \"annotations\", \"dhs\", \"cleanup\".

*****************************************************
"

# defaults:
mapq=30
resolutionsToBuildString="-r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000,500,200,100,50,20,10"
genome_id="hg38"

# multithreading
threads=`parallel --number-of-cores`
threads=$((threads/2))
# adjust for mem usage
tmp=`awk '/MemTotal/ {threads=int($2/1024/1024/2/6-1)}END{print threads+0}' /proc/meminfo 2>/dev/null`
tmp=$((tmp+0))
([ $tmp -gt 0 ] && [ $tmp -lt $threads ]) && threads=$tmp
threadsHic=$threads

## Create temporary directory for parallelization
# tmpdir="HIC_tmp"
# if [ ! -d "$tmpdir" ]; then
#     mkdir "$tmpdir"
#     chmod 777 "$tmpdir"
# fi
# export TMPDIR=${tmpdir}

#staging
first_stage="prep"
last_stage="cleanup"
declare -A stage
stage[prep]=0
stage[hic]=1
stage[annotations]=2
stage[dhs]=3
stage[cleanup]=4

############### HANDLE OPTIONS ###############

while :; do
	case $1 in
		-h|--help)
			echo "$USAGE" >&1
			exit 0
        ;;
## OUTPUT
        -c|--chrom-sizes) OPTARG=$2
            if [ -s $OPTARG ] && [[ $OPTARG == *.chrom.sizes ]]; then
                echo "... -c|--chrom-sizes flag was triggered with $OPTARG value." >&1
                chrom_sizes=$OPTARG
            else
                echo " :( Chrom.sizes file is not found at expected location, is empty or does not have the expected extension. Exiting!" >&2
                exit 1
            fi            
            shift
        ;;
        -g|--genome-id) OPTARG=$2
            genome_id=$OPTARG
            shift
        ;;
        -q|--mapq) OPTARG=$2
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
					echo "... -q|--mapq flag was triggered, will use $OPTARG as mapping quality threshold." >&1
					mapq=$OPTARG
			else
					echo " :( Wrong syntax for mapping quality threshold parameter value. Exiting!" >&2
					exit 1
			fi        	
        	shift
        ;;        
        -r|--resolutions) OPTARG=$2
            resolutionsToBuildString="-r "$OPTARG
            shift
        ;;
        -i|--infer-stats) OPTARG=$2
            echo "... -i|--infer-stats flag was triggered with $OPTARG value." >&1
            infer_stats_from_files_str=$OPTARG
            shift
        ;;
## WORKFLOW
        -t|--threads) OPTARG=$2
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
					echo "... -t|--threads flag was triggered, will try to parallelize across $OPTARG threads." >&1
					threads=$OPTARG
			else
					echo " :( Wrong syntax for thread count parameter value. Exiting!" >&2
					exit 1
			fi        	
        	shift
        ;;
        -T|--threads-hic) OPTARG=$2
        	re='^[0-9]+$'
			if [[ $OPTARG =~ $re ]]; then
					echo "... -T|--threads-hic flag was triggered, will try to parallelize across $OPTARG threads when building hic map." >&1
					threadsHic=$OPTARG
			else
					echo " :( Wrong syntax for hic thread count parameter value. Exiting!" >&2
					exit 1
			fi        	
        	shift
        ;;
        --juicer-dir) OPTARG=$2
            if [ -d $OPTARG ]; then
                echo "... --juicer-dir flag was triggered with $OPTARG." >&1
                juicer_dir=$OPTARG
            else
				exit 1
                echo " :( Juicer folder not found at expected location. Exiting!" >&2
            fi    
            shift
        ;;
		--from-stage) OPTARG=$2
			if [ "$OPTARG" == "prep" ] || [ "$OPTARG" == "hic" ] || [ "$OPTARG" == "annotations" ] || [ "$OPTARG" == "dhs" ] || [ "$OPTARG" == "cleanup" ]; then
        		echo "... --from-stage flag was triggered. Will fast-forward to $OPTARG." >&1
        		first_stage=$OPTARG
			else
				echo " :( Whong syntax for pipeline stage. Please use prep/hic/annotations/dhs/cleanup. Exiting!" >&2
				exit 1
			fi
			shift
        ;;
		--to-stage) OPTARG=$2
			if [ "$OPTARG" == "prep" ] || [ "$OPTARG" == "hic" ] || [ "$OPTARG" == "annotations" ] || [ "$OPTARG" == "dhs" ] || [ "$OPTARG" == "cleanup" ]; then
				echo "... --to-stage flag was triggered. Will exit after $OPTARG." >&1
				last_stage=$OPTARG
			else
				echo " :( Whong syntax for pipeline stage. Please use prep/hic/annotations/dhs/cleanup. Exiting!" >&2
				exit 1			
			fi
			shift
		;;
### utilitarian
        --) # End of all options
			shift
			break
		;;
		-?*)
			echo ":| WARNING: Unknown option. Ignoring: ${1}" >&2
		;;
		*) # Default case: If no more options then break out of the loop.
			break
	esac
	shift
done

## optional TODO: give error if diploid options are invoked without a vcf file

if [[ "${stage[$first_stage]}" -gt "${stage[$last_stage]}" ]]; then
	echo >&2 ":( Please make sure that the first stage requested is in fact an earlier stage of the pipeline to the one requested as last. Exiting!"
	exit 1
fi

[ -z $chrom_sizes ] && { echo >&2 ":( Chrom.sizes file is not optional. Please use the -c flag to point to the chrom.sizes file used when running Juicer. Exiting!"; exit 1; }

[ -z $genome_id ] && { echo >&2 ":| Warning: no genome_id is provided. Please provide a genome_id if using one of the common references to be able to run the motif finder. Ignoring motif finder!"; }

############### HANDLE DEPENDENCIES ###############

## Juicer & Phaser

[ -z $juicer_dir ] && { echo >&2 ":( Juicer directory is not specified. Exiting!"; exit 1; } 

##	Java Dependency
type java >/dev/null 2>&1 || { echo >&2 ":( Java is not available, please install/add to path Java. Exiting!"; exit 1; }

##	GNU Parallel Dependency
type parallel >/dev/null 2>&1 || { echo >&2 ":( GNU Parallel support is set to true (default) but GNU Parallel is not in the path. Please install GNU Parallel or set -p option to false. Exiting!"; exit 1; }
[ $(parallel --version | awk 'NR==1{print $3}') -ge 20150322 ] || { echo >&2 ":( Outdated version of GNU Parallel is installed. Please install/add to path v 20150322 or later. Exiting!"; exit 1; }

## Samtools Dependency
type samtools >/dev/null 2>&1 || { echo >&2 ":( Samtools are not available, please install/add to path. Exiting!"; exit 1; }
ver=`samtools --version | awk 'NR==1{print \$NF}'`
[[ $(echo "$ver < 1.13" |bc -l) -eq 1 ]] && { echo >&2 ":( Outdated version of samtools is installed. Please install/add to path v 1.13 or later. Exiting!"; exit 1; }

## kentUtils Dependency
type bedGraphToBigWig >/dev/null 2>&1 || { echo >&2 ":( bedGraphToBigWig is not available, please install/add to path, e.g. from kentUtils. Exiting!"; exit 1; }

############### HANDLE ARGUMENTS ###############

bam=`echo "${@:1}"`
##TODO: check file extentions

############### MAIN #################
## 0. PREP BAM FILE: {bam}->reads.sorted.bam & reads.sorted.bam.bai

if [ "$first_stage" == "prep" ]; then

	echo "...Extracting unique paired alignments from bams and sorting..." >&1

	# make header for the merged file pipe
	parallel --will-cite "samtools view -H {} > {}_header.bam" ::: $bam
	header_list=`parallel --will-cite "printf %s' ' {}_header.bam" ::: $bam`
	samtools merge --no-PG -f mega_header.bam ${header_list}
	rm ${header_list}

	samtools cat -@ $((threads * 2)) -h mega_header.bam $bam | samtools view -u -d "rt:0" -d "rt:1" -d "rt:2" -d "rt:3" -d "rt:4" -d "rt:5" -@ $((threads * 2)) -F 0x400 -q $mapq - |  samtools sort -@ $threads -m 6G -o reads.sorted.bam
	[ `echo "${PIPESTATUS[@]}" | tr -s ' ' + | bc` -eq 0 ] || { echo ":( Pipeline failed at bam sorting. See stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }
	rm mega_header.bam

	samtools index -@ $threads reads.sorted.bam	
	[ $? -eq 0 ] || { echo ":( Failed at bam indexing. See stderr for more info. Exiting!" | tee -a /dev/stderr && exit 1; }		
	# e.g. will fail with chr longer than ~500Mb. Use samtools index -c -m 14 reads.sorted.bam

	echo ":) Done extracting unique paired alignments from bam and sorting." >&1

	[ "$last_stage" == "prep" ] && { echo "Done with the requested workflow. Exiting after prepping bam!"; exit; }
	first_stage="hic"

fi

## I. HIC: 
###     (a) reads.sorted.bam & reads.sorted.bam.bai -> merged${mapq}.txt
###     (b) $infer_stats_from_files_str -> inter_${mapq}.txt & inter_${mapq}_hists.m
###     (c) merged${mapq}.txt, inter_${mapq}.txt & inter_${mapq}_hists.m -> inter_$mapq.hic

if [ "$first_stage" == "hic" ]; then

	echo "...Generating mega hic file..." >&1
    
    ([ -f reads.sorted.bam ] && [ -f reads.sorted.bam.bai ]) || { echo ":( Files from previous stages of the pipeline appear to be missing. Exiting!" | tee -a /dev/stderr; exit 1; }

### (a) Generate mnd-type input for pre: reads.sorted.bam & reads.sorted.bam.bai -> "merged"${mapq}".txt"

    export SHELL=$(type -p bash)
    doit () {
            mapq=$2
            samtools view -@ 2 -q $mapq reads.sorted.bam $1 | awk -F '\t' -v mapq=$mapq '{ip=0; mp=0; mq=-1; cb=0; chimeric=0; for(i=12;i<=NF;i++){if($i~/^ip:i:/){ip=substr($i,6)+0}else if ($i~/^mp:i:/){mp=substr($i,6)+0}else if ($i~/^MQ:i:/){mq=substr($i,6)+0}else if ($i~/cb:Z:/){cb=i}else if ($i~/SA:Z:/){chimeric=i}}}(mq<mapq){next}$7=="*"{if(!chimeric){next}else{split(substr($chimeric,6),a,","); $7=a[1]}}$7=="="{$7=$3}!cb{next}{cbchr1=gensub("_","","g",$3);cbchr2=gensub("_","","g",$7);testcb1="cb:Z:"cbchr1"_"cbchr2"_"; testcb2="cb:Z:"cbchr2"_"cbchr1"_"}($cb!~testcb1&&$cb!~testcb2){next}(!ip||!mp){next}$7==$3{if(ip>mp){next}else if (ip==mp){keep[$1]=$3" "ip}else{print 0, $3, ip, 0, 0, $7, mp, 1};next}$7<$3{next}{print 0, $3, ip, 0, 0, $7, mp, 1 > "/dev/stderr"}END{for(rd in keep){n=split(keep[rd], a, " "); print 0, a[1], a[2], 0, 0, a[1], a[2], 1}}'
    }

    export -f doit
    
    ### extra variable contains all (small) sequences that are not already in the chrom.sizes file. There is no use in them for generating the hic file, they are generated only for consistensy (and potentially stats)
    extra=`samtools view -H reads.sorted.bam | grep '^@SQ' | sed "s/.*SN:\([^\t]*\).*/\1/g" | awk 'FILENAME==ARGV[1]{drop[$1]=1;next}!($1 in drop)' ${chrom_sizes} - | xargs`

    awk -v extra="$extra" '{print $1}END{print extra}' $chrom_sizes | parallel -j $threads --will-cite --joblog temp.log -k doit {} $mapq >"merged"${mapq}".txt" 2>"merged"${mapq}".tmp.txt"
    exitval=`awk 'NR>1{if($7!=0){c=1; exit}}END{print c+0}' temp.log`
	[ $exitval -eq 0 ] || { echo ":( Pipeline failed at generating the mega mnd file. See stderr for more info. Exiting! " | tee -a /dev/stderr && exit 1; }
    rm temp.log
    sort -k2,2 -k6,6 -S 6G --parallel=${threads} "merged"${mapq}".tmp.txt" >> "merged"${mapq}".txt" && rm "merged"${mapq}".tmp.txt"

### (b) recreate stats: $infer_stats_from_files_str -> inter_${mapq}.txt & inter_${mapq}_hists.m

    touch "inter_"${mapq}".txt" "inter_"${mapq}"_hists.m"
    if [ ! -z ${infer_stats_from_files_str} ]; then
        tmpstr=""
        IFS=',' read -ra STATS <<< "$infer_stats_from_files_str"
        for i in "${STATS[@]}"; do
            [ -s $i ] && tmpstr=$tmpstr" "$i || { echo "One or more of the files listed for shortcutting stats calculation does not exist at expected location or is empty. Continuing without the stats!"; tmpstr=""; break; }
        done
        [ "$tmpstr" != "" ] && java -jar ${juicer_dir}/scripts/merge-stats.jar "inter_"${mapq} ${tmpstr:1} 
    fi

### (c) Run pre: merged${mapq}.txt, inter_${mapq}.txt & inter_${mapq}_hists.m -> inter_$mapq.hic

    ## prep
    if [[ $threadsHic -gt 1 ]]
	then
        ## Same indexing parameter 500000?
	    [[ ! -s "merged"${mapq}"_index.txt" ]] && "${juicer_dir}"/scripts/index_by_chr.awk "merged"${mapq}".txt" 500000 > "merged"${mapq}"_index.txt"
        tempdirPre="HIC"${mapq}"_tmp" && mkdir "${tempdirPre}"
	    threadHicString="--threads $threadsHic -i merged${mapq}_index.txt -t ${tempdirPre}"
    else
        threadHicString=""
	fi

    ## resources TBD for actual mega-jobs: @Muhammad should 24 be $theadsHiC or something?
    export IBM_JAVA_OPTIONS="-Xmx50000m -Xgcthreads24"
    export _JAVA_OPTIONS="-Xmx50000m -Xms50000m"

    ## actually run pre. Passing -q does not do anything, for record.
    "${juicer_dir}"/scripts/juicer_tools pre -n -s inter_${mapq}.txt -g inter_${mapq}_hists.m -q $mapq "$resolutionsToBuildString" "$threadHicString" merged${mapq}.txt inter_${mapq}.hic "$chrom_sizes"
	"${juicer_dir}"/scripts/juicer_tools addNorm --threads $threadsHic inter_${mapq}.hic
    rm -Rf "${tempdirPre}" "inter_${mapq}.hic_cat_outputs.sh"
    ## TODO: check for failure

    echo ":) Done building mega hic files." >&1

	[ "$last_stage" == "hic" ] && { echo "Done with the requested workflow. Exiting after generating the mega.hic file!"; exit; }
	first_stage="annotations"
fi

## II. ANNOTATION PACKAGE. ONLY ADDING HICCUPS AND ARROWHEAD AS EXAMPLE: inter${mapq}.hic -> annotations
if [ "$first_stage" == "annotations" ]; then
	
	echo "...Annotating loops and domains..." >&1

    [ -f inter_${mapq}.hic ] || { echo ":( Files from previous stages of the pipeline appear to be missing. Exiting!" | tee -a /dev/stderr; exit 1; }

    [ -z $genome_id ] && { echo ":( This step requires a genome ID. Please provide (e.g. \"-g hg38\") or skip this step. Exiting!" | tee -a /dev/stderr; exit 1; }

    # Create loop lists file
    "${juicer_dir}"/scripts/juicer_hiccups.sh -j "${juicer_dir}"/scripts/juicer_tools -i inter_${mapq}.hic -m "${juicer_dir}"/references/motif -g "$genome_id"

    ##TODO: check for failure

    "${juicer_dir}"/scripts/juicer_arrowhead.sh -j "${juicer_dir}"/scripts/juicer_tools -i inter_${mapq}.hic

    ##TODO: check for failure

	echo ":) Done annotating loops and domains." >&1

	[ "$last_stage" == "annotations" ] && { echo "Done with the requested workflow. Exiting after generating loop and domain annotations!"; exit; }
	first_stage="dhs"

fi

## III. BUILD ACCESSIBILITY TRACKS: reads.sorted.bam & reads.sorted.bam.bai -> inter_${mapq}.bw
if [ "$first_stage" == "dhs" ]; then

	echo "...Building accessibility tracks..." >&1

    ([ -f reads.sorted.bam ] && [ -f reads.sorted.bam.bai ]) || { echo ":( Files from previous stages of the pipeline appear to be missing. Exiting!" | tee -a /dev/stderr; exit 1; }

    ## figure out platform
    pl=`samtools view -H reads.sorted.bam | grep '^@RG' | sed "s/.*PL:\([^\t]*\).*/\1/g" | sed "s/ILM/ILLUMINA/g;s/Illumina/ILLUMINA/g;s/LS454/454/g" | uniq`
	([ "$pl" == "ILLUMINA" ] || [ "$pl" == "454" ]) || { echo ":( Platform name is not recognized or data from different platforms seems to be mixed. Can't handle this case. Exiting!" | tee -a /dev/stderr && exit 1; }
	[ "$pl" == "ILLUMINA" ] && junction_rt_string="-d rt:2 -d rt:3 -d rt:4 -d rt:5" || junction_rt_string="-d rt:0 -d rt:1"

    export SHELL=$(type -p bash)
    export junction_rt_string=${junction_rt_string}
    doit () {
            mapq=$2
            samtools view -@ 2 ${junction_rt_string} -q $mapq -h reads.sorted.bam $1 | awk 'BEGIN{OFS="\t"}{for (i=12; i<=NF; i++) {if ($i ~ /^ip/) {split($i, ip, ":"); locus[ip[3]]++; break}}}END{for (i in locus) {print $3, i-1, i, locus[i]}}' | sort -k2,2n -S6G
    }
    export -f doit

    # mapq30 accessibility track
    awk '{print $1}' $chrom_sizes | parallel -j $threads --will-cite --joblog temp.log -k doit {} $mapq > tmp.bedgraph
    exitval=`awk 'NR>1{if($7!=0){c=1; exit}}END{print c+0}' temp.log`
	[ $exitval -eq 0 ] || { echo ":( Pipeline failed at building diploid contact maps. See stderr for more info. Exiting! " | tee -a /dev/stderr && exit 1; }
	rm temp.log
    bedGraphToBigWig tmp.bedgraph $chrom_sizes inter_${mapq}.bw
    rm tmp.bedgraph

    echo ":) Done building accessibility tracks." >&1
	
	[ "$last_stage" == "dhs" ] && { echo "Done with the requested workflow. Exiting after building accessibility tracks!"; exit; }
    first_stage="cleanup"
fi

## IX. CLEANUP
	echo "...Starting cleanup..." >&1
	#rm reads.sorted.bam reads.sorted.bam.bai
    #rm inter_${mapq}.txt inter_${mapq}_hists.m
	#rm merged${mapq}.txt
	echo ":) Done with cleanup. This is the last stage of the pipeline. Exiting!"
	exit

}
