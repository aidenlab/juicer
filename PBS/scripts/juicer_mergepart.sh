juiceDir="/lustre1/mzhibo/hic/apps/juicer"
# default queue, can also be set in options via -q
queue="batch"
# default queue time, can also be set in options via -Q
walltime="walltime=24:00:00"
# default long queue, can also be set in options via -l
long_queue="batch"
# default long queue time, can also be set in options via -L
long_walltime="walltime=120:00:00"
read1str="_R1"
read2str="_R2"
roupname="C$(date "+%s"|cut -c 6-11)"
## Default options, overridden by command line arguments
groupname=$(date +"%s" |cut -c 6-11)
# top level directory, can also be set in options
topDir=$(pwd)
# restriction enzyme, can also be set in options
site="DpnII"
# genome ID, default to human, can also be set in options
genomeID="dm6_clean"
# normally both read ends are aligned with long read aligner; 
# if one end is short, this is set                 
shortreadend=0
# description, default empty
about=""
nofrag=0

# Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq/*_R*.fastq"
outputdir=${topDir}"/aligned"
tmpdir=${topDir}"/HIC_tmp"
logdir=${topDir}"/logs"
read1=${splitdir}"/*${read1str}*.fastq"
    timestamp=$(date +"%s" | cut -c 4-10)
    qsub <<ALIGNWRAP
    #PBS -S /bin/bash
    #PBS -q $queue
    #PBS -l $walltime
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=2gb
    #PBS -o ${logdir}/${timestamp}_alnwrap_${groupname}.log
    #PBS -j oe
    #PBS -N AlnWrap_${groupname}
    #PBS -M mzhibo@uga.edu
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
       qsub <<- MRGALL
        #PBS -S /bin/bash
        #PBS -q $queue
        #PBS -l $long_walltime
        #PBS -l nodes=1:ppn=1:AMD
        #PBS -l mem=24gb
        #PBS -M mzhibo@uga.edu
        #PBS -m a
        #PBS -o ${logdir}/\${timestamp}_\${jname}_merge_\${countjobs}_${groupname}.log
        #PBS -j oe
        #PBS -N Mrg_\${countjobs}${groupname}
        #PBS -v name=\${name}
        #PBS -v name1=\${name1}
        #PBS -v name2=\${name2}
        #PBS -v ext=\${ext}
        #PBS -v countjobs=\${countjobs}

        date +"%Y-%m-%d %H:%M:%S"
        export LC_ALL=C

        # sort read 1 aligned file by readname 
        # remove header, add read end indicator to read name
        awk -v name1=\${name1} -v ext=\${ext} -v name=\${name} -f ${juiceDir}/scripts/read1sort.awk \${name1}\${ext}_sort.sam > \${name1}\${ext}_sort1.sam

        #echo awk 'NF >= 11{ \$1 = \$1"/1";print}' \${name1}\${ext}_sort.sam > \${name1}\${ext}_sort1.sam'

        awk -v name1=\${name2} -v ext=\${ext} -v name=\${name} -f ${juiceDir}/scripts/read2sort.awk \${name1}\${ext}_sort.sam > \${name2}\${ext}_sort1.sam
        
        #awk '{if (NF>=11) {\$1=\$1"read2"; print $0}}' \${name1}\${ext}_sort.sam > \${name1}\${ext}_sort1.sam

        #awk 'NF >= 11{ \$1 = \$1"/2";print}' \${name2}\${ext}_sort.sam > \${name2}\${ext}_sort1.sam    
        # merge the two sorted read end files
        sort -T $tmpdir -k1,1 -m \${name1}\${ext}_sort1.sam \${name2}\${ext}_sort1.sam > \${name}\${ext}.sam
        if [ $? -ne 0 ]
        then
            echo "***! Failure during merge of read files"
            echo "Merge of \${name}\${ext}.sam failed"
            exit 1
        else
            # rm \${name1}\${ext}.sa* \${name2}\${ext}.sa* \${name1}\${ext}_sort*.sam \${name2}\${ext}_sort*.sam
            echo "\${name}\$next.sam created successfully."
        fi 
MRGALL
    wait

    jID_3=\$(qstat | grep "Mrg_\${countjobs}${groupname}" | cut -c 1-7 )
    echo "merging align1 and align2 \${coutjobs} id is \${jID_3}"
    echo "starting chimeric step after alignment"
    timestamp=\$(date +"%s" | cut -c 4-10)
    qsub <<- CHIMERIC
    #PBS -S /bin/bash
    #PBS -q $queue
    #PBS -l $walltime
    #PBS -l nodes=1:ppn=1:AMD
    #PBS -l mem=24gb
    #PBS -M mzhibo@uga.edu
    #PBS -m a
    #PBS -o ${logdir}/\${timestamp}_\${jname}_chimeric_\${countjobs}_${groupname}.log
    #PBS -j oe
    #PBS -N Chmr_\${countjobs}${groupname}
    #PBS -W depend=afterok:\${jID_3}
    #PBS -v name=\${name}
    #PBS -v ext=\${ext}
 
    date +"%Y-%m-%d %H:%M:%S"
    export LC_ALL=C    
    # call chimeric_blacklist.awk to deal with chimeric reads; sorted file is sorted by read name at this point
    awk -v "fname1"=\${name}\${ext}_norm.txt -v "fname2"=\${name}\${ext}_abnorm.sam -v "fname3"=\${name}\${ext}_unmapped.sam -f ${juiceDir}/scripts/chimeric_blacklist.awk \${name}\${ext}.sam

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
        awk '{printf("%s %s %s %d %s %s %s %d", \$1, \$2, \$3, 0, \$4, \$5, \$6, 1); for (i=7; i<=NF; i++) {printf(" %s",\$i);}printf("\n");}' \${name}\${ext}_norm.txt > \${name}\${ext}.frag.txt 
    else 
        echo "***! No \${name}\${ext}_norm.txt file created"
        echo "Creation of \${name}\${ext}_norm.txt failed. for results"
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
        rm \${name}\${ext}_norm.txt \${name}\${ext}.frag.txt
    fi    
CHIMERIC
    wait

    jID_4=\$(qstat | grep "Chmr_\${countjobs}${groupname}" | cut -c 1-7)
    echo "chimeric \$countjobs id is \$jID_4"
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
    #PBS -M mzhibo@uga.edu
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
    #PBS -M mzhibo@uga.edu
    #PBS -m a
    #PBS -o ${logdir}/\${timestamp}_alignfailclean_${groupname}.log
    #PBS -j oe
    #PBS -W depend=afternotok\${jobIDstring}
    #PBS -N Alncln_${groupname}

    date +"%Y-%m-%d %H:%M:%S"
    echo "Error with alignment jobs, deleting all remaining jobs of this pipeline."
    RemJob=\$(qstat |grep ${groupname} |grep " Q \| H \| R " | awk 'BEGIN{FS=" "}{print \$1}'| cut -c 1-7)
    qdel \$RemJob
    
CKALIGNFAILCLN

ALIGNWRAP
    #done fastq alignment && alignment jobs failure checking.

