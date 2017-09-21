#!/bin/bash
jID_launch=$( qstat | grep Lnch_${groupname} | cut -c 1-7 )
if [ -z $postproc ]
then
    waitstring3="#PBS -W depend=afterok:${jID_launch}"
fi
echo "waitstring3 is: ${waitstring3}"

timestamp=$(date +"%s" | cut -c 4-10)
qsub <<- POSTPROCWRAP
#PBS -S /bin/bash
#PBS -q $queue
#PBS -l $walltime
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=4gb
#PBS -o ${logdir}/${timestamp}_postproc_wrap_${groupname}.log
#PBS -j oe
#PBS -N PPrWrp${groupname}
#PBS -M ${EMAIL}
#PBS -m a
${waitstring3}

date +"%Y-%m-%d %H:%M:%S"

jID_hic30=\$(qstat | grep hic30_${groupname} |cut -c 1-7)

if [ -z $postproc ]
then
    waitstring4="#PBS -W depend=afterok:\${jID_hic30}" 
fi
echo "waitstring4 is : ${waitstring4}"

timestamp=\$(date +"%s" | cut -c 4-10)
qsub <<POSTPROCESS
	#PBS -S /bin/bash
	#PBS -q $queue
	#PBS -l $long_walltime
	#PBS -l nodes=1:ppn=1:gpus=1:GPU 
	#PBS -l mem=60gb
	#PBS -o ${logdir}/\${timestamp}_postproc_${groupname}.log
	#PBS -j oe
	#PBS -N PProc_${groupname}
	#PBS -M ${EMAIL}
	#PBS -m a
    $waitstring4

    date +"%Y-%m-%d %H:%M:%S"
    $load_java
    $load_cuda
    module list 
    export _JAVA_OPTIONS=-Xmx16384m
    export LC_ALL=en_US.UTF-8

    ${juiceDir}/scripts/juicer_postprocessing.sh -j ${juiceDir}/scripts/juicer_tools -i ${outputdir}/inter_30.hic -m ${juiceDir}/references/motif -g $genomeID

POSTPROCESS
POSTPROCWRAP

jID_postprocwrap=$( qstat |grep PPrWrp${groupname} | cut -c 1-7 )
echo $jID_postprowrap
wait
timestamp=$(date +"%s" | cut -c 4-10)
qsub <<- FINCK
#PBS -S /bin/bash
#PBS -q $queue  
#PBS -l $walltime
#PBS -o ${logdir}/${timestamp}_prep_done_${groupname}.log
#PBS -j oe
#PBS -M ${EMAIL}
#PBS -m a
#PBS -N Pdone_${groupname}
#PBS -W depend=afterok:${jID_postprocwrap}

date +"%Y-%m-%d %H:%M:%S"    
jID_hic30=\$(qstat | grep hic30_${groupname} |cut -c 1-7)
jID_stats0=\$(qstat | grep stats0${groupname} |cut -c 1-7)
jID_stats30=\$(qstat | grep stats30${groupname} |cut -c 1-7)
jID_hic=\$(qstat | grep hic0_${groupname} |cut -c 1-7)
jID_postproc=\$(qstat | grep PProc_${groupname} |cut -c 1-7)

waitstring5="#PBS -W depend=afterok:\${jID_postproc}"
if [ -z $postproc ]
then
    waitstring5="#PBS -W depend=afterok:\${jID_hic30}:\${jID_stats0}:\${jID_stats30}:\${jID_hic}:\${jID_postproc}"
fi
timestamp=\$(date +"%s" | cut -c 4-10)
qsub <<DONE
	#PBS -S /bin/bash
	#PBS -q $queue
	#PBS -l $walltime
	#PBS -l nodes=1:ppn=1:AMD 
	#PBS -l mem=4gb
	#PBS -o ${logdir}/\${timestamp}_done_${groupname}.log
	#PBS -j oe
	#PBS -N done_${groupname}
	#PBS -M ${EMAIL}
	#PBS -m a
	\${waitstring5}

	date +"%Y-%m-%d %H:%M:%S"    
	export splitdir=${splitdir}
	export outputdir=${outputdir}
	${juiceDir}/scripts/check.sh
DONE
FINCK

