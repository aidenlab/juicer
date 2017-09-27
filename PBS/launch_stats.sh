#!/bin/bash
echo "check the variables:"
echo "${juiceDir}/scripts/launch_stats.sh"
echo "genomename: $groupname"
echo "juiceDir: $juiceDir"
echo "about: $about"
echo "site_file: $site_file"
echo "logdir: ${logdir}"
echo "ligation: $ligation"
echo "outputdir: $outputdir"
echo "splitdir: $splitdir"
echo "nofrag: $nofrag"
echo "final: $final"
echo "genomePath: $genomePath"
echo "queue $queue"
echo "resolutions: $resolutions"
echo "walltime: $walltime"
echo "long_walltime: $long_walltime"
echo $load_java
jID_osplit=$( qstat | grep osplit${groupname} | cut -d ' ' -f 1 )

if [ -z $final ]
then
	 waitstring2="#PBS -W depend=afterok:${jID_osplit}"
fi

echo "jID_osplit: $jID_osplit"
echo "waitstring2 is"
echo $waitstring2
timestamp=$(date +"%s" | cut -c 4-10)
qsub <<DOSTAT
#PBS -S /bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=1gb
#PBS -l $walltime
#PBS -o ${logdir}/${timestamp}_launch_${groupname}.log
#PBS -j oe
#PBS -N Lnch_${groupname}
${waitstring2}

date +"%Y-%m-%d %H:%M:%S"
echo "Alignment and merge done, launching the stats job"
#jID_osplit=\$( qstat | grep osplit${groupname} | cut -d ' ' -f 1 )
jID_rmsplit=\$( qstat | grep RmSplt${groupname} | cut -d ' ' -f 1)

if [ -z $final ]
then
    waitstring22="#PBS -W depend=afterok:\${jID_rmsplit}"
fi
echo "waitstring22 is below"
echo ${waitstring22}
echo "jID_rmsplit value is: \${jID_rmsplit}"
echo "this is the value of waitstring22: \${waitstring22}"
timestamp=\$(date +"%s" | cut -c 4-10)
echo "start sbumitting stats job"
qsub <<STATS0
	#PBS -S /bin/bash
	#PBS -q $queue
	#PBS -l nodes=1:ppn=1:AMD
	#PBS -l mem=20gb
	#PBS -l $walltime
	\${waitstring22}
	#PBS -o ${logdir}/\${timestamp}_stats0_${groupname}.log
	#PBS -j oe
	#PBS -N stats0${groupname}

	date +"%Y-%m-%d %H:%M:%S"
	$load_java
	export _JAVA_OPTIONS=-Xmx16384m; 
	export LC_ALL=en_US.UTF-8

	tail -n1 $headfile | awk '{printf"%-1000s\n", \\\$0}' > $outputdir/inter.txt 
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/stats_dups.txt $outputdir/dups.txt
	cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter.txt
	java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt

STATS0

qsub <<STATS30
	#PBS -S /bin/bash
	#PBS -q $queue
	#PBS -l nodes=1:ppn=1:AMD
	#PBS -l mem=20gb
	#PBS -l $walltime
	#PBS -o ${logdir}/\${timestamp}_stats_${groupname}.log
	#PBS -j oe
	#PBS -N stats30${groupname}
	\${waitstring22}
	
	date +"%Y-%m-%d %H:%M:%S"
	$load_java
	export _JAVA_OPTIONS=-Xmx16384m; 
	export LC_ALL=en_US.UTF-8
	echo -e 'Experiment description: $about' > $outputdir/inter_30.txt; 
	tail -n1 $headfile | awk '{printf"%-1000s\n", \\\$0}' > $outputdir/inter_30.txt 
	cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter_30.txt
	java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter_30.txt >> $outputdir/inter_30.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt

STATS30


qsub <<- ABNORMAL
	#PBS -S /bin/bash
	#PBS -q $queue
	#PBS -l nodes=1:ppn=1:AMD
	#PBS -l mem=60gb
	#PBS -l $long_walltime
	#PBS -o ${logdir}/\${timestamp}_abnormal_${groupname}.log
	#PBS -j oe
	#PBS -N abnorm_${groupname}
	\$waitstring22
	
	cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam
	cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam
	awk -f ${juiceDir}/scripts/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt
ABNORMAL

jID_stats0=\$( qstat | grep stats0${groupname} | cut -d ' ' -f 1)

wait
echo "this is the value of jID_stats: \${jID_stats}"
timestamp=\$(date +"%s" | cut -c 4-10)
echo "start submitting hic job"
qsub <<- HICWORK
	#PBS -S /bin/bash
	#PBS -q $queue
	#PBS -l nodes=1:ppn=1:AMD
	#PBS -l mem=60gb
	#PBS -l $long_walltime
	#PBS -o ${logdir}/\${timestamp}_hic0_${groupname}.log
	#PBS -j oe
	#PBS -N hic0_${groupname}
	#PBS -W depend=afterok:\${jID_stats0}
	
	date +"%Y-%m-%d %H:%M:%S"
	echo "finished stats job,now launching the hic job."    
	${load_java}
	export _JAVA_OPTIONS=-Xmx16384m
	if [ \"$nofrag\" -eq 1 ]
	then 
		${juiceDir}/scripts/juicer_tools pre -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath
	else 
		${juiceDir}/scripts/juicer_tools pre -f $site_file -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath
	fi
HICWORK
jID_stats30=\$( qstat | grep stats30${groupname} | cut -d ' ' -f 1)
timestamp=\$(date +"%s" | cut -c 4-10)
echo "start submitting hic30 job"
qsub <<- HIC30WORK
	#PBS -S /bin/bash
	#PBS -q $queue
	#PBS -l nodes=1:ppn=1:AMD
	#PBS -l mem=60gb
	#PBS -l $long_walltime
	#PBS -o ${logdir}/\${timestamp}_hic30_${groupname}.log
	#PBS -j oe
	#PBS -N hic30_${groupname}
	#PBS -W depend=afterok:\${jID_stats30}
	
	date +"%Y-%m-%d %H:%M:%S"
	$load_java
	export _JAVA_OPTIONS=-Xmx16384m
	export LC_ALL=en_US.UTF-8
	echo -e 'Experiment description: $about' > $outputdir/inter_30.txt
	cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/stats_sub.awk >> $outputdir/inter_30.txt
	java -cp ${juiceDir}/scripts/ LibraryComplexity $outputdir inter_30.txt >> $outputdir/inter_30.txt
	${juiceDir}/scripts/statistics.pl -s $site_file -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt
	if [  \"$nofrag\" -eq 1 ]
	then 
	${juiceDir}/scripts/juicer_tools pre -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath
	else 
		${juiceDir}/scripts/juicer_tools pre -f $site_file -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/inter_30.hic $genomePath
	fi
HIC30WORK

DOSTAT
