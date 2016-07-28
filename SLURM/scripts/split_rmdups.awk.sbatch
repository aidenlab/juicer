#!/usr/bin/awk -f
BEGIN{
	tot=0;
	name=0;
}
{
	if (tot >= 1000000) {
		if (p1 != $1 || p2 != $2 || p4 != $4 || p5 != $5 || p8 != $8) {
			sname = sprintf("%s_msplit%04d_", groupname, name);
			sscriptname = sprintf("%s/.%s.slurm", outDir, sname);
			printf("#!/bin/bash -l\n#SBATCH -o %s/dup-split-%s.out\n#SBATCH -e %s/dup-split-%s.err\n#SBATCH -p %s\n#SBATCH -J %s_msplit0\n#SBATCH -t 1440\n#SBATCH -c 1\n#SBATCH --ntasks=1\ndate;awk -f %s/scripts/dups.awk -v name=%s/%s %s/split%04d;\ndate\n", outDir, name, outDir, name, queue, groupname, juicedir, dir, sname, dir, name) > sscriptname;
			sysstring = sprintf("sbatch %s", sscriptname);
			system(sysstring);
			name++;
			tot=0;
		}
	}
	outname = sprintf("%s/split%04d", dir, name);
	print > outname;
	p1=$1;p2=$2;p4=$4;p5=$5;p6=$6;p8=$8;
	tot++;
}
END {
	sname = sprintf("%s_msplit%04d_", groupname, name);
	sscriptname = sprintf("%s/.%s.slurm", outDir, sname);
	printf("#!/bin/bash -l\n#SBATCH -o %s/dup-split-%s.out\n#SBATCH -e %s/dup-split-%s.err\n#SBATCH -p %s\n#SBATCH -J %s_msplit0\n#SBATCH -t 1440\n#SBATCH -c 1\n#SBATCH --ntasks=1\ndate;awk -f %s/scripts/dups.awk -v name=%s/%s %s/split%04d;\ndate\n", outDir, name, outDir, name, queue, groupname, juicedir, dir, sname, dir, name) > sscriptname;
	sysstring = sprintf("sbatch %s", sscriptname);
	system(sysstring);

	sscriptname = sprintf("%s/.%s_msplit.slurm", outDir, groupname);
	printf("#!/bin/bash -l\n#SBATCH -o %s/dup-merge.out\n#SBATCH -e %s/dup-merge.err\n#SBATCH -p %s\n#SBATCH -J %s_msplit0\n#SBATCH -d singleton\n#SBATCH -t 1440\n#SBATCH -c 1\n#SBATCH --ntasks=1\ndate;\necho \"\"; cat %s/%s_msplit*_optdups.txt > %s/opt_dups.txt; cat %s/%s_msplit*_dups.txt > %s/dups.txt; cat %s/%s_msplit*_merged_nodups.txt > %s/merged_nodups.txt;\ndate\n", outDir, outDir, queue, groupname, dir, groupname, dir, dir, groupname, dir, dir, groupname, dir) > sscriptname;
	sysstring = sprintf("sbatch %s", sscriptname);
	system(sysstring);
       
        sscriptname = sprintf("%s/.%s_rmsplit.slurm", outDir, groupname);
        printf("#!/bin/bash -l\n#SBATCH -o %s/dup-rm.out\n#SBATCH -e %s/dup-rm.err\n#SBATCH -p %s\n#SBATCH -J %s_msplit0\n#SBATCH -d singleton\n#SBATCH -t 1440\n#SBATCH -c 1\n#SBATCH --ntasks=1\ndate;\nrm %s/*_msplit*_optdups.txt; rm %s/*_msplit*_dups.txt; rm %s/*_msplit*_merged_nodups.txt;rm %s/split*;\ndate\n", outDir, outDir, queue, groupname, dir, dir, dir, dir) > sscriptname;
        sysstring = sprintf("sbatch %s", sscriptname);
	system(sysstring);

	sscriptname = sprintf("%s/.%s_finalize.slurm", outDir, groupname);
	printf("#!/bin/bash -l\n#SBATCH -o %s/dup-guard-trigger.out\n#SBATCH -e %s/dup-guard-trigger.err\n#SBATCH -p %s\n#SBATCH -J %s_msplit0\n#SBATCH -d singleton\n#SBATCH -t 1440\n#SBATCH -c 1\n#SBATCH --ntasks=1\necho %s %s %s %s;\nsqueue -u %s;\ndate\n", outDir, outDir, queue, groupname, topDir, site, genomeID, genomePath, user, user) > sscriptname;
	sysstring = sprintf("sbatch %s", sscriptname);
        system(sysstring);
        (sysstring | getline maildupid);
        var = gensub("[^0-9]", "", "g", maildupid);
	sysstring = sprintf("scontrol update JobID=%s dependency=afterok:%s", guardjid, var);
	system(sysstring);

        sscriptname = sprintf("%s/.%s_mail.slurm", outDir, groupname);
        printf("#!/bin/bash -l\n#SBATCH -o %s/dup-mail.out\n#SBATCH -e %s/dup-mail.err\n#SBATCH -p %s\n#SBATCH -J %s_msplit0\n#SBATCH -d singleton\n#SBATCH -t 1440\n#SBATCH -c 1\n#SBATCH --ntasks=1\ndate;\necho %s %s %s %s | mail -r Juicer@rice.edu -s \"Juicer pipeline finished successfully @ Rice\" -t %s@rice.edu;\ndate\n", outDir, outDir, queue, groupname, topDir, site, genomeID, genomePath, user) > sscriptname;
        sysstring = sprintf("sbatch %s", sscriptname);
        system(sysstring);
}

