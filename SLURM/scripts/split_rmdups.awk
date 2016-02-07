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

# Dedup script that submits deduping jobs after splitting at known non-duplicate
#
BEGIN{
	tot=0;
	name=0;
}
{
	if (tot >= 1000000) {
		if (p1 != $1 || p2 != $2 || p4 != $4 || p5 != $5 || p8 != $8) {
			sname = sprintf("%s_msplit%04d_", groupname, name);
			sscriptname = sprintf("%s/.%s.slurm", outDir, sname)
			printf("#!/bin/bash\nsrun -o %s/dup.out -e %s/dup.err -p %s -J %s_msplit0 -t 1440 -d singleton awk -f %s/scripts/dups.awk -v name=%s/%s %s/split%04d\n", outDir, outDir, queue, groupname, juicedir, dir, sname, dir, name) > sscriptname;
			sysstring = sprintf("bash %s", sscriptname);
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
	printf("#!/bin/bash\nsrun -o %s/dup-1.out -e %s/dup-1.err -p %s -t 1440 -J %s_msplit0 -d singleton -t 1440 awk -f %s/scripts/dups.awk -v name=%s/%s %s/split%04d\n", outDir, outDir, queue, groupname, juicedir, dir, sname, dir, name) > sscriptname;
	sysstring = sprintf("bash %s", sscriptname);
	system(sysstring);

	sscriptname = sprintf("%s/.%s_msplit.slurm", outDir, groupname);
	printf("#!/bin/bash\nsrun -o %s/dup-2.out -e %s/dup-2.err -p %s -J %s_msplit0 -d singleton -t 1440 echo \"\"; cat %s/%s_msplit*_optdups.txt > %s/opt_dups.txt; cat %s/%s_msplit*_dups.txt > %s/dups.txt; cat %s/%s_msplit*_merged_nodups.txt > %s/merged_nodups.txt\n", outDir, outDir, queue, groupname, dir, groupname, dir, dir, groupname, dir, dir, groupname, dir) > sscriptname;
	sysstring = sprintf("bash %s", sscriptname);
	system(sysstring);
       
	sscriptname = sprintf("%s/.%s_rmsplit.slurm", outDir, groupname);
	printf("#!/bin/bash\nsrun -o %s/dup-3.out -e %s/dup-3.err -p %s -J %s_msplit0 -d singleton -t 1440 rm %s/*_msplit*_optdups.txt; rm %s/*_msplit*_dups.txt; rm %s/*_msplit*_merged_nodups.txt; rm %s/split*\n", outDir, outDir, queue, groupname, dir, dir, dir, dir) > sscriptname;
	sysstring = sprintf("bash %s", sscriptname);
	system(sysstring);
}

