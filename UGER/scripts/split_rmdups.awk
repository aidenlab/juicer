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
# Dedup script that submits deduping jobs after splitting at known 
# non-duplicates 
# Juicer version 1.5
BEGIN{
	tot=0;
	name=0;
	waitstring="";

}
{
	if (tot >= 1000000) {
		if (p1 != $1 || p2 != $2 || p4 != $4 || p5 != $5 || p8 != $8) {
        sname=sprintf("%s_msplit%04d_", groupname, name);
        dupname=dir"/"sname"dups.txt";
        nodupname=dir"/"sname"merged_nodups.txt";
        optname=dir"/"sname"optdups.txt";

        sysstring = sprintf("touch %s %s %s", dupname,nodupname,optname);
        system(sysstring);
        if (justexact) {
            sysstring = sprintf("qsub -o %s -e %s -q %s -l h_rt=%s -N %s -l h_vmem=2g <<- EOF\nawk -f %s/scripts/dups.awk -v nowobble=1 -v name=%s/%s %s/split%04d;\nEOF\n", outfile, errfile, queue, queue_time, sname, juicedir, dir, sname, dir, name, dir, name);
        }
        else {
            sysstring = sprintf("qsub -o %s -e %s -q %s -l h_rt=%s -N %s -l h_vmem=2g <<- EOF\nawk -f %s/scripts/dups.awk -v name=%s/%s %s/split%04d;\nEOF\n", outfile, errfile, queue, queue_time, sname, juicedir, dir, sname, dir, name, dir, name);
        }
        system(sysstring);
        if (name==0) {
            waitstring=sprintf("%s", sname);
        }
        else {
            waitstring=sprintf("%s,%s", waitstring, sname);
        }
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
    sname=sprintf("%s_msplit%04d_", groupname, name);
    dupname=dir"/"sname"dups.txt";
    nodupname=dir"/"sname"merged_nodups.txt";
    optname=dir"/"sname"optdups.txt";

    sysstring = sprintf("touch %s %s %s", dupname,nodupname,optname);
    system(sysstring);
    if (justexact) {
        sysstring = sprintf("qsub -o %s -e %s -q %s -l h_rt=%s -N %s -l h_vmem=2g <<- EOF\nawk -f %s/scripts/dups.awk -v nowobble=1 -v name=%s/%s %s/split%04d;\nEOF\n", outfile, errfile, queue, queue_time, sname, juicedir, dir, sname, dir, name, dir, name);
    }
    else {
        sysstring = sprintf("qsub -o %s -e %s -q %s -l h_rt=%s -N %s -l h_vmem=2g <<- EOF\nawk -f %s/scripts/dups.awk -v name=%s/%s %s/split%04d;\nEOF\n", outfile, errfile, queue, queue_time, sname, juicedir, dir, sname, dir, name, dir, name);
    }
    system(sysstring);
    if (name==0) {
        waitstring=sprintf("%s", sname);
    }
    else {
        waitstring=sprintf("%s,%s", waitstring, sname);
    }
    sysstring = sprintf("qsub -o %s -e %s -q %s -l h_rt=%s -N %s_catsplit -hold_jid %s <<-EOF\nif [ -s %s ]; then echo \"Problem, dedupping did not execute successfully.\"; else  cat %s/%s_msplit*_optdups.txt > %s/opt_dups.txt; cat %s/%s_msplit*_dups.txt > %s/dups.txt;cat %s/%s_msplit*_merged_nodups.txt > %s/merged_nodups.txt; fi;\nEOF\n", outfile, errfile, queue, queue_time, groupname, waitstring, errfile, dir, groupname, dir, dir, groupname, dir, dir, groupname, dir, dir);
	system(sysstring);
	sysstring = sprintf("qsub -o %s -e %s -q %s -l h_rt=%s -hold_jid %s_catsplit -N %s_rmsplit <<- EOF\nif [ -s %s ]; then echo \"Problem, dedupping did not execute successfully.\"; else rm %s/*_msplit*_optdups.txt; rm %s/*_msplit*_dups.txt; rm %s/*_msplit*_merged_nodups.txt; rm %s/split*; fi \nEOF", outfile, errfile, queue, queue_time, groupname, groupname, errfile, dir, dir, dir, dir);
	system(sysstring);

}
