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
# Helper script to determine what dups jobs need to be relaunched.  Call using
# the same flags as Juicer to have the "tmprun" script be correct
flags=$*
ls -l aligned/split* > tmpsplits
ls -l aligned/*msplit* > tmpmsplits
num=`head -n 1 tmpsplits | awk '{size=split($0,a,"t"); print size}'`
sort -tt -k"$num"n,"$num" tmpsplits > tmpsplits2
mv tmpsplits2 tmpsplits 
awk 'NR==1{split($9,a,"_");prev=a[1]a[2]; sum=$5; split(prev, b, "msplit"); num=int(b[2]);} NR>1{split($9, a, "_"); if (prev == a[1]a[2]) { sum += $5}else {split(prev, b, "msplit"); num2=int(b[2]);  while (num+1<num2) {printf("%smsplit%04d 0\n", b[1],num+1); num++}  print prev, sum; sum = $5; prev=a[1]a[2];num=num2; }}END{print prev, sum}' tmpmsplits > tmpmsplits2

#awk 'NR==1{split($9,a,"_");prev=a[1]a[2]; sum=$5;} NR>1{split($9, a, "_"); if (prev == a[1]a[2]) { sum += $5}else {print prev, sum; sum = $5; prev=a[1]a[2]; }}END{print prev, sum}' tmpmsplits > tmpmsplits2
mv tmpmsplits2 tmpmsplits 
res1=`wc -l tmpmsplits | awk '{print $1}'`
res2=`wc -l tmpsplits | awk '{print $1}'`

if [ "$res1" -eq "$res2" ]
    then
    paste tmpsplits tmpmsplits > tmptest
    awk '$5 != $11{print "Total:",$5,"Actual:",$11,"Files:",$9,$10}' tmptest
    failures=$(awk '$5 != $11{count++}END{print count}' tmptest)
    groupname="a"`date +%s`
    if [[ -z $failures ]]
    then
        echo "echo \"cat aligned/*_msplit*_dups.txt > aligned/dups.txt; cat aligned/*_msplit*_merged_nodups.txt > aligned/merged_nodups.txt; cat aligned/*_msplit*_optdups.txt > aligned/opt_dups.txt; rm aligned/*msplit*; rm aligned/split*\" | qsub -o debug/rerundups.out -e debug/rerundups.err -cwd -r y -N ${groupname}done -l h_rt=7200" > tmprun.sh
    else
        awk -v groupname=$groupname 'BEGIN{count=0}$5 != $11{split($10,a,"m"); count++; print  "echo \"awk -f /broad/aidenlab/scripts/dups.awk -v name="a[1]"_m"a[2]"_", $9, "\" | qsub -cwd -o debug/rerundups.out -e debug/rerundups.err -l h_vmem=4g -l h_rt=7200 -r y -N",groupname""count}END{holdname=groupname""1; if(count>0){ for (i=2; i<=count; i++) {holdname=holdname","groupname""i} print "echo \"cat aligned/*_msplit*_dups.txt > aligned/dups.txt; cat aligned/*_msplit*_merged_nodups.txt > aligned/merged_nodups.txt; cat aligned/*_msplit*_optdups.txt > aligned/opt_dups.txt; rm aligned/*msplit*; rm aligned/split*;\" | qsub -hold_jid "holdname" -o debug/rerundups.out -e debug/rerundups.err -cwd -r y -l h_vmem=4g -l h_rt=7200 -N", groupname"done"}}' tmptest > tmprun.sh
    fi

    echo "echo \"/broad/aidenlab/scripts/juicer.sh $flags -S final\" | qsub  -o debug/rerundups.out -e debug/rerundups.err -cwd -r y -hold_jid ${groupname}done" >> tmprun.sh;
    chmod 755 tmprun.sh
    echo "Run ./tmprun.sh if everything looks correct";
#    rm tmpsplits tmpmsplits tmptest tmprun.sh
else
    echo "Number of msplits and splits differ";
fi


#sysstring = sprintf("qsub -o %s -q %s -N %s_catsplit -j y -hold_jid %s <<EOF
#-EOF\ncat %s/%s_msplit*_optdups.txt > %s/opt_dups.txt; cat %s/%s\
#_msplit*_dups.txt > %s/dups.txt;cat %s/%s_msplit*_merged_nodups.txt > %s/merged_nodups.txt; \nEOF\n", outfile, queue, groupname, waitstring,# \
#dir, groupname, dir, dir, groupname, dir, dir, groupname, dir, dir);
#  system(sysstring);
#  sysstring = sprintf("qsub -o %s -j y -q %s -hold_jid %s_catsplit -N %s_rmsplit <<EOF
#- EOF\n rm %s/*_msplit*_optdups.txt; rm %s/*_msplit*_dups.\
#txt; rm %s/*_msplit*_merged_nodups.txt; rm %s/split*;\nEOF",outfile, queue, groupname, groupname, dir, dir, dir, dir);
#  system(sysstring);
#EOF
#EOF
