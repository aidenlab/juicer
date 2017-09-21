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

# Helper script for finding what files successfully aligned and preparing the directory for
# rerunning juicer.sh
# Juicer version 1.5
ls -l splits > ls_splits

awk '($9 ~/_R1/ && $9 ~/fastq$/) || ($9 ~/_R1/ && $9 ~/gz$/) {split($9, a, "_R1");  print a[1]a[2], $9}($9 ~/_R2/ && $9 ~/fastq$/) || ($9 ~/_R2/ && $9 ~/gz$/){split($9,a,"_R2");print a[1]a[2], $9}' ls_splits > fastq.txt
awk '{name="splits/"$1"_norm.txt.res.txt"; if ((getline line < name) > 0){ split(line, a); if (length(a)<6 || a[2] != a[3]+a[4]+a[5]+a[6] || a[2] == 0){print "mv splits/"$2, "not_done; rm -f splits/"$1"*; rm -f splits/"$2"*;"} close(name);}else {print "mv splits/"$2, "not_done; rm -f splits/"$1"*; rm -f splits/"$2"*;"}}' fastq.txt > mv_me.sh
mkdir not_done
chmod 755 mv_me.sh
./mv_me.sh
if [ "$(ls -A not_done)" ]
    then
        if [ -d done_splits ]
        then
            mv splits/* done_splits/.
            rmdir splits
        else
            mv splits done_splits
        fi
        mv not_done splits
else
    rmdir not_done
    echo "Fastqs are aligned, run juicer.sh with -S merge flag";
fi
rm fastq.txt ls_splits mv_me.sh
