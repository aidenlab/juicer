#!/usr/bin/awk -f    
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
# Splits the diploid output file into maternal, paternal, maternal_both,
# paternal_both, diff, and mismatch reads
# maternal_both and paternal_both have both read ends assigned
# diff has paternal assigned to one read end and maternal assigned to the other
# mismatch has at least one SNP that doesn't match either the reference or the
# alternate
# Juicer version 1.5
NF==13{
    if ($12 ~ /mismatch/ || $13 ~ /mismatch/) {
        print >> stem"mismatch.txt";
    }
    else if (($12 ~ /paternal/ && $12 ~ /maternal/) || ($13 ~ /paternal/ && $13 ~ /maternal/)) {
        print >> stem"mismatch.txt";
    }
    else if ($12 ~ /paternal/ && $13 ~ /paternal/) {
        print >> stem"paternal_both.txt";
    }
    else if ($12 ~ /maternal/ && $13 ~ /maternal/) {
        print >> stem"maternal_both.txt";
    }
    else if (($12 ~ /maternal/ && $13 ~ /paternal/) || ($12 ~ /paternal/ && $13 ~ /maternal/)) {
        print >> stem"diff.txt";
    }
    else {
        if ($12 ~ /maternal/ || $13 ~ /maternal/) {
            print >> stem"maternal.txt";
        }
        else if ($12 ~ /paternal/ || $13 ~ /paternal/) {
            print >> stem"paternal.txt";
        }
        else {
            print >> stem"problem.txt";
        }
    }
}
NF!=13 {
    print >> stem"problem.txt";
}
