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
# Script to clean up big repetitive files and zip fastqs. Run after you are 
# sure the pipeline ran successfully.  Run from top directory (HIC001 e.g.).
# Juicer version 2.0

rm aligned/merged0.txt 
rm aligned/merged30.txt 
rm -r splits 
testname=$(ls -l fastq | awk 'NR==1{print $7}')
if [ "${testname: -5}" == ".fastq" ]
then
    for i in fastq/*.fastq
    do
        gzip $i
    done
fi
