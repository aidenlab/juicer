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
# Juicer postprocessing script.
# This will find contact domains with Arrowhead.
# Juicer 1.5

## Read arguments
usageHelp="Usage: ${0} [-h] -j <juicer_tools_file_path> -i <hic_file_path>"

printHelpAndExit() {
    echo "$usageHelp"
    exit $1
}

#set defaults
genomeID="hg19"
hic_file_path="$(pwd)/aligned/inter_30.hic"
juicer_tools_path="/broad/aidenlab/scripts/juicer_tools"

while getopts "h:j:i:" opt; do
    case $opt in
	h) printHelpAndExit 0;;
	j) juicer_tools_path=$OPTARG ;;
	i) hic_file_path=$OPTARG ;;
	[?]) printHelpAndExit 1;;
    esac
done

## Check that juicer tools exists 
if [ ! -e "${juicer_tools_path}" ]; then
  echo "***! Can't find juicer tools in ${juicer_tools_path}";
  exit 1;
fi

## Check that hic file exists    
if [ ! -e "${hic_file_path}" ]; then
  echo "***! Can't find inter_30.hic in ${hic_file_path}";
  exit 1;
fi

echo -e "${juicer_tools_path} is post-processing Hi-C for ${genomeID}\nData read from ${hic_file_path}.\n"
echo -e "ARROWHEAD:\n"
${juicer_tools_path} arrowhead ${hic_file_path} ${hic_file_path%.*}"_contact_domains"
if [ $? -ne 0 ]; then
    echo "***! Problem while running Arrowhead";
    exit 1
else
    echo -e "\n(-: Arrowhead Postprocessing successfully completed (-:"
fi
