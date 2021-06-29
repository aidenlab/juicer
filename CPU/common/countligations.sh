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
#
# Small helper script to count reads with ligation junction
# Juicer version 2.0
export LC_ALL=C
export LC_COLLATE=C

if [ "$usegzip" -eq 1 ]
then 
    if [ "$ligation" = "XXXX" ]
    then
	num1="0"
    else
	if [ $singleend -eq 1 ]
	then
	    num1=$(paste <(gunzip -c $name1$ext) | awk '!((NR+2)%4)' | grep -cE $ligation)
	else
	    num1=$(paste <(gunzip -c $name1$ext) <(gunzip -c $name2$ext) | awk '!((NR+2)%4)' | grep -cE $ligation)
	fi
    fi
    num2=$(gunzip -c ${name1}${ext} | wc -l | awk '{print $1}')
else
    if [ "$ligation" = "XXXX" ]
    then
	num1="0"
    else
	if [ $singleend -eq 1 ]
	then
	    num1=$(paste $name1$ext | awk '!((NR+2)%4)' | grep -cE $ligation)
	else
	    num1=$(paste $name1$ext $name2$ext | awk '!((NR+2)%4)' | grep -cE $ligation)
	fi
    fi
    num2=$(wc -l ${name1}${ext} | awk '{print $1}')
fi
echo -ne "$num1 " > ${name}${ext}_norm.txt.res.txt
echo "$num2" > ${name}${ext}_linecount.txt
