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
# Script to read in individual split outputs and put stats together
# Juicer version 1.6
{
    total += $2;    
    unmapped += $3; 
    normal_paired += $4;
    chimeric_paired += $5;
    chimeric_ambiguous += $6;
    ligations += $1;
}
END{
    printf("Sequenced Read Pairs:  %'d\n Normal Paired: %'d (%0.2f%)\n Chimeric Paired: %'d (%0.2f%)\n Chimeric Ambiguous: %'d (%0.2f%)\n Unmapped: %'d (%0.2f%)\n Ligation Motif Present: %'d (%0.2f%)\nAlignable (Normal+Chimeric Paired): %'d (%0.2f%)\n", total, normal_paired, normal_paired*100/total, chimeric_paired, chimeric_paired*100/total, chimeric_ambiguous, chimeric_ambiguous*100/total, unmapped, unmapped*100/total, ligations, ligations*100/total, normal_paired+chimeric_paired, (normal_paired+chimeric_paired)*100/total);
}
