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
# Helper script for taking in multiple inter.txt and producing sums
# Juicer version 1.5
$1=="Sequenced"{
gsub(/,/,"",$4);
 total=total+$4;
}
$1=="Unmapped:"{
gsub(/,/,"",$2);
 unmapped=unmapped+$2
}
$1=="Normal"{
gsub(/,/,"",$3);
 regular=regular+$3}
$1=="Chimeric" && $2=="Paired:"{
gsub(/,/,"",$3);
 normal=normal+$3}
$2=="Ambiguous:"{
gsub(/,/,"",$3);
 abnorm=abnorm+$3}
$1=="Alignable"{
gsub(/,/,"",$4);
 alignable=alignable+$4}
$1=="Unique" {
gsub(/,/,"",$3);
 dedup=dedup+$3}
$1=="PCR"{
gsub(/,/,"",$3);
 dups=dups+$3}
$1=="Optical"{
gsub(/,/,"",$3);
 optdups=optdups+$3}
END{
 printf("%s %'d\n", "Sequenced Read Pairs:", total);
 printf(" %s %'d (%0.2f%)\n", "Normal Paired:", regular, regular*100/total);
 printf(" %s %'d (%0.2f%)\n", "Chimeric Paired:", normal, normal*100/total);
 printf(" %s %'d (%0.2f%)\n", "Chimeric Ambiguous:", abnorm, abnorm*100/total);
 printf(" %s %'d (%0.2f%)\n", "Unmapped:", unmapped, unmapped*100/total);
 printf(" %s %'d (%0.2f%)\n", "Alignable (Normal+Chimeric Paired):", alignable, alignable*100/total);
 printf("%s %'d\n", "Unique Reads:", dedup);
 printf("%s %'d\n", "PCR Duplicates:", dups);
  printf("%s %'d\n", "Optical Duplicates:", optdups);
}
