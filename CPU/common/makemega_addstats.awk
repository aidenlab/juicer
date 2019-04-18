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
 tot=tot+$4;
}
$1=="Unmapped:"{
gsub(/,/,"",$2);
 unm=unm+$2
}
$1=="Normal"{
gsub(/,/,"",$3);
 norm=norm+$3}
$1=="Chimeric" && $2=="Paired:"{
gsub(/,/,"",$3);
 chim=chim+$3}
$1=="MAPQ" && $2=="0:"{
gsub(/,/,"",$3);
 mapq0=mapq0+$3}
$1=="Collisions:"{
gsub(/,/,"",$2);
 coll=coll+$2}
$3=="Collisions:"{
gsub(/,/,"",$4);
 lowc=lowc+$4}
$1=="Ligation" {
gsub(/,/,"",$4);
 lig=lig+$4}
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
  printf("Sequenced Read Pairs:  %'d\n Normal Paired: %'d (%0.2f%)\n Chimeric Paired: %'d (%0.2f%)\n Collisions: %'d (%0.2f%)\n Low MAPQ Collisions: %'d (%0.2f%)\n Unmapped: %'d (%0.2f%)\n MAPQ 0: %'d (%0.2f%)\n Ligation Motif Present: %'d (%0.2f%)\nAlignable (Normal+Chimeric Paired): %'d (%0.2f%)\n Unique Reads: %'d\n PCR Duplicates: %'d\n Optical Duplicates: %'d\n", tot, norm, norm*100/tot, chim, chim*100/tot, coll, coll*100/tot, lowc, lowc*100/tot, unm, unm*100/tot, mapq0, mapq0*100/tot, lig, lig*100/tot, norm+chim, (norm+chim)*100/tot, dedup, dups, optdups );
}
