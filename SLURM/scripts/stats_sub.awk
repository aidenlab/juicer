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
# Juicer version 2.0
# Number of dups must be sent in via -v dups=# or printed in file fname
BEGIN {
  if (length(fname)>0) {
    while (getline < fname) {
      dups=$1;
    }
  }
}
{
  tot+=$2; # total reads
  unm+=$3; # unmapped
  norm+=$4; # normal paired
  chim+=$5; # chimeric paired
  coll+=$6; # collisions
  lig+=$1; #ligations
  singleton+=$7; #single alignment reads
  insertsize+=$8; #average insert size
}
END{
  if (tot==0) tot=1;
  printf("Sequenced Read Pairs:  %'d\n Normal Paired: %'d (%0.2f%)\n Chimeric Paired: %'d (%0.2f%)\n Chimeric Ambiguous: %'d (%0.2f%)\n Unmapped: %'d (%0.2f%)\n", tot, norm, norm*100/tot, chim, chim*100/tot, coll, coll*100/tot, unm, unm*100/tot);
  if (ligation ~ /XXXX/) {
    printf(" Ligation Motif Present: N/A\n");
  }
  else {
    printf(" Ligation Motif Present: %'d (%0.2f%)\n", lig, lig*100/tot);
  }
  alignable=norm+chim;
  printf(" Single Alignment: %'d (%0.2f%)\n Average insert size: %0.2f\nAlignable (Normal+Chimeric Paired): %'d (%0.2f%)\n",  singleton, singleton*100/tot, insertsize/NR, alignable, alignable*100/tot);
  uniq=alignable-dups;
  if (alignable==0) alignable=1;
  printf("Unique Reads: %'d (%0.2f%, %0.2f%)\nDuplicates: %'d (%0.2f%, %0.2f%)\n", uniq, uniq*100/alignable, uniq*100/tot, dups, dups*100/alignable, dups*100/tot);
}
