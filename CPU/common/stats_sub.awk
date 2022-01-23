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
      singdups=$2;
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
  if (length(singleend)==0) {
      printf("Read type: Paired End\n");
      printf("Sequenced Read Pairs:  %'d\n", tot);
      # in the future, 1 alignment - on one end / substantially overlapping in sum
      printf("No chimera found: %'d (%0.2f%)\n",  unm, unm*100/tot);
      printf(" One or both reads unmapped: %'d (%0.2f%)\n",  unm, unm*100/tot);
      # 1 alignment would go here and be in sum of "no chimera found"
      printf("2 alignments: %'d (%0.2f%)\n", norm+chim, (norm+chim)*100/tot);
      printf(" 2 alignments (A...B): %'d (%0.2f%)\n", norm, norm*100/tot);
      printf(" 2 alignments (A1...A2B; A1B2...B1A2): %'d (%0.2f%)\n", chim, chim*100/tot);
      printf("3 or more alignments: %'d (%0.2f%)\n", coll, coll*100/tot);
     if (ligation ~ /XXXX/) {
       printf("Ligation Motif Present: N/A\n");
     }
     else {
       printf("Ligation Motif Present: %'d (%0.2f%)\n", lig, lig*100/tot);
     }
     printf("Average insert size: %0.2f\n",insertsize/NR);
     alignable=norm+chim;
     uniq=alignable-dups;
     if (alignable==0) alignable=1;
     printf("Total Unique: %'d (%0.2f%, %0.2f%)\nTotal Duplicates: %'d (%0.2f%, %0.2f%)\n", uniq, uniq*100/alignable, uniq*100/tot, dups, dups*100/alignable, dups*100/tot);
  }
  else {
      printf("Read type: Single End\n");
      printf("Sequenced Reads:  %'d\n", tot);
      # in the future, 1 alignment - on one end / substantially overlapping in sum
      printf("No chimera found: %'d (%0.2f%)\n",  unm+singleton, (unm+singleton)*100/tot);
      printf(" 0 alignments: %'d (%0.2f%)\n",  unm, unm*100/tot);
      printf(" 1 alignment: %'d (%0.2f%)\n",  singleton, singleton*100/tot);
      printf("2 alignments: %'d (%0.2f%)\n", chim, chim*100/tot);
      printf("3 or more alignments: %'d (%0.2f%)\n", coll, coll*100/tot);
     if (ligation ~ /XXXX/) {
       printf("Ligation Motif Present: N/A\n");
     }
     else {
       printf("Ligation Motif Present: %'d (%0.2f%)\n", lig, lig*100/tot);
     }
     singuniq=singleton-singdups;
     if (singleton==0) singleton=1;
     printf("1 alignment unique: %'d (%0.2f%, %0.2f%)\n", singuniq, singuniq*100/singleton, singuniq*100/tot);
     printf("1 alignment duplicates: %'d (%0.2f%, %0.2f%)\n", singdups, singdups*100/singleton, singdups*100/tot);
     #Library Complexity Estimate (1 alignment)*
     alignable=norm+chim;
     uniq=alignable-dups;
     if (alignable==0) alignable=1;
     printf("2 alignment unique: %'d (%0.2f%, %0.2f%)\n", uniq, uniq*100/alignable, uniq*100/tot);
     printf("2 alignment duplicates: %'d (%0.2f%, %0.2f%)\n", dups, dups*100/alignable, dups*100/tot);
     #Library Complexity Estimate (2 alignment)*
     printf("Total Unique: %'d (%0.2f%, %0.2f%)\n", uniq+singuniq, (uniq+singuniq)*100/(singleton+alignable), (uniq+singuniq)*100/tot);
     printf("Total Duplicates: %'d (%0.2f%, %0.2f%)\n", dups+singdups, (dups+singdups)*100/(singleton+alignable), (dups+singdups)*100/tot);
     #Library Complexity Estimate (1+2 above)*
  }
}
