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
# Script for handing chimeric reads.
#  
# Produces SAM file with additional fields:
#   - cb: chromosome block: chr1_chr2 plus strand, frag, and external position
#         chromosome block is used for sorting and dedupping
#         hic file creation is indexed on chromosome block
#         chr1 is always <= chr2
#   - rt: read type: normal paired (2 reads), chimeric paired (3 reads),
#         chimeric paired (4 reads), collisions; primary indicated within these
#   - ip, mp: interior position of the read and its mate
# 
# Chimeric reads are treated as follows:
# "Normal" chimeras (with 3 or 4 reads, two with the ligation junction), where the
# third read (w/o ligation junction) is within 20Kbp of one of the ligation
# junction reads, are kept.  All other chimeras (where either there
# are more than 4 reads, or they don't contain the ligation junction, or the
# one non-ligation junction end is far from the others, are collisions.
#
# awk -f chimeric_sam.awk -v stem=output 
# Juicer version 2.0
function bsearch(array,len,target) {  
    low=0;
    high=length(array)-1;
    while(low<=high){
        mid=int((low+high)/2);
        if (array[mid] < target) low=mid+1;
        else if (array[mid] > target) high=mid-1;
        else return mid+1;
    }
    return low;
}
# returns absolute value
function abs(value)
{
  return (value<0?-value:value);
}
# returns minimum of two values
function min(value1,value2)
{
  return (value1<value2?value1:value2);
}
# examines read1 (s1,c1,p1) versus read2 and returns true if
# the first read comes before the second read 
# this is so duplicates can be found after sorting even when the strand and 
# chromosome are the same
function less_than(s1,c1,p1,s2,c2,p2)
{

  if (c1 < c2) return 1;
  if (c1 > c2) return 0;
  # c1 == c2
  if (s1 < s2) return 1;
  if (s1 > s2) return 0;
  # s1 == s2 && c1 == c2
  if (p1 < p2) return 1;
  if (p1 > p2) return 0;
  # all are equal, doesn't matter
  return 1;
}
# adjusts position to internal position based on cigar string
function adjust(pos,strand,cigar,notprimaln)
{ 
  # count Ms,Ds,Ns,Xs,=s for sequence length 
  seqlength=0;
  currstr=cigar;
  # can look like 15M10S20M, need to count all not just first
  if (strand==0 && notprimaln==0)
  { 
    if (cigar ~ /^[0-9]+S/) {
      where = match(cigar,/^[0-9]+S/);
      cigloc = substr(cigar,where,RLENGTH-1) + 0;
      newpos = pos + cigloc;
      currstr=substr(currstr,where+RLENGTH);
    }
    else if (cigar ~ /^[0-9]+H/) {
      where = match(cigar,/^[0-9]+H/);
      cigloc = substr(cigar,where,RLENGTH-1) + 0;
      newpos = pos + cigloc;
      currstr=substr(currstr,where+RLENGTH);
    }
    else {
      newpos = pos;
    }
    where=match(currstr, /[0-9]+[M|D|N|X|=]/);
    while (where>0) {
      seqlength+=substr(currstr,where,RLENGTH-1)+0;
      currstr=substr(currstr,where+RLENGTH);
      where=match(currstr, /[0-9]+[M|D|N|X|=]/);
    }
    newpos = newpos + seqlength - 1;
  }
  else if (strand==16 && notprimaln==0)
  { 
    where=match(currstr, /[0-9]+[M|D|N|X|=]/);
    while (where>0) {
      seqlength+=substr(currstr,where,RLENGTH-1)+0;
      currstr=substr(currstr,where+RLENGTH);
      where=match(currstr, /[0-9]+[M|D|N|X|=]/);
    }
    newpos = pos - seqlength + 1;
    if (cigar ~ /[0-9]+S$/) {
      where = match(cigar,/[0-9]+S$/);
      cigloc = substr(cigar,where,RLENGTH-1) + 0;
      newpos = newpos - cigloc;
    }
    else if (cigar ~ /[0-9]+H$/) {
      where = match(cigar,/[0-9]+H$/);
      cigloc = substr(cigar,where,RLENGTH-1) + 0;
      newpos = newpos - cigloc;
    }
  }
  else {
    newpos = pos;
  }
  if (newpos<0) {
    newpos=0;
  }
  return newpos;
}
BEGIN{
  OFS="\t";
  tottot = -1; # will count first non-group
  innerpairs = 0;
  insertsizesum = 0;
  #this will change to use single stem
  samname=stem".sam";
  fname1=stem".txt";
  # if restriction site file sent in
  if (length(site_file) > 0) {
      fraglen=0;
      while ((getline < site_file) > 0) {
	  gsub(/_/, "", $1);
	  for (i=2; i<=NF; i++) {
	      chromosomes[$1][i-2]=$i;
	  }
	  if (length(NF)>fraglen) {
	      fraglen=length(NF);
	  }
      }
      if (fraglen==0) {
	  print "!** Error while reading site file",site_file,"!" > "/dev/stderr";
	  exit 1;
      }
#      print length(chromosomes["MT"]); 24
#      print chromosomes["MT"][23]; value

  }
  else fragstr = "0_1";
  # fragment str always 0_1 with site="none"
  # seqlength 9 check in header
  chrlen=9;
}
$0 ~ /^@/{
  # print SAM header to SAM file
  print;
  split($3,chrslen,":");
  # save maximum length for 0 padding
  if (length(chrslen[2]) > chrlen) {
      chrlen=length(chrslen[2]);
  }
}
$0 !~ /^@/{
  if (tottot == -1) {
      print "@PG\tID:Juicer\tVN:2.0";
  }
  # input file is sorted by read name.  Look at read name to group 
  # appropriately
  split($1,a,"/");
  if(a[1]==prev){
    # move on to next record.  look below this block for actions that occur
    # regardless.
    count++;
  }
  else {
    # deal with read pair group
    tottot++;
    if (count==3 || count==4) {
      # chimeric read
      for (j=1; j <= count; j++) {
	split(c[j], tmp);
	split(tmp[1],readname,"/");
	# backwards compatibility
	# NB: we only care if the read ends match, not the exact number
	if (length(readname)>1) {
	    read[j] = readname[2];
	}
	else {
	    # first in pair: 64
	    read[j] = (and(tmp[2], 64) > 0);
	}
	name[j] = tmp[1];
	
	# strand; Bit 16 set means reverse strand
	str[j] = and(tmp[2],16);
	# primary alignment or not; assumes -5 bwa flag is used, i.e. primary alignment is left-most
	notprimary[j] = and(tmp[2],256);
	# chromosome
	chr[j] = tmp[3];
	# get rid of "_" in chromosome name, for cb_str
	gsub(/_/, "", chr[j]);
	# position
	pos[j] = tmp[4];
	# mapq score
	m[j] = tmp[5];
	# cigar string
	cigarstr[j] = tmp[6];
	# sequence
	seq[j] = tmp[10];
        #qual[j] = tmp[11];
        
	# get rid of soft clipping to know correct position
	if (str[j] == 0 && tmp[6] ~/^[0-9]+S/ && notprimary[j]==0) {
	  split(tmp[6], cigar, "S");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	# get rid of hard clipping to know correct position
	else if (str[j] == 0 && tmp[6] ~/^[0-9]+H/ && notprimary[j]==0) {
	  split(tmp[6], cigar, "H");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	else if (str[j] == 16) {
	  # count Ms,Ds,Ns,Xs,=s for sequence length 
	  seqlength=0; 
	  currstr=tmp[6];
	  # can look like 15M10S20M, need to count all not just first
	  where=match(currstr, /[0-9]+[M|D|N|X|=]/); 
	  while (where>0) {
	    seqlength+=substr(currstr,where,RLENGTH-1)+0;
	    currstr=substr(currstr,where+RLENGTH);
	    where=match(currstr, /[0-9]+[M|D|N|X|=]/);
	  }
	  pos[j] = pos[j] + seqlength - 1;
	  # add soft clipped bases in for proper position
	  if (tmp[6] ~ /[0-9]+S$/ && notprimary[j]==0) {
	    where = match(tmp[6],/[0-9]+S$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
          else if (tmp[6] ~ /[0-9]+H$/  && notprimary[j]==0) {
	    where = match(tmp[6],/[0-9]+H$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
          }
	}
        # blacklist - if 3rd bit set (=4) it means unmapped
        mapped[j] = and(tmp[2],4) == 0; 
      }
      if (count == 4) {
        # looking for A/B...A/B
	dist[12] = abs(chr[1]-chr[2])*10000000 + abs(pos[1]-pos[2]);
	dist[13] = abs(chr[1]-chr[3])*10000000 + abs(pos[1]-pos[3]);
	dist[14] = abs(chr[1]-chr[4])*10000000 + abs(pos[1]-pos[4]);
	dist[23] = abs(chr[2]-chr[3])*10000000 + abs(pos[2]-pos[3]);
	dist[24] = abs(chr[2]-chr[4])*10000000 + abs(pos[2]-pos[4]);
	dist[34] = abs(chr[3]-chr[4])*10000000 + abs(pos[3]-pos[4]);
	
	# A/B, A/B
	if ((dist[13] < 1000 && dist[24] < 1000 && str[1]!=str[3] && str[2]!=str[4] && read[1]!=read[3] && read[2]!=read[4] && (notprimary[1]+notprimary[3])>0 && (notprimary[2]+notprimary[4])>0)) {
	  if (notprimary[1]==0) {
	    read1 = 1;
	  }
	  else {
	    read1 = 3;
	  }
	  if (notprimary[2]==0) {
            read2 = 2;
	    interiorread2 = 4;
	  }
	  else {
	    read2 = 4;
	    interiorread2 = 2;
          }
	}
	else if ((dist[14] < 1000 && dist[23] < 1000 && str[1]!=str[4] && str[2]!=str[3] && read[1]!=read[4] && read[2]!=read[3] && (notprimary[1]+notprimary[4])>0 && (notprimary[2]+notprimary[3])>0)) {
	  if (notprimary[1]==0) {
	    read1 = 1;
	  }
	  else {
	    read1 = 4;
	  }
          if (notprimary[2]==0) {
	    read2 = 2;
	    interiorread2 = 3;
	  }
	  else {
	    read2 = 3;
	    interiorread2 = 2;
	  }
	}
	else if ((dist[12] < 1000 && dist[34] < 1000 && str[1]!=str[2] && str[3]!=str[4] && read[1]!=read[2] && read[3]!=read[4] && (notprimary[1]+notprimary[2])>0 && (notprimary[3]+notprimary[4])>0)) {
	  if (notprimary[1]==0) {
	    read1 = 1;
	  }
	  else {
	    read1 = 2;
	  }
	  if (notprimary[3]==0) {
	    read2 = 3;
	    interiorread2 = 4;
	  }
	  else {
	    read2 = 4;
	    interiorread2 = 3;
	  }
	}
	else {
	  read1 = 0;
	}
	if (read1 != 0) {
	  if (mapped[read1] && mapped[read2]) {
	    count_norm++;
	    interiorpos1=adjust(pos[read1],str[read1],cigarstr[read1],notprimary[read1]);
	    interiorpos2=adjust(pos[interiorread2],str[interiorread2],cigarstr[interiorread2],notprimary[interiorread2]);

	    if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
	      # keeping for legacy, hope to eliminate
#	      print str[read1],chr[read1],interiorpos1,str[read2],chr[read2],interiorpos2,m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1],"2$",pos[read1],pos[read2] > fname1;
	      
	      if (length(site_file) > 0) {
		  frag1 = bsearch(chromosomes[chr[read1]], length(chromosomes[chr[read1]]),pos[read1]);
		  frag2 = bsearch(chromosomes[chr[read2]], length(chromosomes[chr[read2]]),pos[read2]);
		  fragstr = sprintf("%0" fraglen "d_%0" fraglen "d", frag1,frag2);
	      }

	      # chromosome block string to sort on
	      externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[read1],pos[read2]);

	      cb_str = "cb:Z:"chr[read1]"_"chr[read2]"_"fragstr"_"str[read1]"_"str[read2]"_"externalpos;
	      
	      # assign mate mapping quality, "read type", "interior position", "mate interior position"
	      print c[read1],"MQ:i:"m[read2],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:4",cb_str ;
	      print c[read2],"MQ:i:"m[read1],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:5",cb_str ;
	      for (i in c) {
		  if (i != read1 && i != read2) print c[i],"rt:A:9",cb_str;
	      }
	    }
	    else {
	      # keeping for legacy, hope to eliminate
#	      print str[read2],chr[read2],interiorpos2,str[read1],chr[read1],interiorpos1,m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2],"2$",pos[read2],pos[read1] > fname1;

	      if (length(site_file) > 0) {
		  frag1 = bsearch(chromosomes[chr[read2]], length(chromosomes[chr[read2]]),pos[read2]);
		  frag2 = bsearch(chromosomes[chr[read1]], length(chromosomes[chr[read1]]),pos[read1]);
		  fragstr = sprintf("%0" fraglen "d_%0" fraglen "d", frag1,frag2);
	      }

	      # chromosome block string to sort on
	      externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[read2],pos[read1]);
	      cb_str = "cb:Z:"chr[read2]"_"chr[read1]"_"fragstr"_"str[read2]"_"str[read1]"_"externalpos;
	      # assign mate mapping quality, "read type", "interior position", "mate interior position"
	      print c[read2],"MQ:i:"m[read1],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:4",cb_str;
	      print c[read1],"MQ:i:"m[read2],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:5",cb_str;
	      for (i in c) {
		  if (i != read1 && i != read2) print c[i],"rt:A:9",cb_str;
	      }	      
	    }
	  }
	  else {
	    for (i in c) {
	      print c[i];
	    }	
	    count_unmapped++;
	  }
	}	
	else { 
	  # chimeric read with the 4 ends > 1KB apart
	  for (i in c) {
	      print c[i],"rt:A:8";
	  }
	  count_abnorm++;
	}
      }
      else {
	dist[12] = abs(chr[1]-chr[2])*10000000 + abs(pos[1]-pos[2]);
	dist[23] = abs(chr[2]-chr[3])*10000000 + abs(pos[2]-pos[3]);
	dist[13] = abs(chr[1]-chr[3])*10000000 + abs(pos[1]-pos[3]);
	
	if ((dist[12]<1000&&str[1]!=str[2]&&read[1]!=read[2]&&(notprimary[1]+notprimary[2])>0)||(dist[13]<1000&&str[1]!=str[3]&&read[1]!=read[3]&&(notprimary[1]+notprimary[3])>0)||(dist[23]<1000&&str[2]!=str[3]&&read[2]!=read[3]&&(notprimary[2]+notprimary[3])>0)) {
	  # The paired ends look like A/B...B.  Make sure we take A and B.
	  if (read[1] == read[2]) {
	    # take the unique one "B" for sure
	    read2 = 3;
	    # take the end of "A/B" that isn't close to "B"
	    if (notprimary[1]==0) {
	      read1 = 1; #dist[13] > dist[23] ? 1:2;
	      interiorread2 = 2; #dist[13] > dist[23] ? 2:1;
	    }
	    else {
	      read1 = 2;
	      interiorread2 = 1;
	    }
	  }
	  else if (read[1] == read[3]) {
	    read2 = 2;
	    if (notprimary[1]==0) {
	      read1 = 1; #dist[12] > dist[23] ? 1:3;
	      interiorread2 = 3; #dist[12] > dist[23] ? 3:1;
	    }
	    else {
	      read1 = 3;
	      interiorread2 = 1;
	    }
	  }
	  else if (read[2] == read[3]) {
	    read2 = 1;
	    if (notprimary[2]==0) {
	      read1 = 2; #dist[12] > dist[13] ? 2:3;
	      interiorread2 = 3; #dist[12] > dist[13] ? 3:2;
	    }
	    else {
	      read1 = 3;
	      interiorread2 = 2;
	    }
	  }
	  else {
	    printf("reads strange\n") > "/dev/stderr"
	    exit 1
	  }
	  
	  if (mapped[read1] && mapped[read2]) {
	    count_norm++;
	    interiorpos1=adjust(pos[read1],str[read1],cigarstr[read1],notprimary[read1]);
	    interiorpos2=adjust(pos[interiorread2],str[interiorread2],cigarstr[interiorread2],notprimary[interiorread2]);
	    if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
#	      print str[read1],chr[read1],interiorpos1,str[read2],chr[read2],interiorpos2,m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1],"1$",pos[read1],pos[read2] > fname1;

	      if (length(site_file) > 0) {
		  frag1 = bsearch(chromosomes[chr[read1]], length(chromosomes[chr[read1]]),pos[read1]);
		  frag2 = bsearch(chromosomes[chr[read2]], length(chromosomes[chr[read2]]),pos[read2]);
		  fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag1,frag2);
	      }
	      # chromosome block string to sort on
	      externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[read1],pos[read2]);
	      cb_str = "cb:Z:"chr[read1]"_"chr[read2]"_"fragstr"_"str[read1]"_"str[read2]"_"externalpos;
	      
	      # assign mate mapping quality, "read type", "interior position", "mate interior position"
	      print c[read1],"MQ:i:"m[read2],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:2",cb_str;
	      print c[read2],"MQ:i:"m[read1],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:3",cb_str;
	      for (i in c) {
		  if (i != read1 && i != read2) print c[i],"rt:A:9",cb_str;
	      }
	    }
	    else {
#	      print str[read2],chr[read2],adjust(pos[interiorread2],str[interiorread2],cigarstr[interiorread2],notprimary[interiorread2]),str[read1],chr[read1],adjust(pos[read1],str[read1],cigarstr[read1],notprimary[read1]),m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2],"1$",pos[read2],pos[read1] > fname1;

	      if (length(site_file) > 0) {
		  frag1 = bsearch(chromosomes[chr[read2]], length(chromosomes[chr[read2]]),pos[read2]);
		  frag2 = bsearch(chromosomes[chr[read1]], length(chromosomes[chr[read1]]),pos[read1]);
		  fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag1,frag2);
	      }

	      # chromosome block string to sort on
	      externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[read2],pos[read1]);
	      cb_str = "cb:Z:"chr[read2]"_"chr[read1]"_"fragstr"_"str[read2]"_"str[read1]"_"externalpos;

	      # assign mate mapping quality, "read type", "interior position", "mate interior position"
	      print c[read2],"MQ:i:"m[read1],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:2",cb_str;
	      print c[read1],"MQ:i:"m[read2],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:3",cb_str;
	      for (i in c) {
		  if (i != read1 && i != read2) print c[i],"rt:A:9",cb_str;
	      }	      
	    }
	  }
	  else {
	    for (i in c) {
	      print c[i];
	    }	
	    count_unmapped++;
	  }
	}
	else {
	  # chimeric read with the 3 ends > 1KB apart
	  for (i in c) {
	      print c[i],"rt:A:8";
	  }
	  count_abnorm++;
	}
      }
    }
    else if (count > 3) {
      # chimeric read > 3, too many to deal with
      for (i in c) {
	print c[i],"rt:A:8";
      }
      count_abnorm++;
    }
    else if (count == 2) {
      # code here should be same as above, but it's a "normal" read
      j = 0;
      for (i in c) {
	split(c[i], tmp);
	split(tmp[1],readname,"/");
	str[j] = and(tmp[2],16);
	notprimary[j] = and(tmp[2],256);
	chr[j] = tmp[3];
	# get rid of "_" in chromosome name, for cb_str
	gsub(/_/, "", chr[j]);
	pos[j] = tmp[4];
	m[j] = tmp[5];
	cigarstr[j] = tmp[6];
	seq[j] = tmp[10];
	name[j] = tmp[1];
	
        # blacklist - if 3rd bit set (=4) it means unmapped
        mapped[j] = and(tmp[2],4) == 0; 
	
	if (str[j] == 0 && tmp[6] ~/^[0-9]+S/ && notprimary[j] == 0) {
	  split(tmp[6], cigar, "S");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	else if (str[j] == 0 && tmp[6] ~/^[0-9]+H/ && notprimary[j] == 0) {
	  split(tmp[6], cigar, "H");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	else if (str[j] == 16) {
	  # count Ms,Ds,Ns,Xs,=s for sequence length 
	  seqlength=0; 
	  currstr=tmp[6];
	  where=match(currstr, /[0-9]+[M|D|N|X|=]/); 
	  while (where>0) {
	    seqlength+=substr(currstr,where,RLENGTH-1)+0;
	    currstr=substr(currstr,where+RLENGTH);
	    where=match(currstr, /[0-9]+[M|D|N|X|=]/);
	  }
	  pos[j] = pos[j] + seqlength - 1;
	  if (tmp[6] ~ /[0-9]+S$/ && notprimary[j] == 0) {
	    where = match(tmp[6],/[0-9]+S$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
	  if (tmp[6] ~ /[0-9]+H$/ && notprimary[j] == 0) {
	    where = match(tmp[6],/[0-9]+H$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
	}
	j++;
      }
      if (mapped[0] && mapped[1]) {
	count_reg++;
	interiorpos1=adjust(pos[0],str[0],cigarstr[0],notprimary[0]);
	interiorpos2=adjust(pos[1],str[1],cigarstr[1],notprimary[1]);

	if (length(site_file) > 0) {
	    frag1 = bsearch(chromosomes[chr[0]], length(chromosomes[chr[0]]),pos[0]);
	    frag2 = bsearch(chromosomes[chr[1]], length(chromosomes[chr[1]]),pos[1]);
	    if (printme) {print "HERE"; print frag1,frag2; print length(chromosomes[chr[0]]), length(chromosomes[chr[1]])}
	    # this happens with circular chromosomes (MT)
	    # only adjusting here since it matters for internal position
	    if (frag1 >= length(chromosomes[chr[0]])) 
		frag1 = length(chromosomes[chr[0]])-1;
	    if (frag2 >= length(chromosomes[chr[1]])) 
		frag2 = length(chromosomes[chr[1]])-1;

	    # adjust internal position based on cutting site
	    if (str[0] == 0) {
		interiorpos1 = chromosomes[chr[0]][frag1];
	    }
	    else {
		if (frag1 == 0) {
		    interiorpos1 = 0;
		}
		else {
		    interiorpos1 = chromosomes[chr[0]][frag1-1];
		}
	    }
	    if (str[1] == 0) {
		interiorpos2 = chromosomes[chr[1]][frag2];
	    }
	    else {
		if (frag2 == 0) {
		    interiorpos2 = 0;
		}
		else {
		    interiorpos2 = chromosomes[chr[1]][frag2-1];
		}
	    }
	}
	if (less_than(str[0],chr[0],pos[0],str[1],chr[1],pos[1])) {
#	  print str[0],chr[0],interiorpos1,str[1],chr[1],interiorpos2,m[0],cigarstr[0],seq[0],m[1],cigarstr[1],seq[1],name[0],"0$",pos[0],pos[1] > fname1;
	    
	    if (length(site_file) > 0) {
		fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag1,frag2);
	    }

	    # chromosome block string to sort on
	    externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[0],pos[1]);
	    cb_str = "cb:Z:"chr[0]"_"chr[1]"_"fragstr"_"str[0]"_"str[1]"_"externalpos;
	  
	    # assign mate mapping quality, "read type", "interior position", "mate interior position"
	    print c[1],"MQ:i:"m[1],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:0",cb_str;
	    print c[2],"MQ:i:"m[0],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:1",cb_str;
	    
	    if (str[0]==0&&str[1]==16&&chr[0]==chr[1]&&pos[0]<pos[1]&&pos[1]-pos[0]<1000) {
		innerpairs+=1;
		insertsizesum+=abs(pos[1]-pos[0])
	    }
	}
	else {
#	  print str[1],chr[1],interiorpos2,str[0],chr[0],interiorpos1,m[1],cigarstr[1],seq[1],m[0],cigarstr[0],seq[0],name[1],"0$",pos[1],pos[0] > fname1;

	  if (length(site_file) > 0) {
	      fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag2,frag1);
	  }

	  # chromosome block string to sort on
	  externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[1],pos[0]);
	  cb_str = "cb:Z:"chr[1]"_"chr[0]"_"fragstr"_"str[1]"_"str[0]"_"externalpos;
	      
	  # assign mate mapping quality, "read type", "interior position", "mate interior position"
	  print c[2],"MQ:i:"m[0],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:0",cb_str;
	  print c[1],"MQ:i:"m[1],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:1",cb_str;

	  if (str[1]==0&&str[0]==16&&chr[0]==chr[1]&&pos[1]<pos[0]&&pos[0]-pos[1]<1000) {
		innerpairs+=1;
		insertsizesum+=abs(pos[0]-pos[1]);
	  }
	}
      }
      else {
	for (i in c) {
	  print c[i];
	}	
	count_unmapped++;
      }
    }
    else if (count == 1) {
      # this actually shouldn't happen, but it happens with alternate aligners on occasion
      count_abnorm++;
      for (i in c) {
	print c[i],"rt:A:8";
      }
    }
    # reset variables
    delete c;
    count=1;
    prev=a[1];
  }
  # these happen no matter what, after the above processing
  c[count] = $0;
}
END{
    # deal with read pair group
    tottot++;
    if (count==3 || count==4) {
      # chimeric read
      for (j=1; j <= count; j++) {
	split(c[j], tmp);
	split(tmp[1],readname,"/");
	# backwards compatibility
	# NB: we only care if the read ends match, not the exact number
	if (length(readname)>1) {
	    read[j] = readname[2];
	}
	else {
	    # first in pair: 64
	    read[j] = (and(tmp[2], 64) > 0);
	}
	name[j] = tmp[1];
	
	# strand; Bit 16 set means reverse strand
	str[j] = and(tmp[2],16);
	# primary alignment or not; assumes -5 bwa flag is used, i.e. primary alignment is left-most
	notprimary[j] = and(tmp[2],256);
	# chromosome
	chr[j] = tmp[3];
	# get rid of "_" in chromosome name, for cb_str
	gsub(/_/, "", chr[j]);
	# position
	pos[j] = tmp[4];
	# mapq score
	m[j] = tmp[5];
	# cigar string
	cigarstr[j] = tmp[6];
	# sequence
	seq[j] = tmp[10];
        #qual[j] = tmp[11];
        
	# get rid of soft clipping to know correct position
	if (str[j] == 0 && tmp[6] ~/^[0-9]+S/ && notprimary[j]==0) {
	  split(tmp[6], cigar, "S");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	# get rid of hard clipping to know correct position
	else if (str[j] == 0 && tmp[6] ~/^[0-9]+H/ && notprimary[j]==0) {
	  split(tmp[6], cigar, "H");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	else if (str[j] == 16) {
	  # count Ms,Ds,Ns,Xs,=s for sequence length 
	  seqlength=0; 
	  currstr=tmp[6];
	  # can look like 15M10S20M, need to count all not just first
	  where=match(currstr, /[0-9]+[M|D|N|X|=]/); 
	  while (where>0) {
	    seqlength+=substr(currstr,where,RLENGTH-1)+0;
	    currstr=substr(currstr,where+RLENGTH);
	    where=match(currstr, /[0-9]+[M|D|N|X|=]/);
	  }
	  pos[j] = pos[j] + seqlength - 1;
	  # add soft clipped bases in for proper position
	  if (tmp[6] ~ /[0-9]+S$/ && notprimary[j]==0) {
	    where = match(tmp[6],/[0-9]+S$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
          else if (tmp[6] ~ /[0-9]+H$/  && notprimary[j]==0) {
	    where = match(tmp[6],/[0-9]+H$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
          }
	}
        # blacklist - if 3rd bit set (=4) it means unmapped
        mapped[j] = and(tmp[2],4) == 0; 
      }
      if (count == 4) {
        # looking for A/B...A/B
	dist[12] = abs(chr[1]-chr[2])*10000000 + abs(pos[1]-pos[2]);
	dist[13] = abs(chr[1]-chr[3])*10000000 + abs(pos[1]-pos[3]);
	dist[14] = abs(chr[1]-chr[4])*10000000 + abs(pos[1]-pos[4]);
	dist[23] = abs(chr[2]-chr[3])*10000000 + abs(pos[2]-pos[3]);
	dist[24] = abs(chr[2]-chr[4])*10000000 + abs(pos[2]-pos[4]);
	dist[34] = abs(chr[3]-chr[4])*10000000 + abs(pos[3]-pos[4]);
	
	# A/B, A/B
	if ((dist[13] < 1000 && dist[24] < 1000 && str[1]!=str[3] && str[2]!=str[4] && read[1]!=read[3] && read[2]!=read[4] && (notprimary[1]+notprimary[3])>0 && (notprimary[2]+notprimary[4])>0)) {
	  if (notprimary[1]==0) {
	    read1 = 1;
	  }
	  else {
	    read1 = 3;
	  }
	  if (notprimary[2]==0) {
            read2 = 2;
	    interiorread2 = 4;
	  }
	  else {
	    read2 = 4;
	    interiorread2 = 2;
          }
	}
	else if ((dist[14] < 1000 && dist[23] < 1000 && str[1]!=str[4] && str[2]!=str[3] && read[1]!=read[4] && read[2]!=read[3] && (notprimary[1]+notprimary[4])>0 && (notprimary[2]+notprimary[3])>0)) {
	  if (notprimary[1]==0) {
	    read1 = 1;
	  }
	  else {
	    read1 = 4;
	  }
          if (notprimary[2]==0) {
	    read2 = 2;
	    interiorread2 = 3;
	  }
	  else {
	    read2 = 3;
	    interiorread2 = 2;
	  }
	}
	else if ((dist[12] < 1000 && dist[34] < 1000 && str[1]!=str[2] && str[3]!=str[4] && read[1]!=read[2] && read[3]!=read[4] && (notprimary[1]+notprimary[2])>0 && (notprimary[3]+notprimary[4])>0)) {
	  if (notprimary[1]==0) {
	    read1 = 1;
	  }
	  else {
	    read1 = 2;
	  }
	  if (notprimary[3]==0) {
	    read2 = 3;
	    interiorread2 = 4;
	  }
	  else {
	    read2 = 4;
	    interiorread2 = 3;
	  }
	}
	else {
	  read1 = 0;
	}
	if (read1 != 0) {
	  if (mapped[read1] && mapped[read2]) {
	    count_norm++;
	    interiorpos1=adjust(pos[read1],str[read1],cigarstr[read1],notprimary[read1]);
	    interiorpos2=adjust(pos[interiorread2],str[interiorread2],cigarstr[interiorread2],notprimary[interiorread2]);

	    if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
	      # keeping for legacy, hope to eliminate
#	      print str[read1],chr[read1],interiorpos1,str[read2],chr[read2],interiorpos2,m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1],"2$",pos[read1],pos[read2] > fname1;
	      
	      if (length(site_file) > 0) {
		  frag1 = bsearch(chromosomes[chr[read1]], length(chromosomes[chr[read1]]),pos[read1]);
		  frag2 = bsearch(chromosomes[chr[read2]], length(chromosomes[chr[read2]]),pos[read2]);
		  fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag1,frag2);
	      }

	      # chromosome block string to sort on
	      externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[read1],pos[read2]);
	      cb_str = "cb:Z:"chr[read1]"_"chr[read2]"_"fragstr"_"str[read1]"_"str[read2]"_"externalpos;
	      
	      # assign mate mapping quality, "read type", "interior position", "mate interior position"
	      print c[read1],"MQ:i:"m[read2],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:4",cb_str;
	      print c[read2],"MQ:i:"m[read1],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:5",cb_str;
	      for (i in c) {
		  if (i != read1 && i != read2) print c[i],"rt:A:9",cb_str;
	      }
	    }
	    else {
	      # keeping for legacy, hope to eliminate
#	      print str[read2],chr[read2],interiorpos2,str[read1],chr[read1],interiorpos1,m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2],"2$",pos[read2],pos[read1] > fname1;

	      if (length(site_file) > 0) {
		  frag1 = bsearch(chromosomes[chr[read2]], length(chromosomes[chr[read2]]),pos[read2]);
		  frag2 = bsearch(chromosomes[chr[read1]], length(chromosomes[chr[read1]]),pos[read1]);
		  fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag1,frag2);
	      }

	      # chromosome block string to sort on
	      externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[read2],pos[read1]);
	      cb_str = "cb:Z:"chr[read2]"_"chr[read1]"_"fragstr"_"str[read2]"_"str[read1]"_"externalpos;
	      # assign mate mapping quality, "read type", "interior position", "mate interior position"
	      print c[read2],"MQ:i:"m[read1],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:4",cb_str;
	      print c[read1],"MQ:i:"m[read2],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:5",cb_str;
	      for (i in c) {
		  if (i != read1 && i != read2) print c[i],"rt:A:9",cb_str;
	      }	      
	    }
	  }
	  else {
	    for (i in c) {
	      print c[i];
	    }	
	    count_unmapped++;
	  }
	}	
	else { 
	  # chimeric read with the 4 ends > 1KB apart
	  for (i in c) {
	      print c[i],"rt:A:8";
	  }
	  count_abnorm++;
	}
      }
      else {
	dist[12] = abs(chr[1]-chr[2])*10000000 + abs(pos[1]-pos[2]);
	dist[23] = abs(chr[2]-chr[3])*10000000 + abs(pos[2]-pos[3]);
	dist[13] = abs(chr[1]-chr[3])*10000000 + abs(pos[1]-pos[3]);
	
	if ((dist[12]<1000&&str[1]!=str[2]&&read[1]!=read[2]&&(notprimary[1]+notprimary[2])>0)||(dist[13]<1000&&str[1]!=str[3]&&read[1]!=read[3]&&(notprimary[1]+notprimary[3])>0)||(dist[23]<1000&&str[2]!=str[3]&&read[2]!=read[3]&&(notprimary[2]+notprimary[3])>0)) {
	  # The paired ends look like A/B...B.  Make sure we take A and B.
	  if (read[1] == read[2]) {
	    # take the unique one "B" for sure
	    read2 = 3;
	    # take the end of "A/B" that isn't close to "B"
	    if (notprimary[1]==0) {
	      read1 = 1; #dist[13] > dist[23] ? 1:2;
	      interiorread2 = 2; #dist[13] > dist[23] ? 2:1;
	    }
	    else {
	      read1 = 2;
	      interiorread2 = 1;
	    }
	  }
	  else if (read[1] == read[3]) {
	    read2 = 2;
	    if (notprimary[1]==0) {
	      read1 = 1; #dist[12] > dist[23] ? 1:3;
	      interiorread2 = 3; #dist[12] > dist[23] ? 3:1;
	    }
	    else {
	      read1 = 3;
	      interiorread2 = 1;
	    }
	  }
	  else if (read[2] == read[3]) {
	    read2 = 1;
	    if (notprimary[2]==0) {
	      read1 = 2; #dist[12] > dist[13] ? 2:3;
	      interiorread2 = 3; #dist[12] > dist[13] ? 3:2;
	    }
	    else {
	      read1 = 3;
	      interiorread2 = 2;
	    }
	  }
	  else {
	    printf("reads strange\n") > "/dev/stderr"
	    exit 1
	  }
	  
	  if (mapped[read1] && mapped[read2]) {
	    count_norm++;
	    interiorpos1=adjust(pos[read1],str[read1],cigarstr[read1],notprimary[read1]);
	    interiorpos2=adjust(pos[interiorread2],str[interiorread2],cigarstr[interiorread2],notprimary[interiorread2]);
	    if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
#	      print str[read1],chr[read1],interiorpos1,str[read2],chr[read2],interiorpos2,m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1],"1$",pos[read1],pos[read2] > fname1;

	      if (length(site_file) > 0) {
		  frag1 = bsearch(chromosomes[chr[read1]], length(chromosomes[chr[read1]]),pos[read1]);
		  frag2 = bsearch(chromosomes[chr[read2]], length(chromosomes[chr[read2]]),pos[read2]);
		  fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag1,frag2);
	      }
	      # chromosome block string to sort on
	      externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[read1],pos[read2]);
	      cb_str = "cb:Z:"chr[read1]"_"chr[read2]"_"fragstr"_"str[read1]"_"str[read2]"_"externalpos;
	      
	      # assign mate mapping quality, "read type", "interior position", "mate interior position"
	      print c[read1],"MQ:i:"m[read2],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:2",cb_str;
	      print c[read2],"MQ:i:"m[read1],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:3",cb_str;
	      for (i in c) {
		  if (i != read1 && i != read2) print c[i],"rt:A:9",cb_str;
	      }
	    }
	    else {
#	      print str[read2],chr[read2],adjust(pos[interiorread2],str[interiorread2],cigarstr[interiorread2],notprimary[interiorread2]),str[read1],chr[read1],adjust(pos[read1],str[read1],cigarstr[read1],notprimary[read1]),m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2],"1$",pos[read2],pos[read1] > fname1;

	      if (length(site_file) > 0) {
		  frag1 = bsearch(chromosomes[chr[read2]], length(chromosomes[chr[read2]]),pos[read2]);
		  frag2 = bsearch(chromosomes[chr[read1]], length(chromosomes[chr[read1]]),pos[read1]);
		  fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag1,frag2);
	      }

	      # chromosome block string to sort on
	      externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[read2],pos[read1]);
	      cb_str = "cb:Z:"chr[read2]"_"chr[read1]"_"fragstr"_"str[read2]"_"str[read1]"_"externalpos;

	      # assign mate mapping quality, "read type", "interior position", "mate interior position"
	      print c[read2],"MQ:i:"m[read1],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:2",cb_str;
	      print c[read1],"MQ:i:"m[read2],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:3",cb_str;
	      for (i in c) {
		  if (i != read1 && i != read2) print c[i],"rt:A:9",cb_str;
	      }	      
	    }
	  }
	  else {
	    for (i in c) {
	      print c[i];
	    }	
	    count_unmapped++;
	  }
	}
	else {
	  # chimeric read with the 3 ends > 1KB apart
	  for (i in c) {
	      print c[i],"rt:A:8";
	  }
	  count_abnorm++;
	}
      }
    }
    else if (count > 3) {
      # chimeric read > 3, too many to deal with
      for (i in c) {
	print c[i],"rt:A:8";
      }
      count_abnorm++;
    }
    else if (count == 2) {
      # code here should be same as above, but it's a "normal" read
      j = 0;
      for (i in c) {
	split(c[i], tmp);
	split(tmp[1],readname,"/");
	str[j] = and(tmp[2],16);
	notprimary[j] = and(tmp[2],256);
	chr[j] = tmp[3];
	# get rid of "_" in chromosome name, for cb_str
	gsub(/_/, "", chr[j]);
	pos[j] = tmp[4];
	m[j] = tmp[5];
	cigarstr[j] = tmp[6];
	seq[j] = tmp[10];
	name[j] = tmp[1];
	
        # blacklist - if 3rd bit set (=4) it means unmapped
        mapped[j] = and(tmp[2],4) == 0; 
	
	if (str[j] == 0 && tmp[6] ~/^[0-9]+S/ && notprimary[j] == 0) {
	  split(tmp[6], cigar, "S");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	else if (str[j] == 0 && tmp[6] ~/^[0-9]+H/ && notprimary[j] == 0) {
	  split(tmp[6], cigar, "H");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	else if (str[j] == 16) {
	  # count Ms,Ds,Ns,Xs,=s for sequence length 
	  seqlength=0; 
	  currstr=tmp[6];
	  where=match(currstr, /[0-9]+[M|D|N|X|=]/); 
	  while (where>0) {
	    seqlength+=substr(currstr,where,RLENGTH-1)+0;
	    currstr=substr(currstr,where+RLENGTH);
	    where=match(currstr, /[0-9]+[M|D|N|X|=]/);
	  }
	  pos[j] = pos[j] + seqlength - 1;
	  if (tmp[6] ~ /[0-9]+S$/ && notprimary[j] == 0) {
	    where = match(tmp[6],/[0-9]+S$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
	  if (tmp[6] ~ /[0-9]+H$/ && notprimary[j] == 0) {
	    where = match(tmp[6],/[0-9]+H$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
	}
	j++;
      }
      if (mapped[0] && mapped[1]) {
	count_reg++;
	interiorpos1=adjust(pos[0],str[0],cigarstr[0],notprimary[0]);
	interiorpos2=adjust(pos[1],str[1],cigarstr[1],notprimary[1]);

	if (length(site_file) > 0) {
	    frag1 = bsearch(chromosomes[chr[0]], length(chromosomes[chr[0]]),pos[0]);
	    frag2 = bsearch(chromosomes[chr[1]], length(chromosomes[chr[1]]),pos[1]);

	    # this happens with circular chromosomes (MT)
	    # only adjusting here since it matters for internal position
	    if (frag1 >= length(chromosomes[chr[0]])) 
		frag1 = length(chromosomes[chr[0]])-1;
	    if (frag2 >= length(chromosomes[chr[1]])) 
		frag2 = length(chromosomes[chr[1]])-1;

	    # adjust internal position based on cutting site
	    if (str[0] == 0) {
		interiorpos1 = chromosomes[chr[0]][frag1];
	    }
	    else {
		if (frag1 == 0) {
		    interiorpos1 = 0;
		}
		else {
		    interiorpos1 = chromosomes[chr[0]][frag1-1];
		}
	    }
	    if (str[1] == 0) {
		interiorpos2 = chromosomes[chr[1]][frag2];
	    }
	    else {
		if (frag2 == 0) {
		    interiorpos2 = 0;
		}
		else {
		    interiorpos2 = chromosomes[chr[1]][frag2-1];
		}
	    }
	}
	if (less_than(str[0],chr[0],pos[0],str[1],chr[1],pos[1])) {
#	  print str[0],chr[0],interiorpos1,str[1],chr[1],interiorpos2,m[0],cigarstr[0],seq[0],m[1],cigarstr[1],seq[1],name[0],"0$",pos[0],pos[1] > fname1;

	    if (length(site_file) > 0) {
		fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag1,frag2);
	    }

	    # chromosome block string to sort on
	    externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[0],pos[1]);
	    cb_str = "cb:Z:"chr[0]"_"chr[1]"_"fragstr"_"str[0]"_"str[1]"_"externalpos;
	    
	    # assign mate mapping quality, "read type", "interior position", "mate interior position"
	    print c[1],"MQ:i:"m[1],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:0",cb_str;
	    print c[2],"MQ:i:"m[0],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:1",cb_str;
	    
	    if (str[0]==0&&str[1]==16&&chr[0]==chr[1]&&pos[0]<pos[1]&&pos[1]-pos[0]<1000) {
		innerpairs+=1;
		insertsizesum+=abs(pos[1]-pos[0])
	    }
	}
	else {
#	  print str[1],chr[1],interiorpos2,str[0],chr[0],interiorpos1,m[1],cigarstr[1],seq[1],m[0],cigarstr[0],seq[0],name[1],"0$",pos[1],pos[0] > fname1;

	  if (length(site_file) > 0) {
	      fragstr = sprintf("%0" fraglen "d_%0" fraglen "d",frag2,frag1);
	  }

	  # chromosome block string to sort on
	  externalpos=sprintf("%0" chrlen "d_%0" chrlen "d",pos[1],pos[0]);
	  cb_str = "cb:Z:"chr[1]"_"chr[0]"_"fragstr"_"str[1]"_"str[0]"_"externalpos;
	      
	  # assign mate mapping quality, "read type", "interior position", "mate interior position"
	  print c[2],"MQ:i:"m[0],"ip:i:"interiorpos2,"mp:i:"interiorpos1,"rt:A:0",cb_str;
	  print c[1],"MQ:i:"m[1],"ip:i:"interiorpos1,"mp:i:"interiorpos2,"rt:A:1",cb_str;

	  if (str[1]==0&&str[0]==16&&chr[0]==chr[1]&&pos[1]<pos[0]&&pos[0]-pos[1]<1000) {
		innerpairs+=1;
		insertsizesum+=abs(pos[0]-pos[1]);
	  }
	}
      }
      else {
	for (i in c) {
	  print c[i];
	}	
	count_unmapped++;
      }
    }
    else if (count == 1) {
      # this actually shouldn't happen, but it happens with alternate aligners on occasion
      count_abnorm++;
      for (i in c) {
	print c[i],"rt:A:8";
      }
    }
    resfile=fname1".res.txt";
    if (innerpairs == 0) innerpairs=1;
    avginsertsize=insertsizesum/innerpairs;
    printf("%d %d %d %d %d %f\n", tottot, count_unmapped, count_reg, count_norm, count_abnorm, avginsertsize) >> resfile;
}
