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
# Takes a SAM file, looks for chimeric reads and normal reads.  Outputs  
# only those reads mapped to chromosomes 1-24 with MAPQ > 0.
#
# Output to fname1 is of the form:
# strand1 chromosome1 position1 frag1 strand2 chromosome2 position2 frag2
#  mapq1 cigar1 seq1 mapq2 cigar2 seq2
#   where chr1 is always <= chr2
# 
# Output to fname2 retains SAM format.
#
# Chimeric reads are treated as follows:
# "Normal" chimeras (with 3 reads, two with the ligation junction), where the
# third read (w/o ligation junction) is within 20Kbp of one of the ligation
# junction reads, are sent to fname1.  All other chimeras (where either there
# are more than 3 reads, or they don't contain the ligation junction, or the
# one non-ligation junction end is far from the others, are sent to fname2.
#
# awk -f chimeric_blacklist.awk -v fname1="norm_chimera" fname2="abnorm_chimera" fname3="unmapped"
# Juicer version 1.5 

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
BEGIN{
  OFS="\t";
  tottot = -1; # will count first non-group
}
{
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
    if (count==3) {
      # chimeric read
      for (j=1; j <= 3; j++) {
	split(c[j], tmp);
	split(tmp[1],readname,"/");
	read[j] = readname[2];
	name[j] = tmp[1];
	
	# strand; Bit 16 set means reverse strand
	str[j] = and(tmp[2],16);
	# chromosome
	chr[j] = tmp[3];
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
	if (str[j] == 0 && tmp[6] ~/^[0-9]+S/) {
	  split(tmp[6], cigar, "S");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	# get rid of hard clipping to know correct position
	else if (str[j] == 0 && tmp[6] ~/^[0-9]+H/) {
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
	  if (tmp[6] ~ /[0-9]+S$/) {
	    where = match(tmp[6],/[0-9]+S$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
          else if (tmp[6] ~ /[0-9]+H$/) {
	    where = match(tmp[6],/[0-9]+H$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
          }
	  # Mitochrondria loops around
	  if (chr[j] ~ /MT/ && pos[j] >= 16569) {
	    pos[j] = pos[j] - 16569;
	  }
	}
        # blacklist - if 3rd bit set (=4) it means unmapped
        mapped[j] = and(tmp[2],4) == 0; 
      }
			
      dist[12] = abs(chr[1]-chr[2])*10000000 + abs(pos[1]-pos[2]);
      dist[23] = abs(chr[2]-chr[3])*10000000 + abs(pos[2]-pos[3]);
      dist[13] = abs(chr[1]-chr[3])*10000000 + abs(pos[1]-pos[3]);
      
      if (min(dist[12],min(dist[23],dist[13])) < 1000) {
	# The paired ends look like A/B...B.  Make sure we take A and B.
	if (read[1] == read[2]) {
	  # take the unique one "B" for sure
	  read2 = 3;
	  # take the end of "A/B" that isn't close to "B"
	  read1 = dist[13] > dist[23] ? 1:2;
	}
	else if (read[1] == read[3]) {
	  read2 = 2;
	  read1 = dist[12] > dist[23] ? 1:3;
	}
	else if (read[2] == read[3]) {
	  read2 = 1;
	  read1 = dist[12] > dist[13] ? 2:3;
	}
	else {
	  printf("reads strange\n") > "/dev/stderr"
	  exit 1
	}
	
	if (mapped[read1] && mapped[read2]) {
	  count_norm++;
	  if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
	  #  print str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2],m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1],name[read2],qual[read1],qual[read2] > fname1;
	    print str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2],m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1],name[read2] > fname1;
	  }
	  else {
	  #  print str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos[read1],m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2],name[read1],qual[read2],qual[read1] > fname1;
	    print str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos[read1],m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2],name[read1] > fname1;
	  }
	}
	else {
	  for (i in c) {
	    print c[i] > fname3;
	  }	
	  count_unmapped++;
	}
      }
      else {
	# chimeric read with the 3 ends > 1KB apart
	count_abnorm++;
	for (i in c) {
	  print c[i] > fname2;
	}
      }
    }
    else if (count > 3) {
      # chimeric read > 3, too many to deal with
      count_abnorm++;
      for (i in c) {
	print c[i] > fname2;
      }
    }
    else if (count == 2) {
      # code here should be same as above, but it's a "normal" read
      j = 0;
      for (i in c) {
	split(c[i], tmp);
	split(tmp[1],readname,"/");
	str[j] = and(tmp[2],16);
	chr[j] = tmp[3];
	pos[j] = tmp[4];
	m[j] = tmp[5];
	cigarstr[j] = tmp[6];
	seq[j] = tmp[10];
        #qual[j] = tmp[11];
	name[j] = tmp[1];
	
        # blacklist - if 3rd bit set (=4) it means unmapped
        mapped[j] = and(tmp[2],4) == 0; 
	
	if (str[j] == 0 && tmp[6] ~/^[0-9]+S/) {
	  split(tmp[6], cigar, "S");
	  pos[j] = pos[j] - cigar[1];
	  if (pos[j] <= 0) {
	    pos[j] = 1;
	  }
	}
	else if (str[j] == 0 && tmp[6] ~/^[0-9]+H/) {
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
	  if (tmp[6] ~ /[0-9]+S$/) {
	    where = match(tmp[6],/[0-9]+S$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
	  if (tmp[6] ~ /[0-9]+H$/) {
	    where = match(tmp[6],/[0-9]+H$/);
	    cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	    pos[j] = pos[j] + cigloc;
	  }
	  if (chr[j] ~ /MT/ && pos[j] >= 16569) {
	    pos[j] = pos[j] - 16569;
	  }
	}
	j++;
      }
      if (mapped[0] && mapped[1]) {
	count_reg++;
	if (less_than(str[0],chr[0],pos[0],str[1],chr[1],pos[1])) {
	  # ideally we'll get rid of printing out cigar string at some point
	  #print str[0],chr[0],pos[0],str[1],chr[1],pos[1],m[0],cigarstr[0],seq[0],m[1],cigarstr[1],seq[1],name[0],name[1],qual[0],qual[1] > fname1;
	  print str[0],chr[0],pos[0],str[1],chr[1],pos[1],m[0],cigarstr[0],seq[0],m[1],cigarstr[1],seq[1],name[0],name[1] > fname1;
	}
	else {
	  #print str[1],chr[1],pos[1],str[0],chr[0],pos[0],m[1],cigarstr[1],seq[1],m[0],cigarstr[0],seq[0],name[1],name[0],qual[1],qual[0] > fname1;
print str[1],chr[1],pos[1],str[0],chr[0],pos[0],m[1],cigarstr[1],seq[1],m[0],cigarstr[0],seq[0],name[1],name[0] > fname1;
	}
      }
      else {
	for (i in c) {
	  print c[i] > fname3;
	}	
	count_unmapped++;
      }
    }
    else if (count == 1) {
      # this actually shouldn't happen, but it happens with alternate aligners on occasion
      count_abnorm++;
      for (i in c) {
	print c[i] > fname2;
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
  if (count==3) {
    # chimeric read
    for (j=1; j <= 3; j++) {
      split(c[j], tmp);
      split(tmp[1],readname,"/");
      read[j] = readname[2];
      name[j] = tmp[1];
      
      # strand
      str[j] = and(tmp[2],16);
      # chromosome
      chr[j] = tmp[3];
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
      if (str[j] == 0 && tmp[6] ~/^[0-9]+S/) {
	split(tmp[6], cigar, "S");
	pos[j] = pos[j] - cigar[1];
	if (pos[j] <= 0) {
	  pos[j] = 1;
	}
      }
      else if (str[j] == 0 && tmp[6] ~/^[0-9]+H/) {
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
	if (tmp[6] ~ /[0-9]+S$/) {
	  where = match(tmp[6],/[0-9]+S$/);
	  cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	  pos[j] = pos[j] + cigloc;
	}
	if (tmp[6] ~ /[0-9]+H$/) {
	  where = match(tmp[6],/[0-9]+H$/);
	  cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	  pos[j] = pos[j] + cigloc;
	}
	# Mitochrondria loops around
	if (chr[j] ~ /MT/ && pos[j] >= 16569) {
	  pos[j] = pos[j] - 16569;
	}
      }
      mapped[j] = and(tmp[2],4) == 0;
    }
    
    dist[12] = abs(chr[1]-chr[2])*10000000 + abs(pos[1]-pos[2]);
    dist[23] = abs(chr[2]-chr[3])*10000000 + abs(pos[2]-pos[3]);
    dist[13] = abs(chr[1]-chr[3])*10000000 + abs(pos[1]-pos[3]);
    
    if (min(dist[12],min(dist[23],dist[13])) < 1000) {
      # The paired ends look like A/B...B.  Make sure we take A and B.
      if (read[1] == read[2]) {
	# take the unique one "B" for sure
	read2 = 3;
	# take the end of "A/B" that isn't close to "B"
	read1 = dist[13] > dist[23] ? 1:2;
      }
      else if (read[1] == read[3]) {
	read2 = 2;
	read1 = dist[12] > dist[23] ? 1:3;
      }
      else if (read[2] == read[3]) {
	read2 = 1;
	read1 = dist[12] > dist[13] ? 2:3;
      }
      else {
	printf("reads strange\n") > "/dev/stderr"
	exit 1
      }
      
      if (mapped[read1] && mapped[read2]) {
	count_norm++;
	if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
	  print str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2],m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1],name[read2] > fname1;
	}
	else {
	  print str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos[read1],m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2],name[read1] > fname1;
	}
      }
      else {
	for (i in c) {
	  print c[i] > fname3;
	}	
	count_unmapped++;
      }
    }
    else {
      # chimeric read with the 3 ends > 1KB apart
      count_abnorm++;
      for (i in c) {
	print c[i] > fname2;
      }
    }
  }
  else if (count > 3) {
    # chimeric read > 3, too many to deal with
    count_abnorm++;
    for (i in c) {
      print c[i] > fname2;
    }
  }
  else if (count == 2) {
    # code here should be same as above, but it's a "normal" read
    j = 0;
    for (i in c) {
      split(c[i], tmp);
      split(tmp[1],readname,"/");
      str[j] = and(tmp[2],16);
      chr[j] = tmp[3];
      pos[j] = tmp[4];
      m[j] = tmp[5];
      cigarstr[j] = tmp[6];
      seq[j] = tmp[10];
      #qual[j] = tmp[11];
      name[j] = tmp[1];    
      
      if (str[j] == 0 && tmp[6] ~/^[0-9]+S/) {
	split(tmp[6], cigar, "S");
	pos[j] = pos[j] - cigar[1];
	if (pos[j] <= 0) {
	  pos[j] = 1;
	}
      }
      else if (str[j] == 0 && tmp[6] ~/^[0-9]+H/) {
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
	if (tmp[6] ~ /[0-9]+S$/) {
	  where = match(tmp[6],/[0-9]+S$/);
	  cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	  pos[j] = pos[j] + cigloc;
	}
	if (tmp[6] ~ /[0-9]+H$/) {
	  where = match(tmp[6],/[0-9]+H$/);
	  cigloc = substr(tmp[6],where,RLENGTH-1) + 0;
	  pos[j] = pos[j] + cigloc;
	}
	if (chr[j] ~ /MT/ && pos[j] >= 16569) {
	  pos[j] = pos[j] - 16569;
	}
      }
      mapped[j] = and(tmp[2],4) == 0;
      j++;
    }
    if (mapped[0] && mapped[1]) {
      count_reg++;
      if (less_than(str[0],chr[0],pos[0],str[1],chr[1],pos[1])) {
	# ideally we'll get rid of printing out cigar string at some point
	# qual
	print str[0],chr[0],pos[0],str[1],chr[1],pos[1],m[0],cigarstr[0],seq[0],m[1],cigarstr[1],seq[1],name[0],name[1] > fname1;
      }
      else {
	print str[1],chr[1],pos[1],str[0],chr[0],pos[0],m[1],cigarstr[1],seq[1],m[0],cigarstr[0],seq[0],name[1],name[0] > fname1;
      }
    }
    else {
      for (i in c) {
	print c[i] > fname3;
      }	
      count_unmapped++;
    }
  }
  resfile=fname1".res.txt";
  printf("%d %d %d %d %d\n", tottot, count_unmapped, count_reg, count_norm, count_abnorm) >> resfile;
}
