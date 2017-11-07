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
# only those reads mapped with MAPQ > 0
#
# Output to fname1 is of the form:
# strand1 chromosome1 position1 frag1 strand2 chromosome2 position2 frag2
#  mapq1 cigar1 seq1 mapq2 cigar2 seq2
#   where chr1 is always <= chr2
# 
# SAMs are also outputted.
#
# Chimeric reads are treated as follows:
# "Normal" chimeras (with 3 reads, two with the ligation junction), where the
# third read (w/o ligation junction) is within 20Kbp of one of the ligation
# junction reads, are sent to fname"_norm".  All other chimeras (where either there
# are more than 3 reads, or they don't contain the ligation junction, or the
# one non-ligation junction end is far from the others, are sent to fname"_abnorm".
#
# awk -f chimeric_blacklist.awk -v fname="stem"
# Juicer version 1.6

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
  count_unmapped = -1; # will count first non-group
  norm=fname"_norm.txt"
  alignable=fname"_alignable.sam"
  collisions=fname"_collisions.sam"
  lowqcollisions=fname"_collisions_low_mapq.sam"
  unmapped=fname"_unmapped.sam"
  mapq0=fname"_mapq0.sam"
  dup=fname"_dup.sam" # for after duplicates have been marked
  linenum=0
}
$0 ~ /^@/{
  # print SAM header to SAM files
  linenum++;
  header = header""$0"\n";
}
$0 !~ /^@/{
  if (tottot == -1) {
    linenum++;
    header = header"@PG\tID:Juicer\tVN:1.6";
  }
  # input file is sorted by read name.  Look at read name to group 
  # appropriately
  # we don't need this with paired end but might as well keep it
  split($1,a,"/"); 

  if(a[1]==prev){
    # move on to next record.  look below this block for actions that occur
    # regardless.
    count++;
  }
  else { 
    # deal with read pair group
    tottot++;
    normal=0; # assume not normal, will get set if normal
    # first take care of unmapped reads
    # blacklist - if 3rd bit set (=4) it means unmapped
    sum_mapped = 0;
    for (j=1; j <= count; j++) {
      split(c[j],tmp);
      mapped[j] = and(tmp[2],4) == 0;
      sum_mapped = sum_mapped+mapped[j];
    }
    if (sum_mapped < 2) {
      if (count_unmapped == -1) {
	print header > unmapped;
      }
      # unmapped
      for (j=1; j <= count; j++) {
	print c[j] > unmapped;
      }
      # otherwise first time through will count extra unmapped read
      count_unmapped++;
    }
    else {
      # find if there are at least two mapq > 0 reads
      sum_mapq=0;
      for (j=1; j <= count; j++) {
        split(c[j],tmp);
	if (tmp[5] > 0) {
          sum_mapq=sum_mapq+1;
	}
      }
      if (sum_mapq < 2) {
	if (count_mapq0 == 0) {
          print header > mapq0;
	}
	count_mapq0++;
        for (j=1; j <= count; j++) {
          print c[j] > mapq0;
	}
      }
      else if (count == 1) {
        # count==1 shouldn't happen; occurs with alternate aligners at times
	# otherwise can signal a problem matching the fastq readnames
	if (count_lowqcollisions == 0) {
	  print header > lowqcollisions;
	}	
	print c[1] > lowqcollisions;
	count_lowqcollisions++;
      }
      else if (count > 4) {	  
	# check that MAPQ >= 10 for all reads
	mapq10=1;
	for (j=1; j <= count; j++) {
	  split(c[j],tmp);
	  if (tmp[5] < 10) {
	    mapq10=0;
	  }
	}
	if (mapq10){
	  if (count_collisions == 0) {
	    print header > collisions;
	  }
	  for (j=1; j <= count; j++) {
	    print c[j] > collisions;
	  }
	  count_collisions++;
	}
	else {
	  if (count_lowqcollisions == 0) {
	    print header > lowqcollisions;
	  }
	  for (j=1; j <= count; j++) {
	    print c[j] > lowqcollisions;
	  }
	  count_lowqcollisions++;
	}
      }
      else {
        # count is 2,3,4.
	# 2 will be normal paired 
	# 3 will be chimeric paired or abnormal
        # (chimeric pair, depends on where the ends are)
        # 4 will be chimeric paired or abnormal
        # (chimeric pair, depends on where the ends are) 
	for (j=1; j <= count; j++) {
	  split(c[j],tmp);
	  # first in pair: 64
	  read[j] = (and(tmp[2], 64) > 0);
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
	  } # else if str[j]==16
	} # end for loop, now have all info about each read end
	read1=0;
	read2=0;
	if (sum_mapq == 2) {
	  # take the mapq positive as the reads
	  for (j=1; j <= count; j++) {
	    if (m[j]>0 && !read1) read1=j;
	    else if (m[j]>0 && read1) read2=j;
	  }
	  normal=1;
	}
	else if (sum_mapq == 4) {
	  # looking for A/B...A/B
	  # will always be first in pair, second in pair
	  if (read[1] == read[2] && read[3] == read[4] && read[1] != read[3]) {
	    dist[13] = abs(chr[1]-chr[3])*10000000 + abs(pos[1]-pos[3]);
	    dist[24] = abs(chr[2]-chr[4])*10000000 + abs(pos[2]-pos[4]);
	    dist[14] = abs(chr[1]-chr[4])*10000000 + abs(pos[1]-pos[4]);
	    dist[23] = abs(chr[2]-chr[3])*10000000 + abs(pos[2]-pos[3]);
		
	    # A/B, A/B
	    if ((dist[13] < 1000 && dist[24] < 1000) || (dist[14] < 1000 && dist[23] < 1000)) {
	      read1 = 1;
	      read2 = 2;
	      normal = 1;
	    }
	    else {
	      # abnormal
	      # check that MAPQ >= 10 for all reads
	      mapq10=1;
	      for (j=1; j <= count; j++) {
		split(c[j],tmp);
		if (tmp[5] < 10) {
		  mapq10=0;
		}
	      }
	      if (mapq10) {
		if (count_collisions == 0) {
		  print header > collisions;
		}
		for (j=1; j <= count; j++) {
		  print c[j] > collisions;
		}
		count_collisions++;
	      }
	      else {
		if (count_lowqcollisions == 0) {
		  print header > lowqcollisions;
		}
		for (j=1; j <= count; j++) {
		  print c[j] > lowqcollisions;
		}
		count_lowqcollisions++;
	      }
	    }
	  }
	  else { # count 4, not close together
	    # abnormal
	    # check that MAPQ >= 10 for all reads
	    mapq10=1;
	    for (j=1; j <= count; j++) {
	      split(c[j],tmp);
	      if (tmp[5] < 10) {
		mapq10=0;
	      }
	    }
	    if (mapq10) {
	      if (count_collisions == 0) {
	        print header > collisions;
	      }
	      for (j=1; j <= count; j++) {
		print c[j] > collisions;
	      }
	      count_collisions++;
	    }
	    else {
	      if (count_lowqcollisions == 0) {
		print header > lowqcollisions;
	      }
	      for (j=1; j <= count; j++) {
		print c[j] > lowqcollisions;
	      }
	      count_lowqcollisions++;
	    }
	  }
	} # if sum_mapq == 4
	else { # sum_mapq == 3
	  r1=0; r2=0; r3=0;
	  for (j=1; j <= count; j++) {
	    if (m[j]>0 && !r1 && !r2) r1=j;
	    else if (m[j]>0 && r1 && !r2) r2=j;
	    else if (m[j]>0 && r1 && r2) r3=j;
	  }
          dist[12] = abs(chr[r1]-chr[r2])*10000000 + abs(pos[r1]-pos[r2]);
	  dist[23] = abs(chr[r2]-chr[r3])*10000000 + abs(pos[r2]-pos[r3]);
	  dist[13] = abs(chr[r1]-chr[r3])*10000000 + abs(pos[r1]-pos[r3]);
    
	  if (min(dist[12],min(dist[23],dist[13])) < 1000) {
	    # The paired ends look like A/B...B. Make sure we take A and B
	    if (read[r1] == read[r2]) {
	      # take the unique one "B" for sure
	      read2 = r3;
	      # take the end of "A/B" that isn't close to "B"
	      read1 = dist[13] > dist[23] ? r1:r2;
	    }
	    else if (read[r1] == read[r3]) { # but this shouldn't happen
	      read2 = r2;
	      read1 = dist[12] > dist[23] ? r1:r3;
	      print "This shouldn't happen",name[read2] 
	    }
	    else if (read[r2] == read[r3]) {
	      read2 = r1;
	      read1 = dist[12] > dist[13] ? r2:r3;
	    }
	    else {
	      print "This shouldn't happen"
	    }
	    normal = 1;
	  }
	  else {
	    # chimeric read with the 3 ends > 1KB apart
	    # check that MAPQ >= 10 for all reads
	    mapq10=1;
	    for (j=1; j <= count; j++) {
	      split(c[j],tmp);
	      if (tmp[5] < 10) {
		mapq10=0;
	      }
	    }
	    if (mapq10) {
	      if (count_collisions == 0) {
	        print header > collisions;
	      }
	      for (j=1; j <= count; j++) {
		print c[j] > collisions;
	      }
	      count_collisions++;
	    }
	    else {
	      if (count_lowqcollisions == 0) {
		print header > lowqcollisions;
	      }
	      for (j=1; j <= count; j++) {
		print c[j] > lowqcollisions;
	      }
	      count_lowqcollisions++;
	    }
	  }
	} # sum_mapq == 3
	if (normal) {
	  # still need to check mapq0
	  if (m[read1]==0 || m[read2]==0) {
	    print "This shouldn't happen"
            # mapq0
	    for (j=1; j <= count; j++) {
	      print c[j] > mapq0;
	    }
	    count_mapq0++;
	  }
	  else {
	    if (count_reg == 0 && count_norm == 0) {
	      print header > alignable;
	    }
	    if (sum_mapped == 2 && count==2) {
	      count_reg++;
	    }
	    else count_norm++;
	    prevlinenum=linenum;
	    for (j=1; j <= count; j++) {
	      print c[j] > alignable;
	      linenum++;
	    }

	    if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
	      print str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2],m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1]"$"fname"$"(prevlinenum+1),name[read2]"$"fname"$"linenum > norm;
	    }
	    else {
	      print str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos[read1],m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2]"$"fname"$"(prevlinenum+1),name[read1]"$"fname"$"linenum > norm;
	    }
	  } # else on check mapq0 again
	} # if normal
      } # else (mapq above 0) 
    } # else (sum_mapped >= 2)
    # reset variables
    delete c;
    delete tmp;
    count=1;
    prev=a[1];
  }
  # these happen no matter what, after the above processing
  c[count] = $0;
}
END{
    # deal with read pair group
    tottot++;
    normal=0; # assume not normal, will get set if normal
    # first take care of unmapped reads
    # blacklist - if 3rd bit set (=4) it means unmapped
    sum_mapped = 0;
    for (j=1; j <= count; j++) {
        split(c[j],tmp);
	mapped[j] = and(tmp[2],4) == 0;
	sum_mapped = sum_mapped+mapped[j];
    }

    if (sum_mapped < 2) {
      if (count_unmapped == -1) {
	print header > unmapped;
      }
      # unmapped
      for (j=1; j <= count; j++) {
	print c[j] > unmapped;
      }
      # otherwise first time through will count extra unmapped read
      count_unmapped++;
    }
    else {
      # find minimum mapq; if it's 0, send to mapq0 file
      sum_mapq=0;
      minmapq=1;
      for (j=1; j <= count; j++) {
        split(c[j],tmp);
	if (tmp[5] > 0) {
          sum_mapq=sum_mapq+1;
	}
	else minmapq=0;
      }
      if (sum_mapq < 2) {
	if (count_mapq0 == 0) {
          print header > mapq0;
	}
	count_mapq0++;
        for (j=1; j <= count; j++) {
          print c[j] > mapq0;
	}
      }
      else if (count == 1) {
        # count==1 shouldn't happen; occurs with alternate aligners at times
	# otherwise can signal a problem matching the fastq readnames
	if (count_lowqcollisions == 0) {
	  print header > lowqcollisions;
	}	
	print c[1] > lowqcollisions;
	count_lowqcollisions++;
      }
      else if (count > 4) {	  
	# check that MAPQ >= 10 for all reads
	mapq10=1;
	for (j=1; j <= count; j++) {
	  split(c[j],tmp);
	  if (tmp[5] < 10) {
	    mapq10=0;
	  }
	}
	if (mapq10){
	  if (count_collisions == 0) {
	    print header > collisions;
	  }
	  for (j=1; j <= count; j++) {
	    print c[j] > collisions;
	  }
	  count_collisions++;
	}
	else {
	  if (count_lowqcollisions == 0) {
	    print header > lowqcollisions;
	  }
	  for (j=1; j <= count; j++) {
	    print c[j] > lowqcollisions;
	  }
	  count_lowqcollisions++;
	}
      }
      else {
        # count is 2,3,4.
	# 2 will be normal paired 
	# 3 will be chimeric paired or abnormal
        # (chimeric pair, depends on where the ends are)
        # 4 will be chimeric paired or abnormal
        # (chimeric pair, depends on where the ends are) 
	for (j=1; j <= count; j++) {
	  split(c[j],tmp);
	  # first in pair: 64
	  read[j] = (and(tmp[2], 64) > 0);
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
	  } # else if str[j]==16
	} # end for loop, now have all info about each read end
	read1=0;
	read2=0;
	if (sum_mapq == 2) {
	  # take the mapq positive as the reads
	  for (j=1; j <= count; j++) {
	    if (m[j]>0 && !read1) read1=j;
	    else if (m[j]>0 && read1) read2=j;
	  }
	  normal=1;
	}
	else if (sum_mapq == 4) {
	  # looking for A/B...A/B
	  # will always be first in pair, second in pair
	  if (read[1] == read[2] && read[3] == read[4] && read[1] != read[3]) {
	    dist[13] = abs(chr[1]-chr[3])*10000000 + abs(pos[1]-pos[3]);
	    dist[24] = abs(chr[2]-chr[4])*10000000 + abs(pos[2]-pos[4]);
	    dist[14] = abs(chr[1]-chr[4])*10000000 + abs(pos[1]-pos[4]);
	    dist[23] = abs(chr[2]-chr[3])*10000000 + abs(pos[2]-pos[3]);
		
	    # A/B, A/B
	    if ((dist[13] < 1000 && dist[24] < 1000) || (dist[14] < 1000 && dist[23] < 1000)) {
	      read1 = 1;
	      read2 = 2;
	      normal = 1;
	    }
	    else {
	      # abnormal
	      # check that MAPQ >= 10 for all reads
	      mapq10=1;
	      for (j=1; j <= count; j++) {
		split(c[j],tmp);
		if (tmp[5] < 10) {
		  mapq10=0;
		}
	      }
	      if (mapq10) {
		if (count_collisions == 0) {
		  print header > collisions;
		}
		for (j=1; j <= count; j++) {
		  print c[j] > collisions;
		}
		count_collisions++;
	      }
	      else {
		if (count_lowqcollisions == 0) {
		  print header > lowqcollisions;
		}
		for (j=1; j <= count; j++) {
		  print c[j] > lowqcollisions;
		}
		count_lowqcollisions++;
	      }
	    }
	  }
	  else { # count 4, not close together
	    # abnormal
	    # check that MAPQ >= 10 for all reads
	    mapq10=1;
	    for (j=1; j <= count; j++) {
	      split(c[j],tmp);
	      if (tmp[5] < 10) {
		mapq10=0;
	      }
	    }
	    if (mapq10) {
	      if (count_collisions == 0) {
	        print header > collisions;
	      }
	      for (j=1; j <= count; j++) {
		print c[j] > collisions;
	      }
	      count_collisions++;
	    }
	    else {
	      if (count_lowqcollisions == 0) {
		print header > lowqcollisions;
	      }
	      for (j=1; j <= count; j++) {
		print c[j] > lowqcollisions;
	      }
	      count_lowqcollisions++;
	    }
	  }
	} # if sum_mapq == 4
	else { # sum_mapq == 3
	  r1=0; r2=0; r3=0;
	  for (j=1; j <= count; j++) {
	    if (m[j]>0 && !r1 && !r2) r1=j;
	    else if (m[j]>0 && r1 && !r2) r2=j;
	    else if (m[j]>0 && r1 && r2) r3=j;
	  }
          dist[12] = abs(chr[r1]-chr[r2])*10000000 + abs(pos[r1]-pos[r2]);
	  dist[23] = abs(chr[r2]-chr[r3])*10000000 + abs(pos[r2]-pos[r3]);
	  dist[13] = abs(chr[r1]-chr[r3])*10000000 + abs(pos[r1]-pos[r3]);
    
	  if (min(dist[12],min(dist[23],dist[13])) < 1000) {
	    # The paired ends look like A/B...B. Make sure we take A and B
	    if (read[r1] == read[r2]) {
	      # take the unique one "B" for sure
	      read2 = r3;
	      # take the end of "A/B" that isn't close to "B"
	      read1 = dist[13] > dist[23] ? r1:r2;
	    }
	    else if (read[r1] == read[r3]) { # but this shouldn't happen
	      read2 = r2;
	      read1 = dist[12] > dist[23] ? r1:r3;
	      print "This shouldn't happen",name[read2] 
	    }
	    else if (read[r2] == read[r3]) {
	      read2 = r1;
	      read1 = dist[12] > dist[13] ? r2:r3;
	    }
	    else {
	      print "This shouldn't happen"
	    }
	    normal = 1;
	  }
	  else {
	    # chimeric read with the 3 ends > 1KB apart
	    # check that MAPQ >= 10 for all reads
	    mapq10=1;
	    for (j=1; j <= count; j++) {
	      split(c[j],tmp);
	      if (tmp[5] < 10) {
		mapq10=0;
	      }
	    }
	    if (mapq10) {
	      if (count_collisions == 0) {
	        print header > collisions;
	      }
	      for (j=1; j <= count; j++) {
		print c[j] > collisions;
	      }
	      count_collisions++;
	    }
	    else {
	      if (count_lowqcollisions == 0) {
		print header > lowqcollisions;
	      }
	      for (j=1; j <= count; j++) {
		print c[j] > lowqcollisions;
	      }
	      count_lowqcollisions++;
	    }
	  }
	} # sum_mapq == 3
	if (normal) {
	  # still need to check mapq0
	  if (m[read1]==0 || m[read2]==0) {
	    print "This shouldn't happen"
            # mapq0
	    for (j=1; j <= count; j++) {
	      print c[j] > mapq0;
	    }
	    count_mapq0++;
	  }
	  else {
	    if (count_reg == 0 && count_norm == 0) {
	      print header > alignable;
	    }
	    if (sum_mapped == 2 && count==2) {
	      count_reg++;
	    }
	    else count_norm++;
	    prevlinenum=linenum;
	    for (j=1; j <= count; j++) {
	      print c[j] > alignable;
	      linenum++;
	    }

	    if (less_than(str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2])) {
	      print str[read1],chr[read1],pos[read1],str[read2],chr[read2],pos[read2],m[read1],cigarstr[read1],seq[read1],m[read2],cigarstr[read2],seq[read2],name[read1]"$"fname"$"(prevlinenum+1),name[read2]"$"fname"$"linenum > norm;
	    }
	    else {
	      print str[read2],chr[read2],pos[read2],str[read1],chr[read1],pos[read1],m[read2],cigarstr[read2],seq[read2],m[read1],cigarstr[read1],seq[read1],name[read2]"$"fname"$"(prevlinenum+1),name[read1]"$"fname"$"linenum > norm;
	    }
	  } # else on check mapq0 again
	} # if normal
      } # else (mapq above 0) 
    } # else (sum_mapped >= 2)
    resfile = norm".res.txt"
    printf("%d %d %d %d %d %d %d\n", tottot, count_unmapped, count_reg, count_norm, count_collisions, count_lowqcollisions, count_mapq0) >> resfile;
}
