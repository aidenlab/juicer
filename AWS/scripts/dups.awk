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
# Deduping script, checks that positions are within wobble or not
# Juicer version 1.5
# Returns absolute value of v
function abs(v) {
	return v<0?-v:v;
}
# Examines loc1 and loc2 to see if they are within wobble1
# Examines loc3 and loc4 to see if they are within wobble2
# If both are within wobble limit, they are "tooclose" and are dups
function tooclose(loc1,loc2,loc3,loc4) {
	if (abs(loc1-loc2)<=wobble1 && abs(loc3-loc4)<=wobble2) {
		return 1;
	}
	else if (abs(loc3-loc4)<=wobble1 && abs(loc1-loc2)<=wobble2) {
		return 1;
	}
	else return 0;
}

function optcheck(tile1,tile2,x1,x2,y1,y2) {
	if (tile1==tile2 && abs(x1-x2)<50 && abs(y1-y2)<50){
		return 1;
	}
	else return 0;
}
# Executed once, before beginning the file read
BEGIN {
	i=0;
	wobble1=4;
	wobble2=4;
	# names of output files
	# the variable "name" can be set via the -v flag
	dupname=name"dups.txt";
	nodupname=name"merged_nodups.txt";
	optname=name"optdups.txt"; 
}
# strand, chromosome, fragment match previous; first position (sorted) within wobble1

$1 == p1 && $2 == p2 && abs($3-p3)<=wobble1 && $4 == p4 && $5 == p5 && $6 == p6 && $8 == p8  {
	# add to array of potential duplicates
	line[i]=$0;
	pos1[i]=$3;
	pos2[i]=$7;
	n=split($15,readname,":");
  if (n > 1) {
    tile[i]= readname[3] readname[4] readname[5];
    x[i]= readname[6];
    split(readname[7],y_array,"/");
    y[i]= y_array[1];
  }
  else {
    # can't tell since it's not Illumina.
    # make the tiles different so it's not an optical dup
    tile[i]=i;
  }
	i++;
}
# not a duplicate, one of the fields doesn't match
$1 != p1 || $2 != p2 || $4 != p4 || $5 != p5 || $6 != p6 || $8 != p8 || abs($3-p3)>wobble1 {
	# size of potential duplicate array is bigger than 2
	if (i > 2) {
		for (j=0; j<i; j++) {
			# only consider reads that aren't already marked duplicate
			# (no daisy-chaining)
			if (!(j in dups) && !(j in optdups)) {
				for (k=j+1; k<i; k++) {
					# check each item in array against all the rest 
					if (tooclose(pos1[j],pos1[k],pos2[j],pos2[k])) {
						if (optcheck(tile[j],tile[k],x[j],x[k],y[j],y[k])) {
							optdups[k]++;
						}
						else dups[k]++; #places a 1 at dups[k]
					}
				}
			}
		}
		# print dups out to dup file, non-dups out to non-dup file
		for (j=0; j<i; j++) {
			if (j in dups) {
				print line[j] > dupname
			}	
			else if (j in optdups) {
				print line[j] > optname
			}
			else {
				print line[j] > nodupname
			}
		}
	}
	# size of potential duplicate array is 2, just check the two against each other
	else if (i == 2) {
		print line[0] > nodupname;
		if (tooclose(pos1[0],pos1[1],pos2[0],pos2[1])) {
			if  (optcheck(tile[0],tile[1],x[0],x[1],y[0],y[1])){		
				print line[1] > optname
			}
			else {
				print line[1] > dupname;
			}
		}
		else {
			print line[1] > nodupname;
		}
	}
	# size of potential duplicate array is 1, by definition not a duplicate
	else if (i == 1) {
		print line[0] > nodupname;
	}
	# reset all the potential duplicate array variables
	delete line;
	delete dups;
	delete pos1;
	delete pos2;
	delete tile;
	delete x;
	delete y;
	delete optdups;
	i = 1;
	line[0]=$0;
	pos1[0]=$3;
	pos2[0]=$7;
	n = split($15,readname,":");
  if (n > 1) {
	  tile[0]=readname[3] readname[4] readname[5]; 
	  x[0]=readname[6]; 
	  split(readname[7],y_array,"/");
	  y[0]=y_array[1];
  }
  else {
    tile[0]=0;
  }
}

# always reset the fields we're checking against on each read 
{ p1=$1;p2=$2;p3=$3;p4=$4;p5=$5;p6=$6;p7=$7;p8=$8;}
END {
	if (i > 1) {
		# same code as above, just process final potential duplicate array
		for (j=0; j<i; j++) {
			# only consider reads that aren't already marked duplicate
			# (no daisy-chaining)
			if (!(j in dups) && !(j in optdups)) {
				for (k=j+1; k<i; k++) {					
					if (tooclose(pos1[j],pos1[k],pos2[j],pos2[k])) {
						if (optcheck(tile[j],tile[k],x[j],x[k],y[j],y[k])) {
							optdups[k]++;
						}
						else { 
							dups[k]++;
						}
					}
				}
			}
		}
		for (j=0; j<i; j++) {
			if (j in dups) {
				print line[j] > dupname
			}
			else if (j in optdups) {
				print line[j] > optname
			}
			else {
				print line[j] > nodupname
			}
		}
	}
	else if (i == 1) {
		print line[0] > nodupname
	}
}
