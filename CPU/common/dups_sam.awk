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
# Operates on SAM files with rt and cb flags set
# Juicer version 2.0
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
# Marks the line as duplicate
function markduplicate(input) {
    split(input, output);
    output[2]=or(output[2], 1024);
    str=output[1];
    for (myi=2; myi<=length(output); myi++) {
	str=str"\t"output[myi];
    } 
    print str;
}
# Executed once, before beginning the file read
BEGIN {
    i=0;
    if (length(nowobble)>0) {
	wobble1=0;
	wobble2=0;
    }
    else if (length(wobble1)==0 && length(wobble2)==0) {
	wobble1=4;
	wobble2=4;
    }
    OFS="\t";
    pname="";
}
# header
$0 ~ /^@/{print}
# read we're not considering (chimeric ambiguous or unmapped)
$0 !~ /^@/ && $0 !~ /cb:/{print}
# read we consider
$0 !~ /^@/ && $0 ~ /cb:/ && pname == $1 {
    # save lines associated with readname
    # maximum possible is 4 since maximum chimeric count=4
    # saved1 should already have first line in it
    if (!($1 in saved2)){
	saved2[$1]=$0;
    }
    else if (!($1 in saved3)){ 
	saved3[$1]=$0;
    }
    else {
	saved4[$1]=$0;
    }
}
# handle previous read, this is a new read group
$0 !~ /^@/ && $0 ~ /cb:/ && pname != $1 {
    split(saved1[pname],line);
    # the below happens once, when nothing is saved yet
    if (length(line) == 0) {
	pname=$1;
	saved1[$1]=$0;
	rname[0]=$1;
	next;
    }
    setcb=0;
    for (ind=12; ind<=length(line); ind++) {
      if (line[ind] ~ /^cb:/) {
	split(line[ind], cb_str, ":");
	setcb=1;
      }
    }
    if (!setcb) {
      print "!!Error!! cb field not set" > "/dev/stderr";
      exit;
    }
    split(cb_str[3], cb, "_");
    # check for potential duplicate: chromosome, fragment, strand  match, first position (sorted) within wobble1 
    if (cb[1] == p1 && cb[2] == p2 && int(cb[3]) == p3 && int(cb[4]) == p4 && cb[5] == p5 && cb[6] == p6 && abs(int(cb[7])-p7)<=wobble1) {
	rname[i]=line[1];
	pos1[i]=int(cb[7]);
	pos2[i]=int(cb[8]);
	i++;
    }
    # not a duplicate, one of the fields doesn't match
    else {
	# size of potential duplicate array is bigger than 1
	if (i > 1) {
	    for (j=0; j<i; j++) {
		# only consider reads that aren't already marked duplicate
		# (no daisy-chaining)
		if (!(j in dups)) {
		    for (k=j+1; k<i; k++) {
			# check each item in array against all the rest 
			if (tooclose(pos1[j],pos1[k],pos2[j],pos2[k])) {
			    dups[k]++; #places a 1 at dups[k]
			}
			if (abs(pos1[j]-pos1[k])>wobble1) {
			    break
			}
		    }
		}
	    }
	    # markdups 
	    for (j=0; j<i; j++) {
		if (j in dups) {
		    count_dups++;
		    if (rname[j] in saved1) markduplicate(saved1[rname[j]]);
		    if (rname[j] in saved2) markduplicate(saved2[rname[j]]);
		    if (rname[j] in saved3) markduplicate(saved3[rname[j]]);
		    if (rname[j] in saved4) markduplicate(saved4[rname[j]]);
		}
		else {
		  if (rname[j] in saved1) print saved1[rname[j]];
		  if (rname[j] in saved2) print saved2[rname[j]];
		  if (rname[j] in saved3) print saved3[rname[j]];
		  if (rname[j] in saved4) print saved4[rname[j]];
		}
		delete saved1[rname[j]];
		delete saved2[rname[j]];
		delete saved3[rname[j]];
		delete saved4[rname[j]];
	    }
	}
	# size of potential duplicate array is 1, by definition not a duplicate
	else if (i==1) {
	    if (rname[0] in saved1) print saved1[rname[0]];
	    if (rname[0] in saved2) print saved2[rname[0]];
	    if (rname[0] in saved3) print saved3[rname[0]];
	    if (rname[0] in saved4) print saved4[rname[0]];
	    delete saved1[rname[0]];
	    delete saved2[rname[0]];
	    delete saved3[rname[0]];
	    delete saved4[rname[0]];
	}
	# reset all the potential duplicate array variables
	delete rname;
	delete dups;
	delete pos1;
	delete pos2;
	# reset beginning of potential dups array with current line
	# will be marked not a duplicate
    
	i = 1;
	rname[0]=line[1];
	pos1[0]=int(cb[7]);
	pos2[0]=int(cb[8]);
    } # else not a duplicate
    saved1[$1]=$0;
    pname=$1;
    p1=cb[1]; p2=cb[2]; p3=int(cb[3]); p4=int(cb[4]); p5=cb[5]; p6=cb[6]; p7=int(cb[7]); p8=int(cb[8]); 
    pname=$1;
    saved1[$1]=$0;
} # pname != $1
END {
    split(saved1[pname],line);
    # length of line is 0 if all the reads did not have cb tag set
    # this happens when processing many reads, the first reads in sorted
    # cb order do not have the cb tag set since they are chimeric ambiguous
    if (length(line) != 0) {
	for (ind=12; ind<=length(line); ind++) {
	    if (line[ind] ~ /^cb:/) {
		split(line[ind], cb_str, ":");
		setcb=1;
	    }
	}
	if (!setcb) {
	    print "!!Error!! cb field not set" > "/dev/stderr";
	    exit;
	}
	split(cb_str[3], cb, "_");
	rname[i]=line[1];
	pos1[i]=int(cb[7]);
	pos2[i]=int(cb[8]);
	i++;
	
	# size of potential duplicate array is bigger than 1
	if (i > 1) {
	    for (j=0; j<i; j++) {
		# only consider reads that aren't already marked duplicate
		# (no daisy-chaining)
		if (!(j in dups)) {
		    for (k=j+1; k<i; k++) {
			# check each item in array against all the rest 
			if (tooclose(pos1[j],pos1[k],pos2[j],pos2[k])) {
			    dups[k]++; #places a 1 at dups[k]
			}
			if (abs(pos1[j]-pos1[k])>wobble1) {
			    break
			}
		    }
		}
	    }
	    # markdups 
	    for (j=0; j<i; j++) {
		if (j in dups) {
		    count_dups++;
		    if (rname[j] in saved1) markduplicate(saved1[rname[j]]);
		    if (rname[j] in saved2) markduplicate(saved2[rname[j]]);
		    if (rname[j] in saved3) markduplicate(saved3[rname[j]]);
		    if (rname[j] in saved4) markduplicate(saved4[rname[j]]);
		}
		else {
		    if (rname[j] in saved1) print saved1[rname[j]];
		    if (rname[j] in saved2) print saved2[rname[j]];
		    if (rname[j] in saved3) print saved3[rname[j]];
		    if (rname[j] in saved4) print saved4[rname[j]];
		}
	    }
	}
	# size of potential duplicate array is 1, by definition not a duplicate
	else if (i==1) {
	    if (rname[0] in saved1) print saved1[rname[0]];
	    if (rname[0] in saved2) print saved2[rname[0]];
	    if (rname[0] in saved3) print saved3[rname[0]];
	    if (rname[0] in saved4) print saved4[rname[0]];
	}
	#print count_dups > fname;
    }
}

