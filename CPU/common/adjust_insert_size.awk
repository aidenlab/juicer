#!/usr/bin/awk -f

# returns absolute value
function abs(value)
{
  return (value<0?-value:value);
}
BEGIN{
  while(getline<avgInsertFile>0) {
    insertsize=$7
  }
  saved="";
  OFS="\t";
}
{
    if ($(NF-1) ~ /^rt:/ && ($(NF-1) ~ /:0$/ || $(NF-1) ~ /:1$/)) {
	# first appearance, process second one
	if (length(saved)==0) {
	    saved=$0;
	    saved_length=length($10);
	    next;
	}
	else {
	    split($NF, cb_str, ":");
	    split(cb_str[3], cb, "_");
	    strand1=cb[5];
	    strand2=cb[6];

	    for (i=12; i<=NF; i++) {
		if ($i ~ /^ip/) {
		    split($i, ip, ":");
		}
		else if ($i ~ /^mp/) {
		    split($i, mp, ":");
		}
	    }

	    is = insertsize - length($10) - saved_length;

	    split($(NF-1), rt, ":");
	    if (rt[3]==0) {
		cigar1=$6;
		split(saved, prev);
		cigar2=prev[6];
		pos1=ip[3];
		pos2=mp[3];
	    }
	    else {
		cigar2=$6;
		split(saved, prev);
		cigar1=prev[6];
		pos1=mp[3];
		pos2=ip[3];
	    }
	    cigloc1=0;
	    cigloc2=0;
	    newpos1=pos1;
	    newpos2=pos2;
	    if (strand1==0) {
		if (cigar1 ~ /[0-9]+S$/) {
		    where = match(cigar1,/[0-9]+S$/);
		    cigloc1 = cigloc1 + (substr(cigar1,where,RLENGTH-1) + 0);
		}
		else if (cigar1 ~ /[0-9]+H$/) {
		    where = match(cigar1,/[0-9]+H$/);
		    cigloc1 = cigloc1 + (substr(cigar1,where,RLENGTH-1) + 0);
		}
	    }
	    else if (strand1==16) {
		if (cigar1 ~ /^[0-9]+S/) {
		    where = match(cigar1,/^[0-9]+S/);
		    cigloc1 = cigloc1 + (substr(cigar1,where,RLENGTH-1) + 0);
		}
		else if (cigar1 ~ /^[0-9]+H/) {
		    where = match(cigar1,/^[0-9]+H/);
		    cigloc1 = cigloc1 + (substr(cigar1,where,RLENGTH-1) + 0);
		}
	    }
	    if (strand2==0) {
		if (cigar2 ~ /[0-9]+S$/) {
		    where = match(cigar2,/[0-9]+S$/);
		    cigloc2 = cigloc2 + (substr(cigar2,where,RLENGTH-1) + 0);
		}
		else if (cigar2 ~ /[0-9]+H$/) {
		    where = match(cigar2,/[0-9]+H$/);
		    cigloc2 = cigloc2 + (substr(cigar2,where,RLENGTH-1) + 0);
		}
	    }
	    else if (strand2==16) {
		if (cigar2 ~ /^[0-9]+S/) {
		    where = match(cigar2,/^[0-9]+S/);
		    cigloc2 = cigloc2 + (substr(cigar2,where,RLENGTH-1) + 0);
		}
		else if (cigar2 ~ /^[0-9]+H/) {
		    where = match(cigar2,/^[0-9]+H/);
		    cigloc2 = cigloc2 + (substr(cigar1,where,RLENGTH-1) + 0);
		}
	    }
	    if (cigloc1==0 && cigloc2==0) {
		if (strand1==0) {
		    newpos1 += int(is/2);
		}
		else if (strand1==16) {
		    newpos1 -= int(is/2);
		}
		if (strand2 ==0) {
		    newpos2 += int(is/2);
		}
		else if (strand2 == 16) {
		    newpos2 -= int(is/2);
		}
	    }
	    else if (cigloc1>0&&cigloc2==0) {
		if (strand2==0) {
		    newpos2 = newpos2 + int(is) + cigloc1;
		}
		else if (strand2==16) {
		    newpos2 = newpos2 - int(is) - cigloc1;
		}
	    }
	    else if (cigloc1==0&&cigloc2>0) {
		if (strand1==0) { 
		    newpos1 = newpos1 + int(is) + cigloc2;
		}
		else if (strand1==16) {
		    newpos1 = newpos1 - int(is) - cigloc2;
		}
	    }

	    if (rt[3]==0) {
		p1 = newpos1;
		p2 = newpos2;
	    }
	    else {
		p1 = newpos2;
		p2 = newpos1;
	    }
	    for (i=12; i<=NF; i++) {
		if ($i ~ /^ip/) {
		    $i = "ip:i:"p1;
		}
		else if ($i ~ /^mp/) {
		    $i = "mp:i:"p2;
		}
	    }
		
	    str=prev[1];
	    for (j=2; j<=length(prev); j++) {
		if (prev[j] ~ /^ip/) {
		    prev[j] = "ip:i:"p2;
		}
		else if (prev[j] ~ /^mp/) {
		    prev[j] = "mp:i:"p1;
		}
		str=str"\t"prev[j];
	    }
	    
	    print str;
	    print;
	    saved="";
	}
    }
    # otherwise, not adjusting the insert size, just print (also header, collisions, etc)
    else {
	print;
    }
}
