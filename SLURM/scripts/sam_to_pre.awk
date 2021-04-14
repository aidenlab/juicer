#!/usr/bin/awk -f
BEGIN{
  if (length(mapq)==0){
    mapq=-1;
  }
  if (length(mnd)==0){
    use_external_pos=0;
  }
  else {
    use_external_pos=1;
  }
}
{
    processme=0;
}
$0 ~ /rt:/ {
    for (ind=12; ind<=NF; ind++) {
	if ($(ind) ~ /^rt:/) {
	    split($(ind), rt_str, ":");
	    if (int(rt_str[3]) <= 5) processme=1;
	} 
    }
}
processme != 0 {
    if (length(saved_seq)==0) {
	saved_seq = $10;
	saved_mq = $5;
	saved_cigar = $6;
	saved_chr = $3;
	next;
    }
    else {
	for (ind=12; ind<=NF; ind++) {
	    if ($(ind) ~ /^cb:/) {
		split($(ind), cb_str, ":");
	    }
	}

	split(cb_str[3], cb, "_");

	# fragment and strand from cb string
	# we slightly modify chr name, so take directly from SAM
	frag1 = int(cb[3]);
	frag2 = int(cb[4]);
	str1 = cb[5];
	str2 = cb[6];
	extpos1 = int(cb[7]);
	extpos2 = int(cb[8]);
	
	for (i=12; i<=NF; i++) {
	    if ($i ~ /^ip/) {
		split($i, ip, ":");
	    }
	    else if ($i ~ /^mp/) {
		split($i, mp, ":");
	    }
	}	
	for (ind=12; ind<=NF; ind++) {
	    if ($(ind) ~ /^rt:/) {
		split($(ind), rt, ":");
	    }
	}
	if (rt[3]%2==0) {
	    chr1 = $3;
	    chr2 = saved_chr;
	    pos1 = ip[3];
	    pos2 = mp[3];
	    mq1 = $5;
	    mq2 = saved_mq;
	    cigar1 = $6;
	    cigar2 = saved_cigar;
	    seq1 = $10;
	    seq2 = saved_seq;
	}
	else {
	    chr1 = saved_chr;
	    chr2 = $3;
	    pos1 = mp[3];
	    pos2 = ip[3];
	    mq1 = saved_mq;
	    mq2 = $5;
	    cigar1 = saved_cigar;
	    cigar2 = $6;
	    seq1 = saved_seq;
	    seq2 = $10;
	}

	# merged_nodups
	if (use_external_pos) {
	    pos1 = extpos1;
	    pos2 = extpos2;
	    print str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, mq1, cigar1, seq1, mq2, cigar2, seq2, $1, $1;
	}
	else {
	  if (mq1 >= mapq && mq2 >= mapq) {
	    print str1, chr1, pos1, frag1, str2, chr2, pos2, frag2;
	  }
	}
	saved_seq="";
    }
}
