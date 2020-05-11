# Returns absolute value of v
function abs(v) {
  return v<0?-v:v;
}
# Examines locs to see if they are within wobble
# locs should be pairs of positions to examine
# If all are within wobble limit, they are "tooclose" and are dups
function tooclose(locs) {
  cls=1;
  # abs function checks that distance is within wobble
  # AND function will require all distances to be within wobble
  for (i=1; i<length(locs); i+=2) {
    cls=and(cls,(abs(locs[i]-locs[i+1]) <= wobble));
  }
  return cls;
}
BEGIN{
  size=0;
  wobble=4;
  dupname=name"collisions_dups.txt";
  nodupname=name"collisions_nodups.txt";
  prevNF=0;
}
NF != 0{
  if (NF != prevNF) {
    # not a duplicate; handle dups array
    handledupsarray=1;
  }
  else {
    split(prev, p);
    test=1;
    # this code is testing that strand and chromosomes match 
    # NEVA: make sure we are sorting by strand as well .... prob. in GAWK as well as sort
    for (i=2; i<=NF; i=i+7) {
      test=and(test, ($i == p[i]));
      test=and(test, ($(i+1) == p[i+1]));
    }
    # this looks at the first position being within wobble
    test=and(test, (abs($4-p[4])<= wobble));
    handledupsarray=!test;
  }

  if (!handledupsarray) {
    line[size]=$0;
    size++;
  }
  else {
    if (size >= 2) {
      for (j=0; j<size; j++) {
	split(line[j], check1);
	if (!(j in dups)) {
	  for (k=j+1; k<size; k++) {
	    split(line[k], check2);
	    ind=1;
	    # set up locs array to check distance
	    for (x=4; x<=length(check1); x=x+7) {
	      locs[ind]=check1[x];
	      locs[ind+1]=check2[x];
	      ind=ind+2;
	    }
	    if (tooclose(locs)) {
	      dups[k]++; #places a 1 at dups[k]
	    }
	    delete locs;
	  }
	}
      }
      # print dups out to dup file, non-dups out to non-dup file
      for (j=0; j<size; j++) {
	if (j in dups) {
	  print line[j] > dupname
	}
	else {
	  print line[j] > nodupname
	}
      }
    }
    else { # size=1
      print line[0] > nodupname;
    }
    delete line;
    delete dups;
    size=1;
    line[0]=$0;
  }
  prev=$0;
  prevNF=NF;
}
END {
  if (size >= 2) {
    for (j=0; j<size; j++) {
      split(line[j], check1);
      if (!(j in dups)) {
	for (k=j+1; k<size; k++) {
	  split(line[k], check2);
	  ind=1;
	  # set up locs array to check distance
	  for (x=4; x<=length(check1); x=x+7) {
	    locs[ind]=check1[x];
	    locs[ind+1]=check2[x];
	    ind=ind+2;
	  }
	  if (tooclose(locs)) {
	    dups[k]++; #places a 1 at dups[k]
	  }
	  delete locs;
	}
      }
    }
    # print dups out to dup file, non-dups out to non-dup file
    for (j=0; j<size; j++) {
      if (j in dups) {
	print line[j] > dupname
      }
      else {
	print line[j] > nodupname
      }
    }
  }
  else { # size=1
    print line[0] > nodupname;
  }
}