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
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##########

# Dedup script that submits deduping jobs after splitting at known
# non-duplicates
# Juicer version 2.0

BEGIN{
  tot=0;
  name=0;
}
{
  if (tot >= 1000000 && $0 ~/cb:/) { 
    for (ind=12; ind<=NF; ind++) {
      if ($(ind) ~ /^cb:/) {
        split($(ind), cb_str, ":");
      }
    }
    split(cb_str[3], cb, "_");
    if (p1 != cb[1] || p2 != cb[2] || p3 != cb[3] || p4 != cb[4] || p5 != cb[5] || p6 != cb[6]) {
      print "this is the beginging part of split rmdups awk";
      cmd=sprintf("qstat | grep osplit%s |cut -c 1-8", groupname);
      cmd |& getline jID_osplit;
      sname=sprintf("%s_msplit%04d_", groupname, name);
      if (justexact) {
	sysstring1=sprintf("qsub -l %s -o %s_%s.log -j oe -q %s -N DDuP%s%s -l mem=4gb -l nodes=1:ppn=1 -W depend=afterok:%s <<-EOF\nawk -f  %s/scripts/dups_sam.awk -v nowobble=1 -v fname=%s/%s_count %s/split%04d >  %s/%s;\nEOF\n", walltime, outfile, sname, queue, name, groupname, jID_osplit, juicedir, dir, sname, dir, name, dir, name, dir, sname);
      }
      else {
	sysstring1=sprintf("qsub -l %s -o %s_%s.log -j oe -q %s -N DDuP%s%s -l mem=4gb -l nodes=1:ppn=1 -W depend=afterok:%s <<-EOF\nawk -f  %s/scripts/dups_sam.awk -v fname=%s/%s_count %s/split%04d >  %s/%s;\nEOF\n", walltime, outfile, sname, queue, name, groupname, jID_osplit ,juicedir, dir, sname, dir, name, dir, name, dir, sname);
      }
      system(sysstring1);
      cmd1=sprintf("qstat | grep DDuP%s%s | cut -d ' ' -f 1 | cut -d '.' -f 1-2",name,groupname);
      cmd1 |& getline jID;
      
      close(cmd1);
      if (name==0){
	waitstring=sprintf("%s",jID );
	waitstring2=sprintf("%s",jID)
      }
      else {
	waitstring=sprintf("%s:%s",waitstring,jID);
	waitstring2=sprintf("%s %s",waitstring2,jID)
      }
      
      name++;
      tot=0;
    }
  }
  outname = sprintf("%s/split%04d", dir, name);
  print > outname;
  if ( $0 ~ /cb:/) {
    for (ind=12; ind<=NF; ind++) {
      if ($(ind) ~ /^cb:/) {
        split($(ind), cb_str, ":");
      }
    }
    split(cb_str[3], cb, "_");
    p1=cb[1];p2=cb[2];p3=cb[3];p4=cb[4];p5=cb[5];p6=cb[6];
  }
  tot++;
}
END {
  sname=sprintf("%s_msplit%04d_", groupname, name);
  print "This is the begining of END step, submitting _dedup job in END step";
  if (justexact) {
    sysstring1=sprintf("qsub -l %s -o %s_%s.log -j oe -q %s -N DDuP%s%s -l mem=4gb -l nodes=1:ppn=1 -W depend=afterok:%s <<-EOF\nawk -f  %s/scripts/dups_sam.awk -v nowobble=1 -v fname=%s/%s_count %s/split%04d >  %s/%s;\nEOF\n", walltime, outfile, sname, queue, name, groupname, jID_osplit, juicedir, dir, sname, dir, name, dir, name, dir, sname);
  }
  else {
    sysstring1=sprintf("qsub -l %s -o %s_%s.log -j oe -q %s -N DDuP%s%s -l mem=4gb -l nodes=1:ppn=1 -W depend=afterok:%s <<-EOF\nawk -f  %s/scripts/dups_sam.awk -v fname=%s/%s_count %s/split%04d >  %s/%s;\nEOF\n", walltime, outfile, sname, queue, name, groupname, jID_osplit ,juicedir, dir, sname, dir, name, dir, name, dir, sname);
  }
  
  system(sysstring2);
  sname=sprintf("%s_msplit%04d_", groupname, name);

  cmd2=sprintf("qstat | grep DDuP%s%s | cut -d ' ' -f 1 | cut -d '.' -f 1-2",name,groupname);
  cmd2 |& getline jID;
  close(cmd2);
  if ( name==0 ){
    waitstring=sprintf("%s",jID );
    waitstring2=sprintf("%s",jID);
  }
  else {
    waitstring=sprintf("%s:%s",waitstring,jID);
    waitstring2=sprintf("%s %s",waitstring2,jID);
  }
  cmd3=sprintf("qstat %s",waitstring2);
  cmd3 |& getline waitedjobstatus;
  close(cmd3)
  print "below is the waitingstring1 and waitingstring2";
  print waitstring;
  sysstring = sprintf("qsub -l %s -o %s_catsplit.log -j oe -q %s -N CtSplt%s -W depend=afterok:%s <<-EOF\ncat %s/%s_msplit* >  %s/merged_dedup.sam;  nEOF\n", walltime, outfile, queue, groupname, waitstring, dir, groupname, dir);
  print waitstring2;
  print "below is the waited job status before qsub"
  print waitedjobstatus;
  print "submitting _catsplit jobs";
  system(sysstring);
  cmd4=sprintf("qstat | grep CtSplt%s | cut -d ' ' -f 1 | cut -d '.' -f 1-2",groupname);
  cmd4 |& getline jID_catsplit;
  close(cmd4);
  print "below is the _catsplit id";
  print jID_catsplit;
  sysstring=sprintf("qsub -l %s -o %s_rmsplit.log -j oe -q %s -N RmSplt%s -W depend=afterok:%s <<-EOF\n rm %s/*_msplit*;  rm %s/split*;\nEOF", walltime, outfile, queue, groupname,jID_catsplit, dir, dir);
  print "submitting _rmsplit jobs";
  system(sysstring);
  cmd=sprintf("qstat | grep RmSplt%s | cut -d ' ' -f 1 | cut -d '.' -f 1-2",groupname);
  cmd |& getline jID_rmsplit;
  print "below is the jID_rmsplit:";
  print jID_rmsplit;
}
