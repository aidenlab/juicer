#!/usr/bin/perl 
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
# Perl script to convert to fragment map from infile. The infile should be the
# "intermediate" form: no duplicates, 10 fields, laid out as:
#
# strand1 chr1 pos1 strand2 chr2 pos2 mapq1 sequence1 mapq2 sequence2 
#
# The script also requires a restriction site file, which lists on 
# each line, the sorted locations of the enzyme restriction sites.
#
# Usage:  fragment.pl <infile>
# Juicer 1.5
use POSIX;

# Check arguments
if (scalar(@ARGV) == 3) {
  ($infile,$outfile,$site_file) = @ARGV;
}
else {
  print "Usage: fragment.pl <infile> <outfile> <site file>\n";
  print " <infile>: file in intermediate format to calculate statistics on\n";
  print " <outfile>: output, results of fragment search\n";  
  print " <site file>: list of restriction sites, one line per chromosome\n";
  exit;
}
# Global variables for calculating statistics
my %chromosomes;
my %hindIII;

# read in restriction site file and store as multidimensional array
open FILE, $site_file or die $!;

while (<FILE>) {
  my @locs = split;
  my $key = shift(@locs);
  my $ref = \@locs;
  $chromosomes{$key} = $ref;
	if ($key == "14") {
		$chromosomes{$key."m"} = $ref;  
		$chromosomes{$key."p"} = $ref;
	}
}
close(FILE);

# read in infile and calculate statistics
open INFILE, $infile or die $!;
open OUTFILE,">", $outfile or die $!;

while (<INFILE>) {
  my @record = split;

  # find upper index of position in sites array via binary search
  my $index1 = &bsearch($record[2],$chromosomes{$record[1]});
  my $index2 = &bsearch($record[5],$chromosomes{$record[4]});
  print OUTFILE $record[0] . " " . $record[1] . " " . $record[2] . " " . $index1 . " " . $record[3] . " " . $record[4] . " " . $record[5] . " " . $index2 . " ";
	for (my $i=6; $i < scalar(@record); $i++){
		print OUTFILE $record[$i] . " ";
	}
	print OUTFILE "\n";
}
close(INFILE);
close(OUTFILE);

# Binary search, array passed by reference
# search array of integers a for given integer x
# return index where found or upper index if not found
sub bsearch {
  my ($x, $a) = @_;            # search for x in array a
  my ($l, $u) = (0, @$a - 1);  # lower, upper end of search interval
  my $i;                       # index of probe
  while ($l <= $u) {
    $i = int(($l + $u)/2);
    if ($a->[$i] < $x) {
      $l = $i+1;
    }
    elsif ($a->[$i] > $x) {
      $u = $i-1;
    } 
    else {
      return $i+1; # found
    }
  }
  return $l;         
}
