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
# Parse VCF file into paternal_maternal and chr_pos
# Note that designation is arbitrary as to which comes first
# We'll name "paternal" the first allele and "maternal" the second
#
# Header of VCF
#1. #CHROM
#2. POS
#3. ID
#4. REF
#5. ALT
#6. QUAL
#7. FILTER
#8. INFO
# if GT info present, 9. FORMAT
# 10-end samples
BEGIN {
    OFS="\t";
}
$1 ~ /^#CHROM/ {
    # header line
    if (NF <= 8) {
        print "No genotype information available";
        exit;
    }
    else {
        for (i=10; i<=NF; i++){
            samples[i]=$i;
        }
    }
    prev="notset";
}
$0 !~ /^#/ {
    # grab the sample fields and parse phasing
    for (i=10; i<=NF; i++){
        split($i, a, ":");
        split($5, alt, ",");
        toolong=0;
        for (j=1; j <=length(alt); j++) {
            if (length(alt[j]) > 1) toolong=1;
        }
        # only interested in phased SNPs, not indels
        if (a[1] ~ /\|/ && length($4)==1 && !toolong) {
            # position 0 in phase corresponds to ref
            alt[0]=$4;
            split(a[1], gt, "|");
            # phase is different
            if (gt[1] != gt[2]) {
                # first allele is "paternal" and listed first
                print $1":"$2, alt[gt[1]], alt[gt[2]] >> samples[i]"_paternal_maternal.txt";
                # new chromosome
                if ($1 != prev) {
                    if (prev =="notset") {
                        printf("%s ", $1) >> samples[i]"_chr_pos.txt";
                    }
                    else {
                        printf("\n%s ", $1) >> samples[i]"_chr_pos.txt";
                    }
                    prev=$1;
                }
                printf("%d ", $2) >> samples[i]"_chr_pos.txt";
            }
        }
        else {
            print >> samples[i]"_skipped.vcf";
        }
    }
}
END {
    # add newline to file, if it was created
    if (prev != "notset") {
        for (i in samples) {
            printf("\n") >> samples[i]"_chr_pos.txt";
            print samples[i];
        }
    }
}
