#!/bin/bash
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
# Diploid script. Parallelizes based on chromosome.  Assumes merged_nodups.txt is
# in current working directory.  
# Takes in either the stem, e.g.
# stem=/broad/aidenlab2/neva/diploid/patski/patski_
# where at that location there are files 
# <stem>_chr1_chr_pos.txt <stem>_chr1_paternal_maternal.txt
# <stem>_chr2_chr_pos.txt <stem>_chr2_paternal_maternal.txt
# etc.
# or a phased VCF file.
# If a phased VCF file is sent, the chr_pos and paternal_maternal files will be created
# in the same directory as the VCF file and the stem will be the first sample (column 10)
#
# See the file diploid.pl for an explanation of the chr_pos and paternal_maternal files
#
# This script also requires a genomeID.  Currently only mouse and human genomes are 
# supported.
#
# Usage: diploid.sh -g genomeID -s stem_or_vcf [-h]
# Juicer version 1.5
juicer_version="1.5" 
# cluster specific settings
usePath=/broad/software/scripts/useuse
load_java="use Java-1.7"
# Juicer directory, contains scripts/
juiceDir="/broad/aidenlab"
# unique name for jobs in this run
groupname="d"`date +%s`

usageHelp="Usage: ${0##*/} -g genomeID -s stem_or_vcf [-h]"
genomeHelp="   genomeID must be defined in the script, e.g. \"hg19\" or \"mm10\""
stemHelp="   stem must be the stem where the chr_pos and paternal_maternal files are; alternatively, the full path to the phased vcf"
helpHelp="   -h: print this help and exit"
printHelpAndExit() {
    echo "$usageHelp"
    echo "$genomeHelp"
    echo "$stemHelp"
    echo "$helpHelp"
    exit $1
}

while getopts "hg:s:" opt; do
    case $opt in
        g) genomeID=$OPTARG ;;
        s) stem=$OPTARG ;;
        h) printHelpAndExit 0 ;;
        [?]) printHelpAndExit 1 ;;
    esac
done

if [ -z $genomeID ] || [ -z $stem ]
then
    printHelpAndExit 1
fi

if [ "$genomeID" == "mm10" ] || [ "$genomeID" == "mm9" ]
then
    chromosomes=(chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY)
elif [ "$genomeID" == "hg19" ] || [ "$genomeID" == "hg38" ]
then
    chromosomes=(1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 X Y)
else
    echo "***! Genome ID $genomeID unrecognized"
    echo "$genomeHelp"
    exit 1
fi

if [ ${stem: -4} == ".vcf" ]
then
    # VCF file, check existence, run phasing
    if [ ! -e ${stem} ]
    then
        echo "***! VCF file ${stem} doesn't exist"
        echo "$stemHelp"
        exit 1
    fi
    holdname1="-hold_jid ${groupname}vcf"
    stem1=$(awk '$0~/^#CHROM/{print $10; exit}' $stem)
    qsub -N ${groupname}vcf -o diploid.out -cwd -j y -r y <<-EOF
      #!/bin/bash
      samplename=\$(awk -f ${juiceDir}/scripts/vcftotxt.awk $stem)
      arr=(\$samplename)
      stem=\${arr[0]}"_"

      for i in ${chromosomes[*]}
      do
        awk -v chr=\$i -v stem=\$stem '{split(\$1,a,":"); if (a[1]==chr){print >> stem""chr"_paternal_maternal.txt"}}' \${stem}paternal_maternal.txt
        awk -v chr=\$i -v stem=\$stem '\$1==chr{print >> stem""chr"_chr_pos.txt"}' \${stem}chr_pos.txt
      done
EOF
   stem=${stem1}"_"
fi     

for i in ${chromosomes[*]}
  do
    if [ -z $holdname ]
    then
        holdname="${groupname}diploid${i}"
        catname="cat diploid_${i}.txt "
    else
        holdname="${holdname},${groupname}diploid${i}"
        catname="${catname} diploid_${i}.txt "
    fi
    qsub -N ${groupname}diploid${i} -l h_vmem=8g -o diploid.out -cwd -j y -r y $holdname1 <<-EOF
awk -v chr=$i '\$2==chr && \$6==chr && \$9 >= 10 && \$12 >= 10' merged_nodups.txt | ${juiceDir}/scripts/diploid.pl -s ${stem}${i}_chr_pos.txt -o ${stem}${i}_paternal_maternal.txt > diploid_${i}.txt
EOF
done

qsub -N ${groupname}catdiploid -o diploid.out -cwd -j y -q long -hold_jid $holdname -r y <<-EOF
# list of all chromosomes. do not continue unless all have been created
  #!/bin/bash
  $catname > diploid.txt
  if [ \$? -ne 0 ]                
    then                              
      echo "***! Failure during concatenation"  
      exit 100
    else
      touch cat.done
  fi
EOF
echo "if [ -e cat.done ]; then awk -f ${juiceDir}/scripts/diploid_split.awk diploid.txt; touch split.done; fi" | qsub -hold_jid ${groupname}catdiploid -o diploid.out -cwd -j y -N ${groupname}splitdiploid -r y

qsub -hold_jid ${groupname}splitdiploid -o diploid.out -cwd -j y -N ${groupname}maternal -r y <<-EOF 
  if [ ! -e split.done ]
  then
     echo "***! Error, split didn't work"
     exit 100
  fi
  source $usePath
  $load_java
  sort -k2,2d -m maternal.txt maternal_both.txt | awk '{split(\$11,a,"/"); print a[1], \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,"100","100"}' > mat.txt 
  ${juiceDir}/juicer_tools pre mat.txt maternal.hic $genomeID
EOF
qsub -hold_jid ${groupname}splitdiploid -o diploid.out -cwd -j y -N ${groupname}paternal -r y <<-EOF 
  if [ ! -e split.done ]
  then
     echo "***! Error, split didn't work"
     exit 100
  fi
  source $usePath
  $load_java
  sort -k2,2d -m paternal.txt paternal_both.txt | awk '{split(\$11,a,"/"); print a[1], \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,"100","100"}' > pat.txt 
  ${juiceDir}/juicer_tools pre pat.txt paternal.hic $genomeID
EOF

qsub -hold_jid ${groupname}paternal,${groupname}maternal -o diploid.out -cwd -j y -r y <<-EOF
  rm *.done
  rm diploid_*
  echo "(-: Diploid pipeline successfully completed"
EOF
