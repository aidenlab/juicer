#!/bin/bash

parallel -a $1 --pipepart --will-cite --block 1G --jobs 80% "awk '{if (\$6==\"CG\"){c++}else{h++}}\$4+\$5>0{if (\$6==\"CG\"){if (\$4<\$5){cunmeth++}else{cmeth++} }else { if (\$4<\$5){unmeth++}else{meth++}} }END{print cunmeth, cmeth, unmeth,meth,c,h}' >> $1.tmp"

awk 'BEGIN{cunmeth=0;cmeth=0;unmeth=0;meth=0; }{cunmeth+=$1; cmeth+=$2; unmeth+=$3; meth+=$4; c+=$5; h+=$6}END{print "Cytosine conversion report"; printf("C->T conversion: %0.2f%\n",(unmeth*100)/(unmeth+meth+1)); printf("CHH+CHG coverage: %0.2f%\n", ((unmeth+meth)*100)/h); printf("C->T conversion CpGs: %0.2f%\n",(cunmeth*100)/(cunmeth+cmeth+1)); printf("CpG coverage: %0.2f%\n", ((cunmeth+cmeth)*100)/c)}' $1.tmp

#rm tmp
