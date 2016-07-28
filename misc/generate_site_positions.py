#! /usr/bin/env python
# Generate site positions in genome from given restriction enzyme
# Juicer 1.5
import sys
import subprocess
import re

if len(sys.argv) != 3 and len(sys.argv) != 4:
	print 'Usage: %s <restriction enzyme> <genome> [location]' % (sys.argv[0])
	sys.exit(1)

if sys.argv[1]=='HindIII':
	teststring='AAGCTT'
elif sys.argv[1]=='DpnII':
	teststring='GATC'
elif sys.argv[1]=='MboI':
	teststring='GATC'
else:
	print 'Usage: %s <restriction enzyme> <genome> [location]' % (sys.argv[0])
	print '<restriction enzyme> must be defined in script, such as "HindIII" or "DpnII" '
	sys.exit(1)

genome=sys.argv[2]
if genome == 'hg19':
  filename='/seq/references/Homo_sapiens_assembly19.fasta'
elif genome == 'mm9':
  filename='/seq/references/Mus_musculus_assembly9.fasta'
elif genome == 'mm10':
  filename='/seq/references/Mus_musculus_assembly10.fasta'
elif genome == 'hg18':
  filename='/seq/references/Homo_sapiens_assembly18.fasta'
elif len(sys.argv) == 4:
  filename=sys.argv[3]
else:
  print 'Usage: %s <restriction enzyme> <genome> [location]' % (sys.argv[0])
  print '<genome> not found and [location] not defined'
  sys.exit(1)

check=[]
teststring=teststring.lower()
for i in teststring:
  check.append(i)
f=open(filename,'r')
g=open(genome+"_"+sys.argv[1]+'.txt','w')
chromosomecheck=f.readline()
while (chromosomecheck):
  if (chromosomecheck.startswith('>')):
    firststr=re.split('\s*',chromosomecheck[1:]);
    g.write(firststr[0])
    g.write(' ')
    next=f.read(len(check))
    next=next.lower()
    list=[]
    for i in next:
      list.append(i)
    counter=1
    while len(next)!=0:
      if list==check:
        counterstring=str(counter)
        g.write(counterstring)
        g.write(' ')
      list.pop(0)
      next=f.read(1)
      next=next.lower()
      if next=='\n':
        next=f.read(1)
        next=next.lower()
      if next=='>':
        break
      list.append(next)
      counter+=1
    counter+=len(check)-1
    counterstring=str(counter)
    g.write(counterstring)
    g.write('\n')
    chromosomecheck=next+f.readline()
  else:
    chromosomecheck=f.readline()
print chromosomecheck
            
f.close()
g.close()
            
