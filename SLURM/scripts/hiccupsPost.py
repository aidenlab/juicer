import os
import sys
import math
import numpy as np
from subprocess import check_output

topDir = os.getcwd()
groupname = os.path.split(topDir)[1]
loopsFile = topDir+'/juicer_output/{0}_randomizedpositions_inter_30_loops_hiccups.txt'.format(groupname)
chroms = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX', 'chr4']
loopsDir = topDir+'/juicer_output/loops'
if not os.path.exists(loopsDir):
	os.mkdir(loopsDir)

def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])
rows = wc(loopsFile)

header = np.genfromtxt(loopsFile, dtype=np.str, skip_footer = rows-1)

hiccupsIn = np.loadtxt(loopsFile, dtype=np.str, skiprows=1)

# Get indices of loops smaller than cutoff and remove heterochromatic and poorly assembled chromosomes.
idx = []
for i in range(0, len(hiccupsIn)):
	line = hiccupsIn[i]
	if int(line[5]) - int(line[1]) <= 10000 or line[0] not in chroms:
		idx.append(i)

hiccupsOut = np.delete(hiccupsIn, idx, axis=0)

print 'Saving '+str(len(hiccupsOut))+' loops.'
np.savetxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_loops_hiccups_gutter10.txt'.format(groupname), np.vstack((header, hiccupsOut)), delimiter='\t', fmt='%s')

# Create a bed file of loops.  For loop "regions" use the outermost boundaries of loop anchors.
np.savetxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_loops_hiccups_gutter10.bed'.format(groupname), hiccupsOut[:,[0,1,5]], delimiter='\t', fmt='%s')

# Get loop anchors.  For any anchor less than 5 kb, extend to 5 kb.

anchors1 = np.zeros((len(hiccupsOut), 3), dtype='S10')
for i in range(0, len(hiccupsOut)):
	line = hiccupsOut[i]
	if int(line[2]) - int(line[1]) < 5000:
		anchors1[i] = [line[0], int(round((float(line[2]) + float(line[1]))/2 - 2500)), int(round((float(line[2]) + float(line[1]))/2 + 2500))]
	else:
		anchors1[i] = line[0:3]

anchors2 = np.zeros((len(hiccupsOut), 3), dtype='S10')
for i in range(0, len(hiccupsOut)):
	line = hiccupsOut[i]
	if int(line[5]) - int(line[4]) < 5000:
		anchors2[i] = [line[0], int(round((float(line[5]) + float(line[4]))/2 - 2500)), int(round((float(line[5]) + float(line[4]))/2 + 2500))]
	else:
		anchors2[i] = line[3:6]

np.savetxt(loopsDir+'/{0}_loopAnchors_gutter10.bed'.format(groupname), np.vstack((anchors1, anchors2)), delimiter='\t', fmt='%s')
