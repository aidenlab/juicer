import os
import sys
import math
import numpy as np
from subprocess import check_output

topDir = os.getcwd()
groupname = os.path.split(topDir)[1]
loopsFile = sys.argv[1]
chroms = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX', 'chr4']

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
np.savetxt(loopsFile, np.vstack((header, hiccupsOut)), delimiter='\t', fmt='%s')

#np.savetxt(topDir+'/juicer_output/{0}_inter_30_loops_guttered.bed'.format(groupname), hiccupsOut[:,[0,1,5]], delimiter='\t', fmt='%s')
