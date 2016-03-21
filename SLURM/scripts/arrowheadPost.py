import os
import sys
import math
import numpy as np
from subprocess import check_output

resolutions = [5000, 2000, 1000, 500]

topDir = os.getcwd()
groupname = os.path.split(topDir)[1]
chroms = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX', 'chr4']

def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])
rows = wc(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_res{1}.txt_{1}_blocks'.format(groupname, resolutions[0]))
header = np.genfromtxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_res{1}.txt_{1}_blocks'.format(groupname, resolutions[0]), dtype=np.str, skip_footer = rows-1)

arrowheadIn = np.loadtxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_res{1}.txt_{1}_blocks'.format(groupname, resolutions[0]), dtype=np.str, skiprows=1)
for i in range(1, len(resolutions)):
	a = np.vstack((arrowheadIn, np.loadtxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_res{1}.txt_{1}_blocks'.format(groupname, resolutions[i]), dtype=np.str, skiprows=1)))
	arrowheadIn = a

# Sort by corner score (column 9)
ind = np.argsort(arrowheadIn[:,9])[::-1]
sorted = [arrowheadIn[i] for i in ind]
sorted = np.array(sorted)

# Get the indices of the TADs that conflict in the sorted array.
idx = []	
for i in range(0, len(sorted)):
	if i in idx:
		continue
	line = sorted[i]
	j = i+1
	while j < len(sorted):
		test = sorted[j]
		if line[0] != test[0]: # Not on same chromosome so move on
			j += 1
			continue
		if int(line[2]) - int(line[1]) >= (int(test[2]) - int(test[1])): # Line is bigger
			if int(line[1]) <= int(test[1]) <= int(line[2]):
				idx.append(j)
			if int(line[1]) <= int(test[2]) <= int(line[2]):
				idx.append(j)
		else: # Test is bigger
			if int(test[1]) <= int(line[1]) <= int(test[2]):
				idx.append(j)
			if int(test[1]) <= int(line[2]) <= int(test[2]):
				idx.append(j)
		j += 1

sortedFiltered = np.delete(sorted, np.unique(idx), axis=0)

# Sort by chromosome, TAD start, TAD end	
idx2 = np.lexsort((sortedFiltered[:,5], sortedFiltered[:,2], sortedFiltered[:,4], sortedFiltered[:,1], sortedFiltered[:,3], sortedFiltered[:,0]))
finalTADs = [sortedFiltered[i] for i in idx2]

# Remove heterochromatic and poorly assembled chromosomes.
idx3 = []
for i in range(0, len(finalTADs)):
	if finalTADs[i][0] not in chroms:
		idx3.append(i)
finalTADs2 = np.delete(finalTADs, idx3, axis=0)		

print 'Saving '+str(len(finalTADs2))+' TADs, excluding heterochromatic and poorly mapped chromosomes.'
np.savetxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_Arrowhead.txt'.format(groupname), np.vstack((header, finalTADs2)), delimiter='\t', fmt='%s')

np.savetxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_Arrowhead.bed'.format(groupname), finalTADs2[:,[0,1,2]], delimiter='\t', fmt='%s')
