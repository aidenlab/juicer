import os
import sys
import math
import numpy as np
from subprocess import check_output
import pybedtools

resolutions = [5000, 2000, 1000]

topDir = os.getcwd()
groupname = os.path.split(topDir)[1]
chroms = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX', 'chr4']
chromSizes = {'chr2L': 23011544, 'chr2R': 21146708, 'chr3L': 24543557, 'chr3R': 27905053, 'chrX': 22422827, 'chr4': 1351857}

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
idx2 = np.lexsort((sortedFiltered[:,5].astype(int), sortedFiltered[:,2].astype(int), sortedFiltered[:,4].astype(int), sortedFiltered[:,1].astype(int), sortedFiltered[:,3], sortedFiltered[:,0]))
finalTADs = [sortedFiltered[i] for i in idx2]

# Remove heterochromatic and poorly assembled chromosomes.
idx3 = []
for i in range(0, len(finalTADs)):
	if finalTADs[i][0] not in chroms:
		idx3.append(i)
finalTADs2 = np.delete(finalTADs, idx3, axis=0)	

# Save the TADs in juicebox format.
print 'Saving '+str(len(finalTADs2))+' TADs, excluding heterochromatic and poorly mapped chromosomes.'
np.savetxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_Arrowhead.txt'.format(groupname), np.vstack((header, finalTADs2)), delimiter='\t', fmt='%s')

# Save the TADs in bed format.

names = np.zeros((len(finalTADs2), 3), dtype='|S12')
for i in range(0, len(finalTADs2)):
	names[i] = ['TAD_{0}'.format(i), '0', '.']

np.savetxt(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_Arrowhead.bed'.format(groupname), np.hstack((finalTADs2[:,[0,1,2]], names)), delimiter='\t', fmt='%s')

# Save the boundaries in bed format.
starts = np.zeros((len(finalTADs2), 3), dtype='|S12')
for i in range(0, len(starts)):
	starts[i] = ['TAD_start', '0', '+']
startBoundary = np.hstack((np.hstack((finalTADs2[:,[0,1,]], finalTADs2[:,[1,]].astype(int)+1)), starts))

stops = np.zeros((len(finalTADs2), 3), dtype='|S12')
for i in range(0, len(starts)):
	stops[i] = ['TAD_stop', '0', '-']
endBoundary = np.hstack((np.hstack((np.hstack((finalTADs2[:,[0,]], finalTADs2[:,[2,]])), finalTADs2[:,[2,]].astype(int)+1)), stops))

np.savetxt(topDir+'/juicer_output/TADs/{0}_allBoundary.bed'.format(groupname), np.vstack((startBoundary, endBoundary)), delimiter='\t', fmt='%s')

# Save the interbands in bed format.  Separate into left and right interbands.

TADsbed = pybedtools.BedTool(topDir+'/juicer_output/{0}_randomizedpositions_inter_30_TADs_Arrowhead.bed'.format(groupname))
dm3_major_bed = pybedtools.BedTool('/scratch/PI/kornberg/keagen/referencegenomes/fly/dm3/dm3.major.bed')
TADsbed.complement(g='/scratch/PI/kornberg/keagen/referencegenomes/fly/dm3/dm3.chrom.sizes').intersect(dm3_major_bed, u=True, output=topDir+'/juicer_output/{0}_randomizedpositions_inter_30_interbands_Arrowhead.bed'.format(groupname))

leftInterbands = np.zeros((len(finalTADs2), 3), dtype='|S12')
for i in range(0, len(finalTADs2)):
	if finalTADs2[i-1][0] == finalTADs2[i][0]: # Previous TAD on same chromosome
		leftInterbands[i][0] = finalTADs2[i][0]
		leftInterbands[i][1] = finalTADs2[i-1][2]
		leftInterbands[i][2] = finalTADs2[i][1]
	else: # Previous TAD on different chromosome, interband start is the start of the chromosome
		leftInterbands[i][0] = finalTADs2[i][0]
		leftInterbands[i][1] = 0
		leftInterbands[i][2] = finalTADs2[i][1]

rightInterbands = np.zeros((len(finalTADs2), 3), dtype='|S12')
for i in range(0, len(finalTADs2)):
	if i+1 == len(finalTADs2):
		rightInterbands[i][0] = finalTADs2[i][0]
		rightInterbands[i][1] = finalTADs2[i][2]
		rightInterbands[i][2] = chromSizes[finalTADs2[i][0]]
	elif finalTADs2[i+1][0] == finalTADs2[i][0]: # Next TAD on same chromosome
		rightInterbands[i][0] = finalTADs2[i][0]
		rightInterbands[i][1] = finalTADs2[i][2]
		rightInterbands[i][2] = finalTADs2[i+1][1]
	else: # Next TAD on different chromosome, interband end is the end of the chromosome
		rightInterbands[i][0] = finalTADs2[i][0]
		rightInterbands[i][1] = finalTADs2[i][2]
		rightInterbands[i][2] = chromSizes[finalTADs2[i][0]]

np.savetxt(topDir+'/juicer_output/TADs/{0}_leftInterbands.bed'.format(groupname), np.hstack((leftInterbands, names)), delimiter='\t', fmt='%s')

np.savetxt(topDir+'/juicer_output/TADs/{0}_rightInterbands.bed'.format(groupname), np.hstack((rightInterbands, names)), delimiter='\t', fmt='%s')


