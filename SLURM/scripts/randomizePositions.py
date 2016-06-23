import os
import sys

print 'From "srun python", the python being used is located at' + sys.executable
print 'From "srun python", the PYTHONPATH being used is located at' + os.environ['PYTHONPATH']

import numpy as np
from subprocess import check_output


topDir = os.getcwd()
rank = int(sys.argv[1]) # Which split file to work on
sitesDir = sys.argv[2]
splitsDir = sys.argv[3]

# Load restrictions sites.

def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])
rows = wc(sitesDir)

# Create a dictionary of restriction sites where the keys are the chromosome names.
restrictionSites={}
for i in range(0, rows):
	line = np.genfromtxt(sitesDir, dtype=np.str, skip_header=i, skip_footer=rows - i - 1)
	restrictionSites[line[0]] = line[1:]


# Load merged_nodups, but do this with a reader to iterate over the rows and save memory.
#Columns are strand1, chr1, position1, fragment1, strand2, chr2, position2, fragment2, mapq1, misc1, sequence1, mapq2, sequence2, misc2

merged_nodups = np.loadtxt(splitsDir+'/merged_nodups_{0}.txt'.format("%03d"%rank), dtype=np.str)

# Now get intrafragment read statistics and randomize mapping positions.
# Due to bug in juicer reads mapping beyond the end of the chromosomes are in the last fragment+1 (i.e. there is one more restriction fragment than there should be.

intraFragReads = 0
selfCircles = 0 # Reads point away each other
danglingEnds = 0 # Reads point towards from each other
pairingReads = 0 # Reads in same directions

for i in range(0, len(merged_nodups)):
	line = merged_nodups[i]
	
	# First get intrafragment stats
	if line[1] == line[5] and line[3] == line[7]:
		intraFragReads += 1
	if line[1] == line[5] and line[3] == line[7] and line[0] == '16' and line[4] == '0':
		selfCircles += 1
	if line[1] == line[5] and line[3] == line[7] and line[0] == '0' and line[4] == '16':
		danglingEnds += 1
	if line[1] == line[5] and line[3] == line[7] and line[0] == '0' and line[4] == '0':
		pairingReads += 1
	if line[1] == line[5] and line[3] == line[7] and line[0] == '16' and line[4] == '16':
		pairingReads += 1

	# Next, randomize read positions
	
	chr1 = line[1]
	frag1Idx = int(line[3])
	if frag1Idx == 0:
		site1high = int( restrictionSites[chr1][frag1Idx] )
		randomPos1 = np.random.randint(1, site1high)
	elif frag1Idx >= len(restrictionSites[chr1]): # Accounts for bug in juicer with mapping beyond the end of the chromosome.
		site1high = int( restrictionSites[chr1][frag1Idx - 1] )
		site1low = int( restrictionSites[chr1][frag1Idx - 2] )
		randomPos1 = np.random.randint(site1low, site1high)
	else:
		site1high = int( restrictionSites[chr1][frag1Idx] )
		site1low = int( restrictionSites[chr1][frag1Idx - 1] )
		randomPos1 = np.random.randint(site1low, site1high)

	chr2 = line[5]
	frag2Idx = int(line[7])
	if frag2Idx == 0:
		site2high = int( restrictionSites[chr2][frag2Idx] )
		randomPos2 = np.random.randint(1, site2high)
	elif frag2Idx >= len(restrictionSites[chr2]): # Accounts for bug in juicer with mapping beyond the end of the chromosome.
		site2high = int( restrictionSites[chr2][frag2Idx - 1] )
		site2low = int( restrictionSites[chr2][frag2Idx - 2] )
		randomPos2 = np.random.randint(site2low, site2high)
	else:
		site2high = int( restrictionSites[chr2][frag2Idx] )
		site2low = int( restrictionSites[chr2][frag2Idx - 1] )
		randomPos2 = np.random.randint(site2low, site2high)
	
	np.put(merged_nodups[i], [2, 6], [str(randomPos1), str(randomPos2)])

# Sort the new merged_nodups array.  Keep in mind for np.lexsort the keys to sort by are in reverse order.

ind = np.lexsort((merged_nodups[:,2].astype(int), merged_nodups[:,4].astype(int), merged_nodups[:,0].astype(int), merged_nodups[:,7].astype(int), merged_nodups[:,3].astype(int), merged_nodups[:,5], merged_nodups[:,1]))

merged_nodups_randomize = [merged_nodups[i] for i in ind]

intraFragStats = np.array([intraFragReads, selfCircles, danglingEnds, pairingReads])

np.savetxt(splitsDir+'/intrafragment_statistics_{0}.txt'.format("%03d"%rank), intraFragStats.reshape(1, intraFragStats.shape[0]), fmt='%.1i')

np.savetxt(splitsDir+'/merged_nodups_randomized_{0}.txt'.format("%03d"%rank), merged_nodups_randomize, fmt='%.85s')


