import os
import sys
from subprocess import check_output
import numpy as np

# Just need the number of unique (nodups) reads:
def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])
numUnique = wc(os.getcwd()+'/juicer_aligned/merged_nodups.txt')

dups = np.loadtxt(os.getcwd()+'/juicer_aligned/dups.txt', dtype=np.str, usecols=(0,1,2,4,5,6))
#Columns are strand1, chr1, position1, strand2, chr2, position2

def testDups(array, i, j):
	# Duplicates are on same chromosome, same strand, and within 4 bp of each other.
	return array[i][1] == array[j][1] and array[i][4] == array[j][4] and array[i][0] == array[j][0] and array[i][3] == array[j][3] and abs(float(array[i][2])-float(array[j][2])) <= 4 and abs(float(array[i][5])-float(array[j][5])) <= 4
	
dupCount = np.zeros(len(dups))

# First figure out if last series of entries are duplicates of each other:
x=1
y=1
while testDups(dups, len(dups)-x-1, len(dups)-x) == True:
	testDups(dups, len(dups)-2, len(dups)-1)
	x += 1
	y += 1
dupCount[-x] = y

# Now do the rest	
i=0
j=0
while i < len(dups)-x:
	dupCount[i] += 1
	j = i + 1
	while testDups(dups, i, j) == True:
		dupCount[i] += 1
		j += 1
	i = j

# These counts are the observed read counts in the DUPLICATES list, need to add 1 to every non-zero element in the list since it already also exists in the nodups list.

dupReadCounts = dupCount[np.nonzero(dupCount)]

dupReadCounts += 1

# Merge this with an appropriately size array of 1's to reflect the true read count:

obsReadCounts = np.concatenate((np.ones(int(numUnique) - len(dupReadCounts)), dupReadCounts))

np.savetxt(os.getcwd()+'/preseq_output/obsReadCounts.txt', obsReadCounts)
