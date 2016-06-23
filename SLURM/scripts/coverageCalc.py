import os
import sys
import math
import numpy as np
from scipy.sparse import csr_matrix

resolutions=[5000, 2000, 1000, 500]

topDir = os.getcwd()
dumpdir = sys.argv[1]
chromsize = sys.argv[2]
outputdir = sys.argv[3]

chromSizesIn = np.loadtxt(chromsize, dtype=np.str)
chromSizes = {}
for i in range(0, len(chromSizesIn)):
	chromSizes[chromSizesIn[i][0]] = int(chromSizesIn[i][1])

chroms = chromSizesIn[:,0]
chroms = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX']

outFile = open(outputdir+'/coverage.txt', 'w')
outFile.close()

for res in resolutions:
	for chr in chroms:
		outFile = open(outputdir+'/coverage.txt', 'a')
		statinfo = os.stat(dumpdir+'/{0}_observed_randomizedpositions_{1}.txt'.format(chr, res))
		if statinfo.st_size == 0:
			continue
		data = np.loadtxt(dumpdir+'/{0}_observed_randomizedpositions_{1}.txt'.format(chr, res))
		row = data[:,0]/res
		col = data[:,1]/res
		val = data[:,2]
		length = int(math.ceil(chromSizes[chr]/float(res)))

		mtx = csr_matrix((val, (row, col)), shape=(length, length))

		rowsums = np.zeros(length)
		for i in range(0, length):
			rowsums[i] = np.sum(mtx[i].todense())
		
		outFile.write('At '+str(res)+' bp resolution, for chromosomes '+chr+' the mean, KR normalized coverage is '+str(np.nanmean(rowsums))+' contacts.\n')
		for i in range(10, 100, 10):
			outFile.write('\tAt '+str(res)+' bp resolution, for chromosomes '+chr+' the KR normalized coverage in the {0}th percentile is '.format(i)+str(np.nanpercentile(rowsums, i))+' contacts.\n')

outFile.close()
