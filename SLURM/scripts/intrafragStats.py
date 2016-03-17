import os
import sys
import numpy as np

topDir = os.getcwd()
outputDir = sys.argv[1]

# Load merged intrafragment_statistics file.
# Format is intraFragReads, selfCircles, danglingEnds, pairingReads

if os.path.isfile(outputDir+'/intrafragment_statistics.tmp'):
	intraFragStats = np.loadtxt(outputDir+'/intrafragment_statistics.tmp', dtype=int)

	intraFragReads = np.sum(intraFragStats[:,0])
	selfCircles = np.sum(intraFragStats[:,1])
	danglingEnds = np.sum(intraFragStats[:,2])
	pairingReads = np.sum(intraFragStats[:,3])

	outFile = open(outputDir+'/intrafragment_statistics.txt', 'w')
	outFile.write('There are '+str(intraFragReads)+' intra-fragment reads.\n')
	outFile.write('There are '+str(selfCircles)+' self-ligation (cyclization) reads.\n')
	outFile.write('There are '+str(danglingEnds)+' dangling ends reads.\n')
	outFile.write('There are '+str(pairingReads)+' homolog/sister chromatid pairing reads.\n')
	outFile.close()