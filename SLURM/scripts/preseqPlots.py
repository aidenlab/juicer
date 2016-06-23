import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import bisect

topDir = os.getcwd()
groupname = os.path.split(topDir)[1]

cCurve = np.loadtxt(os.getcwd()+'/preseq_output/{0}_c_curve.txt'.format(groupname), skiprows=1)
lcExtrap = np.loadtxt(os.getcwd()+'/preseq_output/{0}_lc_extrap.txt'.format(groupname), skiprows=1)

# Plot the cCurve.

even = np.arange(0, np.amax(cCurve[:,0]), 10)

percentDups = (1-np.amax(cCurve[:,1])/np.amax(cCurve[:,0]))*100

plt.plot(cCurve[:,0]/1e6, cCurve[:,1]/1e6, 'ro', even/1e6, even/1e6, 'k--')
plt.xlabel('Observed Total Reads (M)')
plt.ylabel('Observed Distinct Reads (M)')
plt.title('Observed Complexity Plot\n{0}% PCR Duplicates (% of aligned reads)'.format("%.2f"%percentDups))
plt.savefig(os.getcwd()+'/preseq_output/{0}_c_curve.png'.format(groupname))
plt.clf()

# Plot the lcExtrap.

a = 0.9*np.amax(lcExtrap[:,1])
b = bisect.bisect_left(lcExtrap[:,1], a)
predict = int(np.amax(lcExtrap[:,1]))
"{:,}".format(predict)

fig, ax1 = plt.subplots()
ax1.plot(lcExtrap[:,0]/1e6, lcExtrap[:,1]/1e6, 'r')
ax1.axis([0, (lcExtrap[:,0][b])/1e6, 0, 0.9*np.amax(lcExtrap[:,1])/1e6])
ax1.set_xlabel('Predicted Total Reads (M)')
ax1.set_ylabel('Predicted Distinct Reads (M)', color='r')
for tl in ax1.get_yticklabels():
    tl.set_color('r')

ax2 = ax1.twinx()
ax2.plot(lcExtrap[:,0]/1e6, lcExtrap[:,1]/lcExtrap[:,0]*100, 'b')
ax2.axis([0, (lcExtrap[:,0][b])/1e6, 0, 100])
ax2.set_ylabel('Predicted Percent Distinct Reads', color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
plt.title('Predicted Complexity Plot\nLibrary Complexity Estimate: {0}'.format("{:,}".format(predict)))
plt.savefig(os.getcwd()+'/preseq_output/{0}_lc_extrap.png'.format(groupname))
plt.clf()
