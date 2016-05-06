import numpy as np
import sys
import matplotlib.pyplot as plt

from astropy.cosmology import wCDM

class hubble:
	def __init__(self):
		self.z=np.random.uniform(0.03, 0.8, 200)

	def plot_hd(self, fname, f='yD', showplot=True, res=False, rowtoskip=2):
		cols=np.loadtxt(fname, usecols=(1,2,3), skiprows=rowtoskip)
		#default cosmology, flat, (om, w) = (0.27, -1)
		wc = wCDM(70, 0.27, 0.73)
		
		if res:
			hres = cols[:,1] - np.mean(cols[:,1])
		else:
			hres = 5*np.log10(wc.luminosity_distance(cols[:,0]).value)+25 + cols[:,1]
		plt.errorbar(cols[:,0], hres, cols[:,2], fmt=f)
		plt.show()
        #if showplot:
        #    plt.show()
