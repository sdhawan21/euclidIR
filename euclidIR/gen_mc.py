"""
Generate a magnitude distribution from the template based on the fits from the SN(oo)Py module

"""
import numpy as np
import matplotlib.pyplot as plt
import sys

from glob import glob 

#imports from within the package
from lsst import discover
from jwst import filters
#simlc is the euclid filters module
from simlc import build_lc


class abs_mag:
	"""
        Absolute magnitude for a given survey strategy
        
        """
	def __init__(self):
		self.this_dir, self.this_file = os.path.split(__file__)
		
        def load_dist(self):
                dist_path = os.path.join(self.this_dir, "templates/", "comb_ir.dat")
                
                dist = np.loadtxt(dist_path, usecols= (1, 2, 3))
                return dist

        
        def final_dist(self, n, band, sys,ep):

                mag_dist = self.load_dist()[:,1]

                jwst_z_dist = z_sne().obs_redshift_dist(n, band, sys, ep)
                euclid_z_dist = 
                
                
                

