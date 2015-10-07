"""
Generate a magnitude distribution from the template based on the fits from the SN(oo)Py module

"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

from glob import glob 

#imports from within the package
from lsst import discover
from jwst import filters, z_sne

#simlc is the euclid filters module
from simlc import build_lc, filtcov

from dist import mod

class abs_mag:
	"""
        Absolute magnitude for a given survey strategy
        
        """
	def __init__(self):
		self.this_dir, self.this_file = os.path.split(__file__)
		
        def load_dist(self):
                """
                Load the peak magnitude distribution 
                """
                dist_path = os.path.join(self.this_dir, "templates/", "comb_ir.dat")
                
                dist = np.loadtxt(dist_path, usecols= (1, 2, 3))
                return dist

        
        def final_z_dist(self, n, band, sys,ep):
                """
                Output redshift distribution for the Euclid Survey
                
                Inputs: Number of SNe observed by JWST
                        
        
                Output: Redshift distribution
                """
                mag_dist = self.load_dist()[:,1]

                jwst_z_dist = z_sne().obs_redshift_dist(n, band, sys, ep)
                euclid_z_dist = build_lc().z_disc_euclid(band, sys, ep)

                return np.concatenate([euclid_z_dist, jwst_z_dist])
                
                
        def dist_err(self, zarr, sig_sys, sig_m=0.02):
                """
                Calculate the error in the distance modulus estimation
                sig_mu = sqrt(sig_sys**2 + (z/zmax)**2 * sig_m**2)
                
                Hence, sig_mu is a function of z and the error on the photometry. 
                The expression and the value for sig_m is taken from Cardone et al. 2012

                See also the expression in Astier et al. 2014
                
                """

                sig_mu = np.sqrt(sig_sys**2 + (zarr/max(zarr))**2 * sig_m**2)
                return sig_mu
                
        def cosmology_array(self, std=[], redfile =  'redshift_distribution_jwst_300.dat'):
                """
                Inputs: redshift distribution file
                
                Output: Array for cosmology
                
                
                """
                
                redpath = os.path.join(self.this_dir, redfile)
                
                z_array = np.loadtxt(redpath)
                
                in_gauss = self.load_dist()

                      
                if std:
                        s_val = std[0]
                        
                else:
                        s_val = np.std(in_gauss[:,1])


                peak_real = np.random.normal(np.mean(in_gauss[:,1]), s_val, len(z_array))

                st = np.std(in_gauss[:,-1])
                
                #err_real =  np.random.uniform(np.mean(in_gauss[:,-1])- 3*st, np.mean(in_gauss[:,-1])+3*st, len(z_array))
                err_real = self.dist_err(z_array, s_val)
                
                d_mod = np.array([mod(ll) for ll in z_array])


                sn = np.array(['sn'+str(k) for k in range(len(z_array)) ])

                vs = np.vstack([sn, z_array, peak_real+d_mod, err_real]).T
                

                return vs

        def cosfit_form(self,name, std=[], path = '/Users/lapguest/workspaces/euclidsims/snfile/', redfile ='redshift_distribution_jwst_300.dat'):

                """

               Given a path to the cosfitter working directory, write the simulation files

                """

                #obtain the array for the final cosmology fit 
                cos = self.cosmology_array(std=std, redfile=redfile)


                sn=cos[:,0]; zer=np.zeros(len(cos))
                #print cos
                st=np.ones(len(cos))
                vs=np.vstack([sn, cos[:,1], cos[:,1], zer, cos[:,2], cos[:,3], st, zer, zer, zer, zer, zer, zer]).T
                assert len(vs[0]) == 13
                np.savetxt(path+name, vs, fmt='%s')
