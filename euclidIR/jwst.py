"""
JWST High-z SNe survey
"""

import sncosmo
import astropy 
import sys
import numpy as np
import matplotlib.pyplot as plt

from simlc import simlc, build_lc

class filters:
   
    def __init__(self):
        self.model = sncosmo.Model(source="Hsiao")
        self.flist = ['f070w', 'f090w','f115w', 'f150w', 'f200w', 'f277w', 'f356w']
    
    def get_jwst_filt(self, indexx, pl='No'):
        filt_val=self.flist[indexx]
        band = sncosmo.get_bandpass(filt_val)
        
        
        if pl=='Yes':
            plt.plot(band.wave, band.trans)
            plt.show()
        else:
            return band

class is_vis:
    """
    To find whether an SN is observed by JWST or not

    """


    def __init__(self):
        #hsiao et al. spectral templates
        self.model = sncosmo.Model(source="Hsiao")
       
        #JWST AB system limiting magnitude

        self.limmag = 28. 
    

        self.flist = ['f090w', 'f150w', 'f200w', 'f277w']

    def out_mod(self, band, z):
        return build_lc().set_params(band, z)
    
    def out_mag(self, band, z, sys, ep, fcosm):
        """
        What magnitude will be observed by JWST in the rest-frame Y band
        """
        mod = self.out_mod(band, z)
        mag_arr = mod.bandmag(fcosm, sys, ep)
        return mag_arr

    def filt_test(self, band):


    def limit_vis(self,  band,z,sys,ep, fcosm):
        """
        Using the JWST magnitude limit, determine whether the SN would be visible for NIRCam
        """
        ma = self.out_mag(band, z, sys,ep, fcosm)
        cond = ma < self.limmag

        vis_arr = ma[cond]
        ep_vis = np.array(ep)[cond]

        if len(ep_vis) > 0:
            print "The SN would be observed by NIRCam"
        elif len(ep_vis) == 0:
            print ""

        return vis_arr, ep_vis
        


