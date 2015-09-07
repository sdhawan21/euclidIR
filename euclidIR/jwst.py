"""
JWST High-z SNe survey
"""

import sncosmo
import astropy 
import sys
import numpy as np
import matplotlib.pyplot as plt

from simlc import simlc

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

class get_sn:
    
    def __init__(self):
        


