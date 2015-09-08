"""
Using LSST for a rest frame Near Infrared survey

Aims:
1. A High-z discovery machine for space-based follow up
2. Optical data to complement sparsely sampled light curves in the YJH

Limits: Filter for the discovery

"""
import sncosmo
import astropy
import numpy as np
import matplotlib.pyplot as plt
import os

from simlc import simlc
from jwst import z_sne, filters


class discover:
    """
    Find whether an SN is discovered by LSST to be followed-up by JWST
    """
    def __init__(self):
        self.filters=['u', 'g','r', 'i', 'z', 'y4']

        self.limits = [26.47, 26.35, 25.96, 25.50, 24.51]

        self.this_dir, self.this_file = os.path.split(__file__)

        self.filt_path = os.path.join(self.this_dir, "filters/lsst/")

         
        
    def obs_filt(self, band):
        
        """
        For a given instrument (in this case NIRCam), test which one has the greatest overlap with the redshifted rest-frame filter (i or Y in most cases)

        Input: rest frame filter, redshift of observation

        Output: Filter with greatest overlap, overlap value
        """

        #use the SNCosmo function for extracting the bandpass
        b = sncosmo.get_bandpass(band)
        
        #obtain the wavelength and transmission values as python readable arrays
        wv = b.wave
        trans = b.trans

        #redshifted wavelength for the rest frame filter 
        wv_red = wv*(1+z)

        #integrate the total flux in the region of the redshifted filter
        tran_int = simps(trans, wv_red)
        
        #define array for filling the filters that have any wavelength overlap

        overlap_array = []
        print "Checking the filter list", self.flist

        for i in self.filters:
            
            #extract the bandpass for LSST
            bp = simlc.simlc().create_LSST_bandpass(i)
            
            wv_obs= bp.wave
            tran_obs = bp.trans

            
            if wv_red[0] > wv_obs[-1]:
                print "The filter being tested is", i
                print "The redshifted filter is very very red"

            elif wv_red[-1] < wv_obs[0]:
                print  "The filter being tested is", i
                print "The redshifted filter is not red enough"

            else:
                print "There is some wavelength overlap with filter", i
                overlap_array.append(i)

        print "The LSST filters which overlap with the redshifted filter are: ", overlap_array
        
        overlap_percent=[]
        for j in overlap_array:

            bp = simlc.simlc().create_LSST_bandpass(i)
            
            wv_obs = bp.wave

            cond = (wv_red > wv_obs[0] ) & (wv_red < wv_obs[-1])
            
            overlap_int=simps(trans[cond], wv_red[cond])

            overlap_percent.append([j, overlap_int*100/tran_int])

        #store the overlap percentage
        overlap_percent=np.array(overlap_percent)


        print "The percentages of the overlap are", overlap_percent
        return overlap_percent[overlap_percent[:,1].astype('float32')==max(overlap_percent[:,1].astype('float32')) ][0]
        
        def is_discover(self, band, z, sys, ep):
            """
            For a given 
            """
            fcosm = self.obs_filt(band)[0]

            mag_arr=filters().out_mag(band,z, sys, ep, fcosm)
            
            filt_arr = np.array(self.filters)
            limmag = np.array(self.limits)[filt_arr == fcosm]
            
            disc_arr = mag_arr[mag_arr < limmag]

            if len(disc_arr) > 0:
                print "SN is discovered by LSST"
                return disc_arr
            else:
                print "No Observation above the threshold"
            
