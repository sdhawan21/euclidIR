"""
JWST High-z SNe survey

Aim: As a 'third' configuration of a high-z SN survey, JWST provides the high-z arm to a low-z survey by Euclid

Dependencies : sncosmo, astropy, scipy, numpy 

"""

import sncosmo
import astropy 
import sys
import numpy as np
import matplotlib.pyplot as plt

from simlc import simlc, build_lc
from scipy.integrate import simps

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

class z_sne:
    """
    To find whether an SN is observed by JWST or not

    """


    def __init__(self):
        #hsiao et al. spectral templates
        self.model = sncosmo.Model(source="Hsiao")
       
        #JWST AB system limiting magnitude

        self.limmag = 28. 
        
        # set of filters on NIRCam (Note, not looking at MIRI)
        ## a few of them (e.g. f356w for rest-frame Y band) will not be required but is nonetheless tested for completeness

        self.flist = ['f070w', 'f090w','f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f444w']

    def out_mod(self, band, z):
        """
        Inline to set the parameters for a given filter and redshift combination
        """
        return build_lc().set_params(band, z)
    
    def out_mag(self, band, z, sys, ep, fcosm):

        """
        What magnitude will be observed by JWST in the rest-frame Y band
        """
        mod = self.out_mod(band, z)
        mag_arr = mod.bandmag(fcosm, sys, ep)
        return mag_arr

    def filt_test(self, band, z):

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

        for i in self.flist:
            
            #extract the bandpass for NIRcam
            bp = sncosmo.get_bandpass(i)
            
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

        print "The NIRcam filters which overlap with the redshifted filter are: ", overlap_array
        
        overlap_percent=[]
        for j in overlap_array:
            bp = sncosmo.get_bandpass(j)
            
            wv_obs = bp.wave

            cond = (wv_red > wv_obs[0] ) & (wv_red < wv_obs[-1])
            
            overlap_int=simps(trans[cond], wv_red[cond])

            overlap_percent.append([j, overlap_int*100/tran_int])

        #store the overlap percentage
        overlap_percent=np.array(overlap_percent)


        print "The percentages of the overlap are", overlap_percent
        return overlap_percent[overlap_percent[:,1].astype('float32')==max(overlap_percent[:,1].astype('float32')) ][0]
        
        

    def limit_vis(self,  band,z,sys,ep):
        """
        Using the JWST magnitude limit, determine whether the SN would be visible for NIRCam
        """
        fcosm = self.filt_test(band,z)[0]
        
        ma = self.out_mag(band, z, sys,ep, fcosm)
        cond = ma < self.limmag

        vis_arr = ma[cond]
        ep_vis = np.array(ep)[cond]

        if len(ep_vis) > 0:
            print "The SN would be observed by NIRCam"
        elif len(ep_vis) == 0:
            print "No epoch of the SN is observed by NIRCam"
        else:
            print "The length is not a strictly positive number"

        return vis_arr, ep_vis
        

    def zdist(self, n, z=[0.8, 1.2]):
        """
        get a redshift distribution for a number of SNe from JWST
        (aiming for ~ 100 or lower)
        """
        time = 200
        area= 20
        
        reds= list(sncosmo.zdist(z[0], z[1], time = time, area =area))
        
        while len(reds) > n:

            time -= 10
            
            reds= list(sncosmo.zdist(z[0], z[1], time=time, area=area))
        
        return reds

    def obs_redshift_dist(self, n, band, sys, ep):
        """

        Use the derived redshift distribution for the redshift range and then truncate based on the limiting magnitude of NIRCam 

        """


        #obtain the redshift distribution 
        expected_z_dist = self.zdist(n)
        
        obs_redshift_dist=[]

        for i in expected_z_dist:
            
            fcosm = self.filt_test(band, i)[0]
            vis_arr, ep_arr = self.limit_vis(band, i, sys, ep)
            if len(vis_arr) == 0 :
                print "List is empty, no detection from JWST"
            else:
                obs_redshift_dist.append(i)

        return np.array(obs_redshift_dist)

    def plot_z_comp(self, n, band, sys,ep):
        """
        Compare the JWST redshift distribution with the expected z-distribution from the rates
        """
        exp = self.zdist(n)

        obs = self.obs_redshift_dist(n, band, sys,ep)
        plt.hist(exp, histtype='step', label='Expected z')
        plt.hist(obs, histtype='step', label='Observed z')
        plt.legend(loc=0)
        plt.show()
        
            

