"""
Simulating Light Curves for the Euclid SN survey in the Deep Fields

Dependencies: astropy, sncosmo

To do: get the right observer frame filter
"""
import sncosmo
import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from astropy.table import Table
from scipy.integrate import simps

class simlc:
    def __init__(self):
        #define the model
        self.model = sncosmo.Model(source="Hsiao")
        
        #define the directory where all the information is stored (e.g. filters, templates, magnitudes)
        #self.sourcedir = sourcedir

    #use the bandpass function to create a readable filter transmission function
    def create_bandpass(self, filt):
        """
        Inputs: Filter (e.g. Y, J), NOTE: has to be written in the ../filters directory as a .trans file
        Outputs: SNCOSMO readable filter
        """

        try:
            this_dir, this_file = os.path.split(__file__)
            data_path = os.path.join(this_dir, "filters/", filt+".trans")
            wv,fl = np.loadtxt(data_path, unpack=True)
            band = sncosmo.Bandpass(wv, fl, name=filt+'_euclid')
            return band

        except:
            print "The filter not in the Euclid bandpasses list"
            raise

    def create_CSP_bandpass(self, filt):
        """
        Use the filter set for CSP-I from the observatories webpage (mainly i and Y) and convert to sncosmo readable format
        """
        try:
            this_dir, this_file = os.path.split(__file__)
            data_path = os.path.join(this_dir, "filters/", filt+".trans")
            wv,fl = np.loadtxt(data_path, unpack=True)
            wv, fl = np.loadtxt(self.sourcedir+'filters/'+filt+'_CSP.dat', unpack=True)
            band = sncosmo.Bandpass(wv, fl, name=filt+'_CSP')
            return band
        except:
            print "Not a CSP filter"
            return 0

    def redshifts(self, area, tmin, tmax, zmax):
        """
        Get the redshifts of the supernovae in the survey
        Inputs: maximum redshift,
        """
        reds = list(sncosmo.zdist(0., zmax, time=(tmax-tmin), area = area))
        return reds

    
    def obs(self, taxis, band):
        o=Table({'time':taxis, 'band':band, 'gain':[1.], 'skynoise':[191.27], 'zp':[24.], 'zpsys':['ab']})
        return o

    def params(self):
        """
        Model parameters for the survey
        """
        p = {'z':0.5, 't0':0}
        return p

    def reals(self, taxis, band):
        lcs  = sncosmo.realize_lcs(self.obs(taxis, band), self.model, [self.params()])
        return lcs

class build_lc:
    """
    Unlike the previous class, this one only calculates a single light curve for a given bandpass and redshift (and model)
    """

    def modeldef(self):
        #source = sncosmo.get_source('hsiao', version='2.0')
        model=sncosmo.Model('Hsiao')
        return model

    def set_params(self, band, z, zpsys='ab'):
        """
        set the model parameters, most importantly the absolute magnitude scale for the tempalte spectra.
        """
        model=self.modeldef()
        model.set(z=z)
        model.set_source_peakabsmag(-18.4, band, zpsys)
        return model


class filtcov:
    """
    Class for filter coverage

    Determine the overlap between a rest -frame filter (redshifted) and an observer frame filter from the Euclid YJH set)

    """
    def __init__(self, z):
        self.z = z
        self.y = simlc().create_bandpass('Y')
        self.j = simlc().create_bandpass('J')
        self.h = simlc().create_bandpass('H')

    
        
    def frac(self, filt1):
        """
        Fractional coverage of the filter in observer frame with the filters on board Euclid
        """
        f1 = simlc().create_bandpass(filt1)
        reds_f1 = np.vstack([f1.wave, f1.trans]).T
        
        reds_f1[:,0]*=(1+self.z)
        
        totfl = simps(reds_f1[:,1], reds_f1[:,0])

        if reds_f1[0][0] > self.y.wave[-1]:
            print 'the redshifted filter is redder than Y'
        
        else:
            cond = (reds_f1[:,0] > self.y.wave[0]) & (reds_f1[:,0] < self.y.wave[-1]) 
            t1 = reds_f1[:,1][cond]
            w1 = reds_f1[:,0][cond]
            s1 = simps(t1, w1)
            if s1 / totfl < .75:
                print "Not sufficient overlap"

            else:
                ofilt = 'Y'

      
        if reds_f1[0][0] > self.j.wave[-1]:
            print 'the redshifted filter is redder than J'
        
        else:
            cond = (reds_f1[:,0] > self.j.wave[0]) & (reds_f1[:,0] < self.j.wave[-1]) 
            t1 = reds_f1[:,1][cond]
            w1 = reds_f1[:,0][cond]
            s1 = simps(t1, w1)
            if s1 / totfl < .75:
                print "Not sufficient overlap"

            else:
                ofilt = 'J'
        
        if reds_f1[0][0] > self.h.wave[-1]:
            print 'the redshifted filter is redder than H'
        
        else:
            cond = (reds_f1[:,0] > self.h.wave[0]) & (reds_f1[:,0] < self.h.wave[-1]) 
            t1 = reds_f1[:,1][cond]
            w1 = reds_f1[:,0][cond]
            s1 = simps(t1, w1)

            if s1 / totfl < .75:
                print "Not sufficient overlap"

            else:
                ofilt = 'H'
        return ofilt

