"""
Simulating Light Curves for the Euclid SN survey in the Deep Fields

Dependencies: astropy, sncosmo

euclid discovery in the deep drilling fields

Discovery: The peak magnitude is drawn from N(u, sigma) where mu ~ -18.47 and sigma ~ 0.13 mag (very, very crude approximation of the template drawing method)

Note: When plotting z-distributions define binsize using the arange function.

"""

import sncosmo
import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from astropy.table import Table
from scipy.integrate import simps

class simlc:
    """
    Class to create bandpasses for Euclid
    """
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
            #data_path = os.path.join(this_dir, "filters/", filt+".trans")
            #wv,fl = np.loadtxt(data_path, unpack=True)
            wv, fl = np.loadtxt(self.sourcedir+'filters/'+filt+'_CSP.dat', unpack=True)
            band = sncosmo.Bandpass(wv, fl, name=filt+'_CSP')
            return band
        except:
            print "Not a CSP filter"
            return 0

    def create_LSST_bandpass(self, filt):
        """
        Use the filter set for LSST
        """
        try:
        
            this_dir, this_file = os.path.split(__file__)
            #data_path = os.path.join(this_dir, "filters/lsst/LSST_", filt+".dat")
        
            #wv,fl = np.loadtxt(data_path, unpack=True)
           
            
            wv, fl = np.loadtxt(os.path.join(this_dir, 'filters/lsst/', 'LSST_'+filt+'.dat'), unpack=True)
            
            band = sncosmo.Bandpass(wv, fl, name=filt+'_LSST')
            return band
        except:
            print "Not an LSST filter"
            return 0
            #sys.exit(1)
            
    def redshifts(self, area, tmin, tmax, zmax):
        """
        Get the redshifts of the supernovae in the survey
        Inputs: maximum redshift,
        """
        reds = list(sncosmo.zdist(0., zmax, time=(tmax-tmin), area = area))
        return reds

    
    def obs(self, taxis, band):

        """
        Setup the observation parameters into an astropy table
        """

        o=Table({'time':taxis, 'band':band, 'gain':[1.], 'skynoise':[191.27], 'zp':[24.], 'zpsys':['ab']})
        return o

    def params(self):

        """
        Model parameters for the survey
        """
        p = {'z':0.5, 't0':0}
        return p

    def reals(self, taxis, band):
        """
        Generate simulated light curves for a given parameter setting

        """

        assert len(taxis) == len(band)

        lcs  = sncosmo.realize_lcs(self.obs(taxis, band), self.model, [self.params()])
        return lcs

class build_lc:
    """
    Unlike the previous class, this one only calculates a single light curve for a given bandpass and redshift (and model)
    """

    def __init__(self):

        self.filters=['Y', 'J', 'H']

        #from Table 1 Astier et al. 2014 (the H-band is 24 mag (AB) not 24.74)
        self.limits =['24.03', '24.08', '24.00']

        #from the Deep Survey (DESIRE)
        self.deep_limits = ['25.51', '25.83', '26.08']

    def modeldef(self):
        #define the source of the template, e.g. SALT2, Hsiao et al., Nugent et al. 
        #in this case, using Hsiao for the NIR extension

        #source = sncosmo.get_source('hsiao', version='2.0')
        model=sncosmo.Model('Hsiao')
        return model

    def set_params(self, band, z, peakmag=-18.4, zpsys='ab'):
        """
        set the model parameters, most importantly the absolute magnitude scale for the tempalte spectra.
        """
        model=self.modeldef()
        model.set(z=z)
        
        try:
            model.set_source_peakabsmag(peakmag, band, zpsys)
        except:
            model.set_source_peakabsmag(peakmag, simlc().create_bandpass(band), zpsys)

        return model


    def is_discover(self, band, z, sys, ep, peakmag=-18.4,  deep='No'):

            """
            INPUTS: Filter (rest frame), Redshift, Magnitude System, Epochs of observation 
            OPTIONS: Absolute Peak magnitude

            Outputs: array of observed magnitudes that are above the detection limit

            """
            
            input_filter = filtcov(z).obs_filt(band, z)[0]
            

            fcosm = simlc().create_bandpass(input_filter[0])

            mod = self.set_params(band, z, peakmag=peakmag)

            
            mag_arr=mod.bandmag(fcosm, sys, ep)

            filt_arr = np.array(self.filters)
            #if the deep fields limit is set use, else use dedicated survey limits from Astier+2014
            if deep == 'Yes':
                limarr = np.array(self.limits)
            else:
                limarr = np.array(self.deep_limits)
            
            limmag = limarr[filt_arr == input_filter[0]]

            print limmag, mag_arr

            disc_arr = mag_arr[mag_arr < float(limmag[0])]
            
            
            disc_arr = list(disc_arr)
            if not disc_arr:
                print "No Observation above the threshold"
            
                return []
            else:
                print "SN is discovered by Euclid"
                return list(disc_arr)

    def snrate_perrett(self, z, r0 =0.21e-4, a = 1.70):
	snr = r0*pow((1+z), a)
	return snr

    def expected_z_dist(self, z=[0., 0.8], t=200, area=20):

        """
        For a 200d 20 square degree survey, the redshift distribution of expected supernovae (no magnitude cuts)

        Zmax  set by the maximum redshift of the rest frame filters

        SN rates are taken from the Perrett et al. 2012 paper for the SNe from SNLS surveyx
        """

                        
        return sorted(list(sncosmo.zdist(z[0], z[1], time=t, area=area, ratefunc=self.snrate_perrett)))

    def z_disc_euclid(self, band, sys,ep, z=[0., 0.8], t=200, area=20, deep='No', peakmag=-18.4, stdmag=0.13):

        """
        From the expected distribution, which SNe are discovered
        """
        
        #start with the expected distribution for a given zmax (given from the filter coverage of the satellite)

        expected_z= self.expected_z_dist(z=z, t=t, area=area)

        obs_z_arr=[]

        obs_mag_arr=[]
        for i, z_val in enumerate(expected_z):
            
            mag_val = peakmag#np.random.normal(peakmag, stdmag) 
            disc_arr =self.is_discover(band,z_val,sys,ep, peakmag=mag_val, deep=deep)
            
            #disc_arr = np.array(disc_arr)
            #disc_arr =list(disc_arr)

            if not disc_arr:
                print "No observations"
            else:
                obs_z_arr.append(z_val)                
                obs_mag_arr.append(mag_val)

        return np.array(obs_z_arr), np.array(obs_mag_arr)


class filtcov:

    """
    Class for filter coverage

    Determine the overlap between a rest -frame filter (redshifted) and an observer frame filter from the Euclid YJH set
    """

    def __init__(self, z):
        self.z = z
        self.y = simlc().create_bandpass('Y')
        self.j = simlc().create_bandpass('J')
        self.h = simlc().create_bandpass('H')

        self.filters = ['Y', 'J', 'H']
        

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


    def obs_filt(self, band ,z):
        
        """
        For a given instrument (in this case NIRCam), test which one has the greatest overlap with the redshifted rest-frame filter (i or Y in most cases)

        Input: rest frame filter, redshift of observation

        Output: Filter with greatest overlap, overlap value
        """

        #use the SNCosmo function for extracting the bandpass
        try:
            b = sncosmo.get_bandpass(band)
        except:
            b = simlc().create_bandpass(band)
            
        #obtain the wavelength and transmission values as python readable arrays
        wv = b.wave
        trans = b.trans

        #redshifted wavelength for the rest frame filter 
        wv_red = wv*(1+z)

        #integrate the total flux in the region of the redshifted filter
        tran_int = simps(trans, wv_red)
        
        #define array for filling the filters that have any wavelength overlap

        overlap_array = []
        print "Checking the filter list", self.filters

        for i in self.filters:
            
            #extract the bandpass for LSST
            bp = simlc().create_bandpass(i)
            
            tran_obs = bp.trans

            trans_thresh = max(tran_obs)/1e5

            wv_obs = bp.wave[bp.trans > trans_thresh]
            print wv_red[0], wv_obs[0], wv_red[-1], wv_obs[-1]
            if wv_red[0] > wv_obs[-1]:
                print "The filter being tested is", i
                print "The redshifted filter is very very red"

            elif wv_red[-1] < wv_obs[0]:
                print  "The filter being tested is", i
                print "The redshifted filter is not red enough"

            else:
                print "There is some wavelength overlap with filter", i
                overlap_array.append(i)

        print "The Euclid filters which overlap with the redshifted filter are: ", overlap_array
        
        overlap_percent=[]

        for j in overlap_array:

            bp = simlc().create_bandpass(j)
            
            trans_thresh = max(bp.trans)/1e5
            
            
            wv_obs = bp.wave[bp.trans > trans_thresh]

            cond = (wv_red > wv_obs[0] ) & (wv_red < wv_obs[-1])
            
            overlap_int=simps(trans[cond], wv_red[cond])

            overlap_percent.append([j, overlap_int*100/tran_int])

        #store the overlap percentage
        overlap_percent=np.array(overlap_percent)


        print "The percentages of the overlap are", overlap_percent

        wave_eff_arr =[]
        
        eff_wave_rf = b.wave_eff
        eff_wave_obs = eff_wave_rf *(1+z)

        for k in overlap_percent:

            if len(np.unique(overlap_percent[:,1])) < len(overlap_percent):
                
                bp = simlc().create_bandpass(k[0])
                
                wave_eff_arr.append([k[0], abs(bp.wave_eff-eff_wave_obs)])

        print "The difference between the effective wavelength for the LSST filters and the redshifted rest frame filter is:", wave_eff_arr

   

        #deal with unique and non-unique cases separately.

        if len(wave_eff_arr) > 0:
            print "In case of similar overlapping values, effective wavelengths were used to decide which filter to use"
            
            wave_eff_arr = np.array(wave_eff_arr)

    
            return wave_eff_arr[wave_eff_arr[:,1].astype('float32') == min(wave_eff_arr[:,1].astype('float32'))][0]
        else:
            print "The values for the overlap were all unique"
            return overlap_percent[overlap_percent[:,1].astype('float32')==max(overlap_percent[:,1].astype('float32')) ][0]


class redshift_distribution:
    """
    calculate the expected redshift distribution for a given restframe filter. Main constraints:
    1. filter cutoff
    2. depth
    """
    def __init__(self):
        """
        Initialise the Euclid filters
        """
        
        filtclass = simlc()
        buildclass = build_lc()
        self.euclidy = filtclass.create_bandpass('Y')
        self.euclidj = filtclass.create_bandpass('J')
        self.euclidh = filtclass.create_bandpass('H')
        self.time_period = 100
        self.z_expect = buildclass.expected_z_dist(t=self.time_period)
        self.effwave_arr = np.array([self.euclidy.wave_eff, self.euclidj.wave_eff, self.euclidh.wave_eff])
        self.filtarr = np.array([self.euclidy, self.euclidj, self.euclidh])
        
    def filtcheck(self, bandpass, z):
        """
        check if the redshifted effective wavelength is redder than the effective
        wavelength of the reddest filter (yes, its a complicated sentence)
        Input is a bandpass (as a string) and redshift
        
        """
        bp_rest = sncosmo.get_bandpass(bandpass)

        if bp_rest.wave_eff*(1+z)  > max(self.effwave_arr):
            filtfun=self.filtarr[self.effwave_arr==max(self.effwave_arr)][0]

            #condition: check what fraction of the redshifted filter is 
            cond = bp_rest.wave*(1+z) < max(filtfun.wave[filtfun.trans > 0])
            
            simp_prod = simps(bp_rest.trans, bp_rest.wave)
            simp_prod_cond = simps(bp_rest.trans[cond],bp_rest.wave[cond])
            if simp_prod_cond/simp_prod > .75:
                return 1
            else:
                print "rest-frame filter is too red for observation"
                return 0
        else:        
            return 1
    def filt_cons_redshift(self, bandpass):
        zarr = self.z_expect
        
        truth_arr=[]
        for i, zval in enumerate(zarr):
            retval=self.filtcheck(bandpass, zval)
            truth_arr.append(retval)

        truth_arr=np.array(truth_arr)
        print truth_arr
        return np.array(zarr)[truth_arr>0.]
