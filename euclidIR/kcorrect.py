"""
A routine to calculate the dispersions in K-correction; using SN(oo)Py heavily

Has several potentially useful functions for convolving the spectrum with a filter transmission curve



"""

import numpy as np
import scipy.interpolate
import os
import sys
try:
    import pyfits
    have_FITS=0  #variable naming is kept parallel to the SNPy convention
except ImportError:
    sys.stderr.write('Error:  You need pyfits to run snpy.  You can get it\n')
    sys.stderr.write('        from:  http://www.stsci.edu/resources/'+\
                       'software_hardware/pyfits/\n')
    raise ImportError
#import the filters directory from SNooPY
#import filters
#from mangle_spectrum import mangle_spectrum2   #currently not exploiting the mangle function

debug=0
h = 6.626068e-27
c = 2.997925e18
ch = c * h

# This converts dm15 into a stretch that can be used to
#    stretch the Hsiao template.  This is a first correction
#    to "warp" the SED template
def dm152s(dm15):
   return 2.13 - 2.44*dm15 + 2.07*dm15**2 - 0.7*dm15**3

#define the base directory and the spectrum storage directory
base = os.path.dirname(globals()['__file__'])
spec_base = os.path.join(base,'typeIa')

#load Eric Hsiao's optical+NIR "uberspectrum" (term lifted from Chris Burns' code)
f = pyfits.open(os.path.join(spec_base, "Hsiao_SED_V3.fits"))
h3_sed = f[0].data
head = f[0].header
h3_wav = head['CRVAL1']+(np.arange(head['NAXIS1'],dtype=np.float32) - \
            head['CRPIX1'] + 1)*head['CDELT1']
f.close()

subsample = 1


def read_input(fname="input.kcorr"):
    ls = []
    empt_dict = {}
    fout=open(fname, 'r')
    for row in fout:
        ls.append(row.split())
    for i in ls:
        empt_dict[i[0]]=i[1]
    return empt_dict
        

#draw the template SED
def get_template_SED(day):
    if day < -19 or day > 70:
        return (None, None)
    else:
        return (h3_wav, h3_sed[day+20,:])

def response(wave, spec, filt, z=0, zeropad=0, photons=1, interp_method="spline", integ_method="simpsons"):
    """
    Calculate the response for spectrum through a filter function
    --> wave: 1-D array of wavelengths (angstroms)
    --> spec: 1-D array of flux values
    --> filt: 2-D array of filter wavelength and transmission
    """
    wvfilt = filt[:,0]
    transfilt = filt[:,1]
    
    if z > 0:
        swave = wave*(1.+z)
    elif z < 0:
        swave = wave/(1.+z)
    else:
        swave = wave
    if (min(wvfilt) < swave[0] or max(wvfilt) > swave[-1]) and not zeropad:
           return(-1.0)

    # Now figure out the limits of the integration:
    x_min = np.minimum.reduce(wave)
    x_max = np.maximum.reduce(wave)

    try:
        i_min = np.nonzero(num.greater(swave - x_min, 0))[0]
    except:
        i_min = 0
    try:
        i_max = np.nonzero(num.greater(swave - x_max, 0))[0]
    except:
        i_max = len(swave)-1
   
    if i_min >= 5:
        i_min -= 5
    else:
        i_min = 0
    if i_max <= len(swave)-6:
        i_max += 5
    else:
        i_max = len(swave) - 1

    trim_spec = spec[i_min:i_max+1:subsample]
    trim_wave = swave[i_min:i_max+1:subsample]
    # Now, we need to resample the response wavelengths to the spectrum:
    if interp_method == "spline":
        tck = scipy.interpolate.splrep(wvfilt, transfilt, k=1, s=0)
        fresp_int = scipy.interpolate.splev(trim_wave, tck)
    else:
        mint = scipy.interpolate.interp1d(wvfilt, transfilt, kind=interp_method)
        fresp_int = mint(trim_wave)
    # Zero out any places beyond the definition of the filter:
    fresp_int = np.where(np.less(trim_wave, x_min), 0, fresp_int)
    fresp_int = np.where(np.greater(trim_wave, x_max), 0, fresp_int)

    integrand = fresp_int*trim_spec
    if photons:
        integrand = integrand*trim_wave/ch

    if integ_method=='simpsons':
         result = scipy.integrate.simps(integrand, x=trim_wave, even='avg')
    elif integ_method=='trapz':
         result = scipy.integrate.trapz(integrand, x=trim_wave)
    else:
        result = (trim_wave[-1] - trim_wave[0])/(len(trim_wave)-1)*sum(integrand)

    return(result)
    
def K(wv, spec, f1, f2, z, zpt1, zpt2, photons=1):
    """
    define the expression for the K-correction term
    """
    # The zero-points
    #zpt1 = f1.zp
    #zpt2 = f2.zp

    # compute the response through each filter
    ## uses the response function defined above instead of the more sophisticated formalism of SN(oo)Py by C. Burns
    
    f1flux_0 = response(wv, spec, f1, photons=photons)
    f2flux_z = response(wv, spec, f2, z=z, photons=photons)
    if f1flux_0 < 0 or f2flux_z <= 0:
    # Something clearly went wrong
        return (0.0, 0)
    else:
    # Finally calculate the cross-band K Correction
    # Since we blueshift the spectrum (instead of redshift the filter)
    # the sign of the 2.5 is now positive
      kf1f2 =  2.5*np.log10(1+z) + 2.5*np.log10(f1flux_0/f2flux_z) - \
               zpt1 + zpt2
      return (kf1f2, 1)

    
        
def main():
    inputval=read_input()
    filtdir = inputval['filtdir']
    obs_dir = inputval['obsdir']
    filtname = inputval['filtinp']
    zval = float(inputval['z'])
    phasein = list(inputval['phase'])
    filtname2 = inputval['filtinp2']

    
    zpfile = np.loadtxt(filtdir+obs_dir+'filters.dat', dtype='string')
    filt1 = np.loadtxt(filtdir+obs_dir+filtname+'.dat')
    filt2 = np.loadtxt(filtdir+obs_dir+filtname2+'.dat')

    
    if filtname+'.dat' in zpfile[:,1]:
        zpf1 = float(zpfile[zpfile[:,1] == filtname+'.dat'][0][2])
    else:
        print "No value for zeropoint given, terminating script"
        sys.exit()
    kcorr_arr = []
    print filt1

    phase = int(sys.argv[1])
    #for day in phasein:
    wv, fl = get_template_SED(phase)
    kval = K(wv, fl, filt1, filt2, zval, zpf1, zpf1)
    print fl
 

    #kcorr_arr = np.array(kcorr_arr)
    #compare with the data
    datfile = sys.argv[2]
    spec_source = np.loadtxt(datfile); wv_source = spec_source[:,0]; fl_source = spec_source[:,1]
    #print fl_source
    #sys.exit()
    ksrc = K(wv_source, fl_source, filt1, filt1, zval, zpf1, zpf1)
    print kval, ksrc

    
if __name__ == "__main__":
    main()
