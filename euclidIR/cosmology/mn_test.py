"""

Likelihood function to evaluate with PyMultiNest for parameter estimation

"""

import pymultinest as pmn
import numpy as np
import pandas as pd
import sys
from scipy.integrate import quad

#load the data from a file in the snfile/ directory within euclid sims

sim_data = pd.read_csv('/Users/lapguest/workspaces/euclidsims/snfile/'+sys.argv[1], usecols=(0, 1, 4 ,5), delim_whitespace=True, header = None, names=['sn', 'z', 'mag', 'error'])

#intrinsic scatter in the sample to evaluate the chi square

sig_int = 0.12



def wa_func(z, om, w0, wa):
	"""

	Function to evaluate the distance 

	Carroll et al. 1992, C code by B. Leibundgut

	"""
	a= 1./(1+z)
	w = w0 + wa*(1-a)
	alpha = 3.*(1+w)
	res = 1./np.sqrt(pow((1.+z), 3)*om+(1-om)*pow((1+z), alpha))	
	return res

def w_func(z, om, w):
	a= 1./(1+z)
	alpha = 3.*(1+w)
	res = 1./np.sqrt(pow((1.+z), 3)*om+(1-om)*pow((1+z), alpha))	
	return res

c= 2.99e5	#check precision

def m_func(z, om, w0, wa, h0):
	"""
	Final function to evaluate the peak apparent magnitude
	
	"""
	
	out = np.empty_like(z)
        for i, z_value in enumerate(z):
                out[i] = quad(wa_func, 0, z_value, args=(om, w0, wa))[0]	
                 
        dpm = out*c/h0
	ret = dpm * (1+z)
	return 5*np.log10(ret) + 25 - 18.4

def m_func_const(z, om, w, h0):
	out = np.empty_like(z)
        for i, z_value in enumerate(z):
                out[i] = quad(w_func, 0, z_value, args=(om, w))[0]	
                 
        dpm = out*c/h0
	ret = dpm * (1+z)
	return 5*np.log10(ret) + 25 - 18.4

def llhood_w_const(model_param, ndim, nparam):
        """
	Likelihood evaluation for MultiNest to compute 
	
	Model parameters: H_0, Omega_m, w_0, w_a
	"""

	h0, om, w = [model_param[i] for i in xrange(3)]

	#model function for a given redshift array. 
        model = m_func_const(sim_data.z.values, om, w, h0)
	return -0.5 * np.sum(((sim_data.mag.values - model)**2/(sim_data.error.values**2 + sig_int**2)))
        
def llhood(model_param, ndim, nparam):
        """
	Likelihood evaluation for MultiNest to compute 
	
	Model parameters: H_0, Omega_m, w_0, w_a
	"""

	h0, om, w0, wa = [model_param[i] for i in xrange(4)]

	#model function for a given redshift array. 
        model = m_func(sim_data.z.values, om, w0, wa, h0)
	return -0.5 * np.sum(((sim_data.mag.values - model)**2/(sim_data.error.values**2 + sig_int**2)))

def prior_transform(cube, ndims, nparams):
	cube[0] = (cube[0]*50)+50
	cube[1] = cube[1]*2	
	cube[2] = cube[2]*5 - 5
	cube[3] = (cube[3]*10) - 5

def prior_transform_const(cube, ndims, nparams):
	cube[0] = (cube[0]*50)+50
	cube[1] = cube[1]*2	
	cube[2] = cube[2]*5 - 5
	#cube[3] = (cube[3]*10) - 5
	
