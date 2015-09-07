import numpy as np
import matplotlib.pyplot as plt
import sncosmo

from simlc import simlc 
from astropy.table import Table

class lcfit:
	
	def __init__(self):
		self.model = sncosmo.Model(source='Hsiao')
	
        def obs_tab(self, taxis, bandarr):
                """
                Create the astropy table with the input variables 

                
                """
                obs = Table({
                'time':taxis,
                'band':bandarr,
                'gain':[1., 1., 1.],
                'skynoise':[191, 147, 160],
                'zp': [30., 30., 30.],
                'zpsys':['ab', 'ab', 'ab']
                })

                return obs        

	def fit_lc(self, taxis, band):	
		
		data = simlc().reals(taxis, band)[0]
		mod = self.model
		res, fitted_model = sncosmo.fit_lc(data, mod, ['t0'])
		
		return res, fitted_model

		
		


