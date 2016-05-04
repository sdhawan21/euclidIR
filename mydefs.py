import numpy as np
import os
import sncosmo

cwd = os.getcwd()
cwd+='/euclidIR'
y_file = np.loadtxt(cwd+'/filters/Y.trans')
j_file = np.loadtxt(cwd+'/filters/J.trans')
h_file = np.loadtxt(cwd+'/filters/H.trans')

y_lsst_file = np.loadtxt(cwd+'/filters/lsst/LSST_y4.dat')
u_lsst_file = np.loadtxt(cwd+'/filters/lsst/LSST_u.dat')
g_lsst_file = np.loadtxt(cwd+'/filters/lsst/LSST_g.dat')


band = sncosmo.Bandpass(y_file[:,0], y_file[:,1], name='euclidY')
sncosmo.registry.register(band)

band1 = sncosmo.Bandpass(j_file[:,0], j_file[:,1], name='euclidJ')
sncosmo.registry.register(band1)

band2 = sncosmo.Bandpass(h_file[:,0], h_file[:,1], name='euclidH')
sncosmo.registry.register(band2)


lband = sncosmo.Bandpass(y_lsst_file[:,0], y_lsst_file[:,1], name='lssty4')
sncosmo.registry.register(lband)

lband1 = sncosmo.Bandpass(u_lsst_file[:,0], u_lsst_file[:,1], name='lsstu')
sncosmo.registry.register(lband1)

lband2 = sncosmo.Bandpass(g_lsst_file[:,0], g_lsst_file[:,1], name='lsstg')
sncosmo.registry.register(lband2)




