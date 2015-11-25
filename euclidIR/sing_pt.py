"""
Single point Euclid Simulations (the namespaces are not well-defined, this is old code)

Aim: Single point observations for N_SN (set by the redshift distribution) supernovae from the Euclid satellite to constrain the cosmology

the errors on the photometry are taken as a normal distribution centered at the mean and sigma of the CSP-I SNe



"""

from snpy import *
from matplotlib.pyplot import *
from scipy.interpolate import interp1d
from random import *
from numpy import *
from glob import glob
from sys import argv
from time import time
from calc_kcorr import calc_kcor

start=time()
#define the speed of light 
c=2.997924562e5


def lum_dist(z, ol, om, h0):		#expression for luminosity distance in lambda CDM

	q0=(om/2)-ol
	a=z*(1-q0)/(sqrt(1+2*q0*z)+1+q0*z )
	dl=(1+a)*c*z/h0
	return dl		
	
#yh=loadtxt('dat_files/YH_col_mean_curve.dat', usecols=(0, 1))
#jh=loadtxt('dat_files/JH_col_mean_curve.dat', usecols=(0, 1))	

pt2='/share/splinter/ug_msci2/mean_dat/paper_red/'
pt = '/share/splinter/ug_msci2/'
temp=get_sn(pt+'mcmc_files/J_mean_curve.dat')
t1=get_sn(pt+'mcmc_files/H_mean_curve.dat')
t2=get_sn(pt+'mcmc_files/Y_mean_curve.dat')


#lc=temp.get_mag_table('J')
#lc1=t1.get_mag_table('H')
#lc2=t2.get_mag_table('Y')
t2 = t2.get_mag_table('Y')

#j=loadtxt('mean_dat/J_mean_curve_5.dat')
#h=loadtxt('mean_dat/H_mean_curve_5.dat')
#y=loadtxt('mean_dat/Y_mean_curve_st_ran.dat')	#y-band mean curve is changed to see the different effects fo the dm15 range
#cad=int(argv[3])

zfile = argv[3]
z=loadtxt(pt+'mcmc_files/'+zfile)
lim=float(argv[5])
w=arange(0.01, 0.8, 0.01)
dist=[round(i, 2) for i in w]
direc=argv[6]
#print min(j[:,2])
print z
j = np.loadtxt(pt2+'Y_mean_curve_nodered.dat')
#z=z[z>float(argv[7])]

#j = np.vstack([t2['MJD'], t2['Y'], t2['e_Y']]).T

#draw light curves for each redshift
for m in range(len(z)):
	ind_arr=[]
	#d=list(z).count(z[m])
	#mu=5*log10(lum_dist(z[m], 0.73, 0.27, 70))+25
	#cad1=round(cad/float(1+z[m]))
	#mjd=lc1['MJD']
	mjd = j[:,0]
	while len(set(ind_arr))<int(argv[2]):
		i=uniform(0, 1)
		ind=int(len(mjd)*i)
		if mjd[ind]>=int(argv[1]) and mjd[ind]<=int(argv[4]): 	
			indexx = ind
			ind_arr.append(indexx)
			#while ind<len(mjd):
			#	ind_arr.append(ind)
			#	ind+=cad1				
	fout=open('/share/splinter/ug_msci2/sn_files/'+direc+'randlc'+str(m+1)+'.dat', 'w')
	"""
	ind_arr=sorted(set(ind_arr))
	col=jh[ind_arr]
	if z[m]<=0.5:					#a rather silly approach to k-corrections. The corrections give out the right values upto 0.5, however, cant be trusted for z beyond. This will be modified with high redshift kcorr from spectra and possibly bruno's code
		kc=calc_kcor('J', z[m], 'JH', col[:,1])
	else:
		kc=calc_kcor('J', 0.5, 'JH', col[:,1])
	
	"""
	ph=j[:,0][indexx]#=lc['MJD'][ind_arr]	
	mag=j[:,1][indexx]#=lc['J'][ind_arr]
	mag2=j[:,2][indexx]#=lc['e_J'][ind_arr]
	
	err = j[:,3][indexx]
	err2 = j[:,4][indexx]
	
	gau=random.normal(mag, mag2)
	e_gau = random.uniform(err-err2, err+err2)
	#gau1=gau[gau-35+mu<lim]
	
	
	fout.write('SNmean'+'\t'+str(0.00008)+'\t'+str(3.02)+'\t'+str(-0.9)+'\n')
	
		
	fout.write('filter Y'+'\n')
	fout.write(str(ph)+'\t'+str(gau)+'\t'+str(e_gau)+'\n')
	
	"""
	#for rn in range(len(gau1)):
	#		fout.write(str(ph[rn])+'\t'+str(gau1[rn])+'\t'+str(mag2[rn])+'\n')		
	col=jh[ind_arr]
	if z[m]<=0.5:
		kc=calc_kcor('H', z[m], 'JH', col[:,1])
	else:
		kc=calc_kcor('H', 0.5, 'JH', col[:,1])
	ph1=h[:,0][ind_arr]#=lc1['MJD'][ind_arr]	
	mag1=h[:,1][ind_arr]#=lc1['H'][ind_arr]
	mag2=h[:,2][ind_arr]#=lc1['e_H'][ind_arr]
	gau=random.normal(mag1, mag2)
	gau1=gau[gau-35+mu+kc<lim]
	if len(gau1)>1:
		fout.write('filter H'+'\n')
		for rn in range(len(gau1)):
			fout.write(str(ph[rn])+'\t'+str(gau1[rn])+'\t'+str(mag2[rn])+'\n')
	col=yh[ind_arr]
	if z[m]<=0.5:
		kc=calc_kcor('Y', z[m], 'YH', col[:,1])
	else:
		kc=calc_kcor('Y', 0.5, 'YH', col[:,1])
	ph=y[:,0][ind_arr]	
	mag=y[:,1][ind_arr]
	mag2=y[:,2][ind_arr]
	gau=random.normal(mag, mag2)
	gau1=gau[gau-35+mu+kc<lim]
	if len(gau1)>1:
		fout.write('filter Y'+'\n')
		for rn in range(len(gau1)):
			fout.write(str(ph[rn])+'\t'+str(gau1[rn])+'\t'+str(mag2[rn])+'\n')
	"""		
	fout.close()
end=time()
print 'it took'+str(end-start)+'seconds' 


