from scipy.interpolate import interp1d

import numpy as np

rn=np.random.normal
c=2.997924562e5
def h0_sne(M, zp):
	"""
	Hubble constant from SNe observed in different passbands

	"""
	val=M+25-zp
	h0=pow(10, 0.2*val)
	return h0

def h0_mc(M, zp, n):
	"""
	Monte carlo to estimate H0 given the absolute magnitude the SNae zero point in a given filter for an Sn (with n realizations )

	"""

	arr=np.array([h0_sne(rn(M[0], M[1]), rn(zp[0], zp[1])) for k in range(n)])
	return np.mean(arr), np.std(arr)

def h0_withcosm(dl, z, om=0.27, ol=0.73):
	"""
	
	Using the complete expression for h0 with values of omega_m and omega_lambda, for low z  this approaches the nocosm value 
	om=0.27, ol=0.73 default
	"""
	q0=(om/2)-ol
	a=z*(1-q0)/(np.sqrt(1+2*q0*z)+1+q0*z )
	h0=(1+a)*c*z/dl
	return h0

def h0_nocosm(dl, z):
	return c*z/dl 
def arr_crt(fil1, fil2):
	"""	
	From two files create arrays of the parameter values for SN present in both
	"""
	a1=np.loadtxt(fil1, dtype='string')
	a2=np.loadtxt(fil2, dtype='string')

	arr1=[float(i[1]) for i in a1 if i[0] in a2[:,0]]
	arr2=[float(a2[a2[:,0]==i[0]][0][1]) for  i in a1  if i[0] in a2[:,0] ]

	return arr1, arr2

def lum_dist(z, om=0.27, ol=0.73, h0=70):
	q0=(om/2)-ol
	a=z*(1-q0)/(np.sqrt(1+2*q0*z)+1+q0*z )
	dl=(1+a)*c*z/h0
	return dl



















WM=0.27
WV=0.73
WK=1-WM-WV
H0=70
WR=0
c = 299792.458
def mod(z):

	az = 1.0/(1+1.0*z)

	age = 0
	n=1000        
	

	DTT = 0.0
	DCMR = 0.0


	for i in range(n):
  		a = az+(1-az)*(i+0.5)/n
  		adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
 		DTT = DTT + 1./adot
  		DCMR = DCMR + 1./(a*adot)

	DTT = (1.-az)*DTT/n
	DCMR = (1.-az)*DCMR/n

	DCMR_Mpc = (c/H0)*DCMR
	ratio = 1.00
	x = np.sqrt(abs(WK))*DCMR
	if x > 0.1:
		if WK > 0:
      			ratio =  0.5*(np.exp(x)-np.exp(-x))/x 
    		else:
     	
  
      			ratio = np.sin(x)/x
	else:
    		y = x*x
    	if WK < 0: y = -y
    	ratio = 1. + y/6. + y*y/120.
	DCMT = ratio*DCMR
	DA = az*DCMT	
	DL = DA/(az*az)
	DL_Mpc = (c/H0)*DL
	mu=5*np.log10(DL_Mpc)+25
	return mu



"""
All array creation functions below

"""













def rev_arr_crt(fil1, fil2, n):
	"""	
	From two files create arrays of the parameter values for SN present in both
	"""
	a1=np.loadtxt(fil1, dtype='string', delimiter='&')
	a2=np.loadtxt(fil2, dtype='string', delimiter='&')

	nm=np.array([i[0][:-1] for i in a2])

	arr1=[i[n][2:6] for i in a1 if 'SN'+i[0][:-1] in nm]
	
	arr2=[float(a2[nm=='SN'+i[0][:-1]][0][3]) for  i in a1  if 'SN'+i[0][:-1] in nm ]
	
	f1=[]; f2=[]
	for k in range(len(arr1)):
		try:
			f1.append(float(arr1[k]))
			f2.append(float(arr2[k]))
		except:
			arr1[k]	
	return np.array(f1), np.array(f2)


def self_arr_crt(fil1, fil2, n, n1):
	"""	
	From two files create arrays of the parameter values for SN present in both
	"""
	a1=np.loadtxt(fil1, dtype='string', delimiter='&')
	a2=np.loadtxt(fil2, dtype='string', delimiter='&')

	#nm=np.array([i[0][:-1] for i in a2])

	arr1=[i[n][2:6] for i in a1 if i[0] in a2[:,0]]
	
	arr2=[a2[a2[:,0]==i[0]][0][n1][2:6] for  i in a1  if i[0] in a2[:,0] ]
	
	f1=[]; f2=[]
	for k in range(len(arr1)):
		try:
			finp=float(arr1[k])
			finp2=float(arr2[k])
			f1.append(finp)
			f2.append(finp2)
		except:
			arr1[k]	
	return np.array(f1), np.array(f2)	


def dm_arr_crt(fil1, fil2, n, n1):
	"""	
	From two files create arrays of the parameter values for SN present in both
	"""
	a1=np.loadtxt(fil1, dtype='string', delimiter='&')
	a2=np.loadtxt(fil2, dtype='string', delimiter='&')
	
	
	
	dm=np.loadtxt('/home/sdhawan/workspaces/nir_ref_report/files/dist_tab.tex', dtype='string', delimiter='&')
	#print dm[:,0], a2[:,0]
	#nm=np.array([i[0][:-1] for i in a2])

	arr1=[]; arr2=[]
	for i in a1:
		if i[0] in a2[:,0]:
			try:
				if 'SN'+i[0][:-1]+'\t' in dm[:,0]:
					arr1.append([i[n][2:6], dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][5][1:-1]])
					arr2.append([a2[a2[:,0]==i[0]][0][n1][2:6], dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][5][1:-1]])
				elif 'SN'+i[0][:-2]+'\t\t' in dm[:,0]:
					arr1.append([i[n][2:6], dm[dm[:,0]=='SN'+i[0][:-2]+'\t\t'][0][5][1:-1]])
					arr2.append([a2[a2[:,0]==i[0]][0][n1][2:6], dm[dm[:,0]=='SN'+i[0][:-2]+'\t\t'][0][5][1:-1]])
			except:
				i[0]
	#arr1=[i[n][2:6] for i in a1 if i[0] in a2[:,0] ]
	#arr1=[[i[n][2:6], dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][5][1:-1]] for i in a1 if i[0] in a2[:,0] and 'SN'+i[0][:-1]+'\t' in dm[:,0]]
	
	
	#arr2=[[a2[a2[:,0]==i[0]][0][n1][2:6], dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][5][1:-1]] for  i in a1  if i[0] in a2[:,0] and 'SN'+i[0][:-1]+'\t' in dm[:,0]]
	
	f1=[]; f2=[]
	for k in range(len(arr1)):
		try:
			finp=float(arr1[k][0])-float(arr1[k][1])
			finp2=float(arr2[k][0])-float(arr2[k][1])
			f1.append(finp)
			f2.append(finp2)
		except:
			arr1[k]	
	return np.array(f1), np.array(f2)
	
		
def mix_arr_crt(fil1, fil2, n, n1):
	"""	
	From two files create arrays of the parameter values for SN present in both
	"""
	a1=np.loadtxt(fil1, dtype='string', delimiter='&')
	a2=np.loadtxt(fil2, dtype='string', delimiter='&')

	
	dm=np.loadtxt('/home/sdhawan/workspaces/nir_ref_report/files/dist_tab.tex', dtype='string', delimiter='&')
	#print dm[:,0]
	#nm=np.array([i[0][:-1] for i in a2])
	
	
	
	if n % 2 == 1:
		arr1=[]; arr2=[]
		for i in a1:
			if i[0] in a2[:,0]:
				try:
					if 'SN'+i[0][:-1]+'\t' in dm[:,0]:
						arr1.append([i[n][2:7], dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][5][1:-1]])
						
					elif 'SN'+i[0][:-2]+'\t\t' in dm[:,0]:
						arr1.append([i[n][2:7], dm[dm[:,0]=='SN'+i[0][:-2]+'\t\t'][0][5][1:-1]])
					arr2.append(a2[a2[:,0]==i[0]][0][n1][2:6])
				except:
					i[0]
		
		#arr1=[[i[n][2:6], dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][5][1:-1]] for i in a1 if i[0] in a2[:,0] and 'SN'+i[0][:-1]+'\t' in dm[:,0]]
	
		#arr2=[a2[a2[:,0]==i[0]][0][n1][2:6] for  i in a1  if i[0] in a2[:,0] and 'SN'+i[0][:-1]+'\t' in dm[:,0]]
		
	elif n % 2 == 0:
		arr1=[]; arr2=[]
		for i in a1:
			if i[0] in a2[:,0]:
				try:
					if 'SN'+i[0][:-1]+'\t' in dm[:,0]:
						arr2.append([i[n][2:7], dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][5][1:-1]])
						
					elif 'SN'+i[0][:-2]+'\t\t' in dm[:,0]:
						arr2.append([i[n][2:7], dm[dm[:,0]=='SN'+i[0][:-2]+'\t\t'][0][5][1:-1]])
					arr1.append(a2[a2[:,0]==i[0]][0][n1][2:6])
				except:
					i[0]
		#arr2=[[i[n][2:6], dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][5][1:-1]] for i in a1 if i[0] in a2[:,0] and 'SN'+i[0][:-1]+'\t' in dm[:,0]]
	
		#arr1=[a2[a2[:,0]==i[0]][0][n1][2:6] for  i in a1  if i[0] in a2[:,0] and 'SN'+i[0][:-1]+'\t' in dm[:,0]]	
	


	f1=[]; f2=[]
	for k in range(len(arr1)):
		try:
			if n % 2 == 1:
				
				finp=float(arr1[k][0])-float(arr1[k][1])
				
				finp2=float(arr2[k])
			elif n % 2 == 0:
				finp=float(arr1[k])
			finp2=float(arr2[k][0])-float(arr2[k][1])
			f1.append(finp)
			f2.append(finp2)
		except:
			arr1[k]	
	return np.array(f1), np.array(f2)
	
def dm1_arr_crt(fil1,  n):
	"""	
	From two files create arrays of the parameter values for SN present in both
	"""
	a1=np.loadtxt(fil1, dtype='string', delimiter='&')
	#a2=np.loadtxt(fil2, dtype='string', delimiter='&')
	dm=np.loadtxt('/home/sdhawan/workspaces/nir_ref_report/files/dist_tab.tex', dtype='string', delimiter='&')
	arr1=[]; arr2=[]
	for i in a1:
		
		try:
			if 'SN'+i[0][:-1]+'\t' in dm[:,0]:
				arr2.append([float(i[n][2:7]), float(dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][3][1:-1])])
					
			elif 'SN'+i[0][:-2]+'\t\t' in dm[:,0]:
				arr2.append([float(i[n][2:7]), float(dm[dm[:,0]=='SN'+i[0][:-2]+'\t\t'][0][3][1:-1])])
				#arr1.append(a2[a2[:,0]==i[0]][0][n1][2:6])
		except:
			i[0]
	
	
	#print dm[:,0]
	#nm=np.array([i[0][:-1] for i in a2])
	
	
	
	
	arr2=np.array(arr2)
	return arr2[:,0].astype('float32'), arr2[:,1].astype('float32')
	
	
def tl_arr_crt(n, par):
	a1=np.loadtxt('../files/lira_tab.tex', dtype='string', delimiter='&')
	arr2=[]
	print par
	if par[0] == 't':
		
		ff=np.loadtxt('../files/mod_tabJ-exp.tex', dtype='string', delimiter='&')
	
		for i in a1:
		
			try:
				if i[0][:-1]+' ' in ff[:,0]:
					print ff[ff[:,0]==i[0][:-1]+' '][0][n][2:6]
					arr2.append([float(i[3][1:]), float(ff[ff[:,0]==i[0][:-1]+' '][0][n][2:6])])
					
				elif i[0][:-2]+'  ' in ff[:,0]:
					print ff[ff[:,0]==i[0][:-2]+'  '][0][n][1:-1]
					arr2.append([float(i[3][1:]), float(ff[ff[:,0]==i[0][:-2]+'  '][0][n][2:6])])
					
			except:
				i[0]
	

	elif par == 'Dm15':
		dm=np.loadtxt('/home/sdhawan/workspaces/nir_ref_report/files/dist_tab.tex', dtype='string', delimiter='&')
		#arr2=[]
		for i in a1:
		
			try:
				#print i[0], dm[0][0]
				if 'SN'+i[0][:-1]+'\t' in dm[:,0]:
					
					arr2.append([float(i[3]), float(dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][3][1:-1])])
					
				elif 'SN'+i[0][:-2]+'\t\t' in dm[:,0]:
					arr2.append([float(i[3]), float(dm[dm[:,0]=='SN'+i[0][:-2]+'\t\t'][0][3][1:-1])])
				#arr1.append(a2[a2[:,0]==i[0]][0][n1][2:6])
			except:
				i[0]
	
	else:
		arr2=np.ones([10, 2])
	arr2=np.array(arr2)
	
	
	return arr2[:,0], arr2[:,1]

def inv_tl_arr_crt(n, par):
	a1=np.loadtxt('../files/lira_tab.tex', dtype='string', delimiter='&')
	arr2=[]
	print par
	if par[0] == 't':
		
		ff=np.loadtxt('../files/mod_tabJ-exp.tex', dtype='string', delimiter='&')
	
		for i in a1:
		
			try:
				if i[0][:-1]+' ' in ff[:,0]:
					 
					arr2.append([float(i[3][1:]), float(ff[ff[:,0]==i[0][:-1]+' '][0][n][2:6])])
					
				elif i[0][:-2]+'  ' in ff[:,0]:
					print ff[ff[:,0]==i[0][:-2]+'  '][0][n][1:-1]
					arr2.append([float(i[3][1:]), float(ff[ff[:,0]==i[0][:-2]+'  '][0][n][2:6])])
					
			except:
				i[0]
	

	elif par == 'Dm15':
		dm=np.loadtxt('/home/sdhawan/workspaces/nir_ref_report/files/dist_tab.tex', dtype='string', delimiter='&')
		#arr2=[]
		for i in a1:
		
			try:
				#print i[0], dm[0][0]
				if 'SN'+i[0][:-1]+'\t' in dm[:,0]:
					
					arr2.append([float(i[3]), float(dm[dm[:,0]=='SN'+i[0][:-1]+'\t'][0][3][1:-1])])
					
				elif 'SN'+i[0][:-2]+'\t\t' in dm[:,0]:
					arr2.append([float(i[3]), float(dm[dm[:,0]=='SN'+i[0][:-2]+'\t\t'][0][3][1:-1])])
				#arr1.append(a2[a2[:,0]==i[0]][0][n1][2:6])
			except:
				i[0]
	
	arr2=np.array(arr2)
	
	
	return arr2[:,1], arr2[:,0]



