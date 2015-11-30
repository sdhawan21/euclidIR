A package to simulate supernova surveys in the rest frame Near Infrared

Currently used for cosmology with Euclid, however, it can be extended to any filter-set (JWST,WFIRST being some of the first examples)

Three main Strategies:

1. Euclid DDF-only (with optical data)
2. Euclid dedicated survey
3. Euclid DDF + JWST for high-z anchor 

Dependencies
astropy, sncosmo

Installation:
	git clone https://github.com/sdhawan21/euclidIR.git

	cd euclidIR


	python setup.py install

	Installation for cosmology likelihood computation:
		Depends on pymultinest
		install the python wrapper
		pip install pymultinest (--user if no root privileges)

		Multinest installation
		Download source from Johannes Buchner's github page
	
		cd Multinest
		
		mkdir build
		
		cd build
			
		sudo cmake .. && sudo make

		run demo files to check

Redshift Distributions (a short guide):
	1. EuclidDeep are the deep drilling fields redshift distribution, they do not inlcude a low- redshift anchor unless explicitly specific with 'lowz' in the title
	2. EuclidDed is the dedicated NIR i-band survey with Euclid
	3. EuclidDDF+JWST is a combination of low-z with Euclid and higher-z with JWST	
	4. redshift_uniform is a uniform distribution U(0, zmax) as a 'fake' test scenario
	5. JLA_z is the Betoule et al. 2014 distribution
	6. 'z_euclid_100d_lowz.dat' is a shorter survey with a low-z anchor from the ground

	A. euclid + 'survey time' + 'lowz anchor redshift'+'N_SN' in the anchor


Filter set sources:
	LSST: SNANA, R.Kessler
	Euclid: P. Astier, LPNHE Supernova
	JWST : Built-in to SNCosmo (K. Barbary)
