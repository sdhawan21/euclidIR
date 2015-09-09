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

Filter set sources:
	LSST: SNANA, R.Kessler
	Euclid: P. Astier, LPNHE Supernova
	JWST : Built-in to SNCosmo (K. Barbary)
