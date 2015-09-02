from setuptools import setup, find_packages

from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
	long_description = f.read()

setup(
	name='euclidIR',
	
	version='1.0.0',

	description='Euclid IR project',
	long_description = long_description,

	author='Suhail Dhawan',
	author_email='sdhawan@eso.org',

	#url still to be added

	#choose your license
	license='MIT',

	#see the list of classifiers
	#how mature is this code

	classifiers=[
		'Development Status :: 3 - Alpha',

		#intended audience
		'Intended Audience :: Astronomers',
		
		'License :: OSI Approved :: MIT License',
	
	
], 


)
