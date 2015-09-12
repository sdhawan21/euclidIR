import pymultinest as pmn
import pandas as pd
import triangle

from mn_test import *


pmn.run(llhood_w_const, prior_transform_const, 3, verbose=True, n_live_point = 150)


data = pd.read_csv('chains/1-post_equal_weights.dat', names=['H0', 'om', 'w', 'lhood'], delim_whitespace=True, header=False)

triangle.corner(data.values[:,:-1], labels = data.columns[:-1], smooth=1)


