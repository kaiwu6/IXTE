# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 19:31:19 2016
@author: wk
"""
#from functools import reduce
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import rc
from scipy.optimize import curve_fit
from const import K_ARRAY
#from BSEquation import *

cos = np.cos
def band_expansion(k, *t):
  '''hopping terms'''
  x, y = k
  f_test = t[0] \
    + 2.0*t[1]*(cos(x)+cos(y)+cos(x+y))\
    + 2.0*t[2]*(cos(2.0*x)+cos(2.0*y)+cos(2.0*(x+y)))\
    + 2.0*t[3]*(cos((x+2.0*y))+cos((2.0*x+y))+cos((x-y))) \
    + 2.0*t[4]*(cos(2.0*(x+2.0*y))
                + cos(2.0*(y-2.0*(x+y)))+cos(2.0*(-(x+y)+2.0*x)))
  return f_test

# def band_terms(k_points):
KLIST = np.transpose(K_ARRAY)
def band_fitting(data_source):
  '''curve fitting to get the hopping'''
  if isinstance(data_source, str):
    ka_kb_eng = np.loadtxt(data_source)
    eng_1d = ka_kb_eng[:, 2]
  else:
    eng_1d = data_source
  #eng_2d = eng_1d.reshape(KNUMB, KNUMB)
  weight = data_source + 1e-5
  popt, _ = curve_fit(band_expansion, KLIST, eng_1d,
                      p0=np.zeros(5), sigma=weight)
  return popt
  