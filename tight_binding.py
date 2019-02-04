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
  k_x, k_y = k
  f_test = t[0] \
    + 2.0*t[1]*(cos(k_x)+cos(k_y)+cos(k_x+k_y))\
    + 2.0*t[2]*(cos(2.0*k_x)+cos(2.0*k_y)+cos(2.0*(k_x+k_y)))\
    + 2.0*t[3]*(cos((k_x+2.0*k_y))+cos((2.0*k_x+k_y))+cos((k_x-k_y))) \
    + 2.0*t[4]*(cos(2.0*(k_x+2.0*k_y))
                + cos(2.0*(k_y-2.0*(k_x+k_y)))+cos(2.0*(-(k_x+k_y)+2.0*k_x)))
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

def main():
  '''main function'''
  print('test')

if __name__ == '__main__':
  main()
