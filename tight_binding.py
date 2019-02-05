# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 19:31:19 2016
@author: wk
"""
import numpy as np
from scipy.optimize import curve_fit
from const import K_ARRAY

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

def band_fitting(data_source):
  '''curve fitting to get the hopping'''
  k_list = np.transpose(K_ARRAY)
  eng_1d = data_source
  weight = data_source + 1e-5
  popt, _ = curve_fit(band_expansion, k_list, eng_1d,
                      p0=np.zeros(5), sigma=weight)
  return popt

def main():
  '''use method to load data from file'''
  bands = np.loadtxt('bands.dat')
  band1 = bands[:, 2]
  band2 = bands[:, 3]
  t_1 = band_fitting(band1)
  t_2 = band_fitting(band2)
  print('t_1: ', t_1)
  print('t_2: ', t_2)

if __name__ == '__main__':
  main()
