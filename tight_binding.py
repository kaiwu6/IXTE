# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 19:31:19 2016
@author: wk
"""
from __future__ import print_function
#from functools import reduce
import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import rc
from const import pi, KNUMB
#from BSEquation import *
from scipy.optimize import curve_fit

cos = np.cos
def band_expansion(k, *t):
  '''hopping terms'''
  x, y = k
  f_test = t[0]+2.0*t[1]*(cos(2.0*pi*x)+cos(2.0*pi*y)+cos(2.0*pi*(x+y)))\
    + 2.0*t[2]*(cos(2.0*pi*2.0*x)+cos(2.0*pi*2.0*y)+cos(2.0*pi*2.0*(x+y)))\
    + 2.0*t[3]*(cos(2.0*pi*(x+2.0*y))+cos(2.0*pi*(2.0*x+y))+cos(2.0*pi*(x-y))) \
    + 2.0*t[4]*(cos(2.0*pi*(2.0*x+4.0*y))
                + cos(2.0*pi*(2.0*y-4.0*(x+y)))+cos(2.0*pi*(-2.0*(x+y)+4.0*x)))
  return f_test

# def band_terms(k_points):
k_list = []
for j in range(KNUMB):
  for i in range(KNUMB):
    k_list.append([float(i)/KNUMB, float(j)/KNUMB])
klistT = np.transpose(k_list)

def band_fitting(data_source):
  '''curve fitting to get the hopping'''
  if isinstance(data_source, str):
    ka_kb_eng = np.loadtxt(data_source)
    eng_1d = ka_kb_eng[:, 2]
  else:
    eng_1d = data_source
  #eng_2d = eng_1d.reshape(KNUMB, KNUMB)
  weight = data_source + 1e-5
  popt, _ = curve_fit(band_expansion, klistT, eng_1d,
                      p0=np.zeros(5), sigma=weight)
  return popt

'''
data_source_test = '/Users/wk/Dropbox/2015rp/TMDC_Exciton/Tight-binding Model/exband.dat'
'''
'''
def Ft(a, b, i, j):
    kexp = np.exp((float(a*i)/KNUMB + float(b*j)/KNUMB)*2.0*pi*1j)
    return kexp


kplottest = []
for i in range(int(KNUMB/2)):  # Gamma to M
    dim = i
    kplottest.append([float(i)/KNUMB, 0])
for j in range(int(KNUMB/6)):  # M to K
    ka = int(KNUMB/2) - j
    kb = j*2
    dim = ka + kb * KNUMB
    kplottest.append([float(ka)/KNUMB, float(kb)/KNUMB])
for l in range(int(KNUMB/3)):  # K to Gamma
    ka = int(KNUMB/3) - l
    kb = int(KNUMB/3) - l
    dim = ka + kb * KNUMB
    kplottest.append([float(ka)/KNUMB, float(kb)/KNUMB])

kplottest = np.transpose(kplottest)

kplotlist = []
for i in range(int(KNUMB/2)):  # Gamma to M
    dim = i
    kplotlist.append(dim)
for j in range(int(KNUMB/6)):  # M to K
    ka = int(KNUMB/2) - j
    kb = j*2
    dim = ka + kb * KNUMB
    kplotlist.append(dim)
for l in range(int(KNUMB/3)):  # K to Gamma
    ka = int(KNUMB/3) - l
    kb = int(KNUMB/3) - l
    dim = ka + kb * KNUMB
    kplotlist.append(dim)


testingband = [eng_1d[i] for i in kplotlist]

newtesting = band_fit(kplottest, *popt)
fig, axes = plt.subplots()
plt.title('Exciton band at high symmetry lines')
plt.ylabel('$\epsilon_{ex}$(eV)', fontsize=20)
axes.get_xaxis().set_visible(False)
plt.plot(newtesting, 'b')
plt.plot(testingband, 'ro')
plt.axis([0, 48, 0.0, 0.55])
plt.show()
'''
