# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 19:31:19 2016

@author: wk
"""
from __future__ import print_function
import numpy as np
from functools import reduce
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import rc
from BSEmodule import *
from transport_module import *
from const import *

knumb = KNUMB

# def band_terms(k_points):

# generate finite k grid in BZ
k_points = k_grid_gen(knumb)


klist = []
for j in range(knumb):
    for i in range(knumb):        
        klist.append([float(i)/knumb, float(j)/knumb])

dir = '/Users/wk/Dropbox/2015rp/TMDC_Exciton/Tight-binding Model/'
ka_kb_eng = np.loadtxt(dir+'exband.dat')
eng_1 = ka_kb_eng[:,2]
xeng_2d = eng_1.reshape((knumb, knumb))

def Ft(a,b,i,j):
    kexp = np.exp((float(a*i)/knumb + float(b*j)/knumb)*2.0*pi*1j)
    return kexp
  
'''
tft = []   
tshort = [] 
for a in range(int(knumb)):
    for b in range(int(knumb)):
        teff = 0
        for i in range(knumb):
            for j in range(knumb):
                teff += xeng_2d[i][j]*Ft(a,b,i,j)
        teff = teff/Nt
        tft.append(teff.real)
        if abs(teff) >0.005:
            tshort.append([a,b,teff.real])
#short range with 00, 01, 10, 11, 12, 21, 24, 42? Is this enough?
#print(tshort)
'''

klistT = np.transpose(klist)


def band_fit(k, *t):
    x,y=k
    f = t[0]+2.0*t[1]*(cos(2.0*pi*x)+cos(2.0*pi*y)+cos(2.0*pi*(x+y)))\
        +2.0*t[2]*(cos(2.0*pi*2.0*x)+cos(2.0*pi*2.0*y)+cos(2.0*pi*2.0*(x+y)))\
        +2.0*t[3]*(cos(2.0*pi*(x+2.0*y))+cos(2.0*pi*(2.0*x+y))+cos(2.0*pi*(x-y))) \
        +2.0*t[4]*(cos(2.0*pi*(2.0*x+4.0*y))+cos(2.0*pi*(2.0*y-4.0*(x+y)))+cos(2.0*pi*(-2.0*(x+y)+4.0*x))) \
        #+2.0*t[5]*(cos(2.0*pi*(x+2.0*(x+y)))+cos(2.0*pi*(y-2.0*x))+cos(2.0*pi*(-(x+y)-2.0*y))
        #+cos(2.0*pi*(2.0*x+(x+y))+cos(2.0*pi*(2.0*y-x))+cos(2.0*pi*(-2.0*(x+y)-y))))\
    return f
       
weight = eng_1

popt, pcov = curve_fit(band_fit, klistT, eng_1, p0 =np.zeros(5), sigma = weight)
np.savetxt('tbt.dat', popt)


kplottest = []
for i in range(int(knumb/2)): #Gamma to M
    dim = i
    kplottest.append([float(i)/knumb,0])
for j in range(int(knumb/6) ): #M to K
    ka = int(knumb/2) - j
    kb = j*2
    dim = ka + kb * knumb
    kplottest.append([float(ka)/knumb,float(kb)/knumb])
for l in range(int(knumb/3)): #K to Gamma
    ka = int(knumb/3) - l
    kb = int(knumb/3) - l
    dim = ka + kb *knumb
    kplottest.append([float(ka)/knumb,float(kb)/knumb])

kplottest = np.transpose(kplottest)

kplotlist = []
for i in range(int(knumb/2)): #Gamma to M
    dim = i
    kplotlist.append(dim)
for j in range(int(knumb/6) ): #M to K
    ka = int(knumb/2) - j
    kb = j*2
    dim = ka + kb * knumb
    kplotlist.append(dim)
for l in range(int(knumb/3)): #K to Gamma
    ka = int(knumb/3) - l
    kb = int(knumb/3) - l
    dim = ka + kb *knumb
    kplotlist.append(dim)


testingband =[eng_1[i] for i in kplotlist]

newtesting = band_fit(kplottest,*popt)
fig, axes= plt.subplots()
plt.title('Exciton band at high symmetry lines')
plt.ylabel('$\epsilon_{ex}$(eV)', fontsize = 20)
axes.get_xaxis().set_visible(False)
plt.plot(newtesting, 'b')
plt.plot(testingband, 'ro')
plt.axis([0,48,0.0,0.55])
plt.show()