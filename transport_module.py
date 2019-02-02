# -*- coding: utf-8 -*-
"""Created on Sat Jan  2 19:10:56 2014

@author: wk
"""
import numpy as np
from const import GS, CONDUCT_UNIT, SEEBECK_UNIT, THERMAL_CONDUCT_UNIT
from const import TEMP, EV, KB
from basic_tools import gen_energy_terms, gen_velocity_term

# from numpy import linalg as LA
# from pylab import pcolor
# from functools import reduce
# from scipy.optimize import curve_fit

def get_band_data(t_hop):
  '''energy and velocity'''
  energy_terms = gen_energy_terms()
  velocity_terms = gen_velocity_term()
  energy = sum(map(lambda x, y: x*y, t_hop, energy_terms))
  velocity = sum(map(lambda x, y: x*y, t_hop, velocity_terms))
  return energy, velocity

class Statistics():
  @staticmethod
  def bose_(eng_grid, mu, temp):
    exp_temp = np.exp(-(eng_grid - mu) / temp)
    n = GS * exp_temp / (1.0 - exp_temp)
    dn = GS * (1 / temp) * exp_temp / ((1.0 - exp_temp) * (1.0 - exp_temp))
    return n, dn

  @staticmethod
  def fermi_(eng_grid, mu, temp):
    exp_temp = np.exp(-(eng_grid - mu) / temp)
    n = GS * exp_temp / (1.0 + exp_temp)
    dn = GS * (1 / temp) * exp_temp / ((1.0 + exp_temp) * (1.0 + exp_temp))
    return n, dn

def correlation(eng_mu, velocity, d_density):
  '''caculate correlations'''
  l0_grid = velocity * velocity * d_density
  l1_grid = eng_mu * l0_grid
  l2_grid = eng_mu * l1_grid
  cor0 = l0_grid.sum()
  cor1 = l1_grid.sum()
  cor2 = l2_grid.sum()
  return cor0, cor1, cor2

def thermal_result(cor0, cor1, cor2, kappa_ph):
  '''get final results'''
  Lorenz = (cor2 * cor0 - cor1 * cor1) / (cor0 * cor0 * TEMP**2)
  conduct = cor0 * CONDUCT_UNIT
  Seebeck = cor1 / cor0 * SEEBECK_UNIT
  tconduct = (cor2 - (cor1 * cor1) / cor0) * THERMAL_CONDUCT_UNIT
  zT = conduct * Seebeck ** 2 * TEMP / \
  (tconduct + 2.0 * kappa_ph) * EV / KB
  power_factor = conduct * Seebeck ** 2
  return [Lorenz, conduct, Seebeck, tconduct, zT, power_factor]
