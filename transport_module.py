# -*- coding: utf-8 -*-
"""Created on Sat Jan  2 19:10:56 2014

@author: wk
"""
import numpy as np
from const import * #pi, KNUMB, GS, CONDUCT_UNIT, SEEBECK_UNIT, THERMAL_CONDUCT_UNIT 
# from numpy import linalg as LA
# from pylab import pcolor
# from functools import reduce
# from scipy.optimize import curve_fit

cos = np.cos
sin = np.sin

class EnergyVelocity():  # need to improve to mesh_grid computation
  ''' exciton band structure '''
  def __init__(self, k_grid, t):
    self.k_grid = k_grid
    self.t = t

  def cal_energy(self):  # take meshgrid as input
    kx = self.k_grid[0]
    ky = self.k_grid[1]
    term0 = np.full(kx.shape, 1.0)
    term1 = 2.0 * (cos(kx) + cos(ky) + cos(kx + ky))
    term2 = 2.0 * (cos(2.0 * kx) + cos(2.0 * ky) + cos(2.0 * (kx + ky)))
    term3 = 2.0 * (cos(kx + 2.0 * ky) + cos(2.0 * kx + ky) + cos(kx - ky))
    term4 = 2.0 * (cos(2.0 * kx + 4.0 * ky) + cos(2.0 * ky - 4.0 * (kx + ky)) + cos(2.0 * (kx - ky)))
    hopping_terms = [term0, term1, term2, term3, term4]
    return sum(map(lambda x, y: x * y, self.t, hopping_terms))

  def cal_velocity(self):
    kx = self.k_grid[0]
    ky = self.k_grid[1]
    term0 = np.zeros(kx.shape)
    term1 = 2.0 * (sin(kx) + 0.5 * sin(kx + ky) - 0.5 * sin(ky))
    term2 = 2.0 * (sin(2.0 * kx) + 0.5 * (sin(2.0 * (kx + ky)) - sin(2.0 * ky)))
    term3 = 2.0 * (sin(2.0 * kx + ky) + sin(kx - ky))
    term4 = 2.0 * (sin(4.0 * kx + 2.0 * ky) + sin(2.0 * (kx - ky)))
    velocity_terms = [term0, term1, term2, term3, term4]
    return sum(map(lambda x, y: x * y, self.t, velocity_terms))

  def get_results(self):
    band = self.cal_energy()
    velocity = self.cal_velocity()
    return band, velocity

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
  l0_grid = velocity * velocity * d_density
  l1_grid = eng_mu * l0_grid
  l2_grid = eng_mu * l1_grid
  l0 = sum(sum(l0_grid)) #double check thisi
  print("l0: ", l0.size)
  l1 = sum(sum(l1_grid))
  l2 = sum(sum(l2_grid))
  return l0, l1, l2

def thermal_result(l0_, l1_, l2_):
  Lb = (l0_ * l0_ - l1_ * l1_) / (l0_ * l0_ * TEMP ** 2)
  conduct = l0_ * CONDUCT_UNIT
  seebeck = l1_ / l0_ * SEEBECK_UNIT
  thermal_conduct = (l2_ - (l1_ * l1_) / l0_) * THERMAL_CONDUCT_UNIT
  zT = conduct * seebeck ** 2 * TEMP / \
  (thermal_conduct + 2.0 * KAPPA_PH) * EV / KB
  PF = conduct * seebeck ** 2 * 1.e3
  return [Lb, conduct, seebeck, thermal_conduct, zT, PF]
