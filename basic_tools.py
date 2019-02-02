'''some general bascis too'''
import math
import numpy as np
from const import KNUMB_TRANS

PI = math.pi
cos = np.cos
sin = np.sin

def gen_grid():
  '''generate grid'''
  lines = np.linspace(0, 2.0*PI, KNUMB_TRANS)
  k_x, k_y = np.meshgrid(lines, lines)
  return k_x, k_y

def gen_energy_terms():
  '''energy term grid'''
  k_grid = gen_grid()
  kx, ky = k_grid
  term0 = np.full(kx.shape, 1.0)
  term1 = 2.0 * (cos(kx) + cos(ky) + cos(kx + ky))
  term2 = 2.0 * (cos(2.0 * kx) + cos(2.0 * ky) + cos(2.0 * (kx + ky)))
  term3 = 2.0 * (cos(kx + 2.0 * ky) + cos(2.0 * kx + ky) + cos(kx - ky))
  term4 = 2.0 * (cos(2.0 * kx + 4.0 * ky) + cos(2.0 * ky - 4.0 * (kx + ky))
                 + cos(2.0 * (kx - ky)))
  return [term0, term1, term2, term3, term4]

def gen_velocity_term():
  '''velocity term grid'''
  k_grid = gen_grid()
  kx, ky = k_grid
  term0 = np.zeros(kx.shape)
  term1 = -2.0 * (sin(kx) + 0.5 * sin(kx + ky) - 0.5 * sin(ky))
  term2 = -2.0 * (2.0 * sin(2.0 * kx) + sin(2.0 * (kx + ky)) - sin(2.0 * ky))
  term3 = -3.0 * (sin(2.0 * kx + ky) + sin(kx - ky))
  term4 = 3.0 * (sin(-4.0*kx-2.0*ky) - sin(2.0*(kx - ky)))
  return [term0, term1, term2, term3, term4]