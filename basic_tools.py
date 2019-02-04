'''some general bascis too'''
import math
import numpy as np
from const import KNUMB_TRANS, KNUMB

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
  k_x, k_y = k_grid
  term0 = np.full(k_x.shape, 1.0)
  term1 = 2.0 * (cos(k_x) + cos(k_y) + cos(k_x + k_y))
  term2 = 2.0 * (cos(2.0 * k_x) + cos(2.0 * k_y) + cos(2.0 * (k_x + k_y)))
  term3 = 2.0 * (cos(k_x + 2.0 * k_y) + cos(2.0 * k_x + k_y) + cos(k_x - k_y))
  term4 = 2.0 * (cos(2.0 * k_x + 4.0 * k_y) + cos(2.0 * k_y - 4.0 * (k_x + k_y))
                 + cos(2.0 * (k_x - k_y)))
  return [term0, term1, term2, term3, term4]

def gen_velocity_term():
  '''velocity term grid'''
  k_grid = gen_grid()
  k_x, k_y = k_grid
  term0 = np.zeros(k_x.shape)
  term1 = -2.0 * (sin(k_x) + 0.5 * sin(k_x + k_y) - 0.5 * sin(k_y))
  term2 = -2.0 * (2.0 * sin(2.0*k_x) + sin(2.0*(k_x+k_y)) - sin(2.0*k_y))
  term3 = -3.0 * (sin(2.0 * k_x + k_y) + sin(k_x - k_y))
  term4 = 3.0 * (sin(-4.0*k_x-2.0*k_y) - sin(2.0*(k_x - k_y)))
  return [term0, term1, term2, term3, term4]

def plot_line_short():
  '''symmetry plot: Gamma -> M -> K'''
  plot_line = [i for i in range(int(KNUMB / 2))]  # Gamma to M
  plot_line_2 = list(map(lambda x: 2*KNUMB*x + (int(KNUMB/2)-x),
                         range(int(KNUMB/6))))  # M to K
  plot_line_3 = list(map(lambda y: (KNUMB+1) * (int(KNUMB/3)-y),
                         range(int(KNUMB/3))))  # K to Gamma
  plot_line.extend(plot_line_2)
  plot_line.extend(plot_line_3)
  plot_line = list(map(int, plot_line))
  return plot_line

def plot_line_long():
  '''long symmetry line: Gamma -> M -> K -> Gamma -> K ->M'''
  plot_line = [i for i in range(int(KNUMB / 2))]  # Gamma to M
  plot_line_2 = list(map(lambda x: 2*KNUMB*x + (int(KNUMB/2)-x),
                         range(int(KNUMB/6))))  # M to K
  plot_line_3 = list(map(lambda y: (KNUMB+1) * (int(KNUMB/3)-y),
                         range(int(KNUMB/3) + 1)))  # K to Gamma
  plot_line_4 = list(map(lambda z: (KNUMB + 1) * (KNUMB - z - 1),
                         range(int(KNUMB/3))))  #Gamma to K'
  plot_line_5 = list(map(lambda x: KNUMB * (KNUMB * 2/3 - x - 1)
                         +KNUMB * 2/3 - 4 * (x + 1), range(int(KNUMB / 6))))
                         #K' to M
  plot_line.extend(plot_line_2)
  plot_line.extend(plot_line_3)
  plot_line.extend(plot_line_4)
  plot_line.extend(plot_line_5)
  plot_line = list(map(int, plot_line))
  print(plot_line)
  return plot_line
