'''transport only'''
import numpy as np
from transport_module import TightBindingModel, Get_Stats(), correlation
from transport_module import thermal_result
from const import pi, NT_TRANS, KNUMB_TRANS, TEMP
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

EVALUATION_POINTS = 20


grid_num = NT_TRANS
def gen_grid(grid_num):
  '''generate the meshgrid'''
  lines = np.linspace(0, 2 * pi, grid_num)
  k_x, k_y = np.meshgrid(lines, lines)
  return [k_x, k_y]

def laod_data_from_file(file):
  return 0

def prepare_mu_points(points):
  mu0 = []
  for point in points:
    mu0.append(point)
  return mu0

def main():
  '''load t from files and create the band structure'''
  file1 = 'file1'
  file2 = 'file2'
  t_1 = laod_data_from_file(file1)
  t_2 = laod_data_from_file(file2)

  k_grid = gen_grid(grid_num)
  band_1 = TightBindingModel(k_grid, t_1)
  band_2 = TightBindingModel(k_grid, t_2)
  ix_band_1, ix_velocity_1 = band_1.get_results()
  ix_band_2, ix_velocity_2 = band_2.get_results()

  mu_points = prepare_mu_points(EVALUATION_POINTS)
  for mu_ in mu_points:
    print("mu: ", mu_)
    ix_n_1, ix_dn_1 = GetStats().bose_(ix_band_1, mu_, TEMP)
    ix_n_2, ix_dn_2 = GetStats().bose_(ix_band_2, mu_, TEMP)
    ix_density = (ix_n_1.sum() + ix_n_2.sum()) / NT_TRANS
    print("density: ", ix_density)
    eng_mu_1 = ix_band_1 - mu_
    eng_mu_2 = ix_band_2 - mu_

  ix_correlation_1 = correlation(eng_mu_1, ix_velocity_1, ix_dn_1)
  ix_correlation_2 = correlation(eng_mu_2, ix_velocity_2, ix_dn_2)

  ix_correlation = list(map(lambda x, y: x + y, ix_correlation_1, ix_correlation_2))

  results = map(thermal_result, ix_correlation)
  [Lorenz, conduct, seebeck, thermal_conduct, zT, PF] = list(zip(*results))[:6]
  print("Lb: ", Lorenz)
  print("zT: ", zT)
  print("seebeck: ", seebeck)
  print("conduct: ", conduct)
  print("thermal: ", thermal_conduct)
  print("compute finished")

  plt.figure(1)
  plt.title('Bose Lorenz number')
  plt.xlabel('density(per unit cell)')
  plt.ylabel('Lorenz number')
  plt.xlim([1.e-3, 1.2e-1])
  plt.semilogx(ix_density, Lorenz, 'bs')
  plt.show()

  print("figure 1 done")
  _, ax1 = plt.subplots()
  ax1.set_xlim([1.e-3, 1.2e-1])
  ax1.set_xlabel('density(per unit cell)')
  ax1.semilogx(ix_density, seebeck, 'g^')
  ax1.set_ylabel('seebeck($"\"mu$V/K)')

  ax2 = ax1.twinx()
  ax2.set_xlim([1.e-3, 1.2e-1])
  ax2.semilogx(ix_density, conduct, 'r-', ix_density, thermal_conduct, 'b-')
  ax2.set_ylabel('conduct')
  plt.show()

  fig4, ax1 = plt.subplots()
  ax1.set_xlim([1.e-3, 1.2e-1])
  ax1.set_xlabel('density(per unit cell)')
  ax1.semilogx(ix_density, zT, 'ro')
  ax1.set_ylabel('zT')

  ax2 = ax1.twinx()
  ax2.set_xlim([1.e-3, 1.2e-1])
  ax2.semilogx(ix_density, PF, 'b-')
  ax2.set_ylabel('PF')
  plt.show()



