'''main'''
#import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import BSEquation as BSE
from tight_binding import band_fitting
from transport_module import EnergyVelocity, Statistics, correlation
from transport_module import thermal_result
from const import pi, KNUMB, EVALUATION_POINTS, TEMP, NT
matplotlib.use('TkAgg')

def generate_k_grid():
  '''Generate the meshgrid'''
  lines = np.linspace(0, 2 * pi, KNUMB) 
  #not necessary to be KNUMB
  k_x, k_y = np.meshgrid(lines, lines)
  return [k_x, k_y]

def main():
  ''' main program: calculate the transportation'''
  print('link start')
  dir_mos2 = '/Users/wk/Dropbox/2015rp/TMDC_Exciton/Tight-binding Model/MoS2/'
  seedname = 'MoS2'
  # todo: improve for the folder to dump the results.
  # question: where is the best place to output the plot?
  ix_band_1, ix_band_2 = BSE.compute_ixbands(dir_mos2, seedname)
  t_1 = band_fitting(ix_band_1)
  t_2 = band_fitting(ix_band_2)
  np.savetxt("t_1.dat", t_1)
  np.savetxt("t_2.dat", t_2)

  k_grid = generate_k_grid()
  band_1 = EnergyVelocity(k_grid, t_1)
  band_2 = EnergyVelocity(k_grid, t_2)
  ex_band_1, ex_velocity_1 = band_1.get_results()
  ex_band_2, ex_velocity_2 = band_2.get_results()
  print(ex_band_1.shape)
# Prepare evaluation points for the thermal coefficient:
  fugacity = list(map(lambda x: 1.0 * (1.0 + x) / 100, range(EVALUATION_POINTS)))
  mu_points = list(map(lambda x: TEMP * np.log(x), fugacity))  # -np.log(fugacity)*temp
  print('mu_points:', mu_points)
  ex_L0 = []
  ex_L1 = []
  ex_L2 = []
  ex_density = []
  for mu_ in mu_points:
    ex_n_1, ex_dn_1 = Statistics.bose_(ex_band_1, mu_, TEMP)
    ex_n_2, ex_dn_2 = Statistics.bose_(ex_band_2, mu_, TEMP)
    density = (sum(ex_n_1) + sum(ex_n_2)) / NT
    eng_mu_1 = ex_band_1 - mu_
    eng_mu_2 = ex_band_2 - mu_

    ex_l0_1, ex_l1_1, ex_l2_1 = correlation(eng_mu_1, ex_velocity_1, ex_dn_1)
    ex_l0_2, ex_l1_2, ex_l2_2 = correlation(eng_mu_2, ex_velocity_2, ex_dn_2)
    print('ex_l0', ex_l0_1)
    ex_L0.append(ex_l0_1 + ex_l0_2)
    ex_L1.append(ex_l1_1 + ex_l1_2)
    ex_L2.append(ex_l2_1 + ex_l2_2)
    ex_density.append(density)

  # results = list(map(thermal_result, ex_L0, ex_L1, ex_L2))
  results = map(thermal_result, ex_L0, ex_L1, ex_L2)

  [Lb, conduct, seebeck, thermal_conduct, zT, PF] = list(zip(*results))[:6]
  print(Lb)
  print("compute finished")

  plt.figure(1)
  plt.title('Bose Lorenz number')
  plt.xlabel('density(per unit cell)')
  plt.ylabel('Lorenz number')
  plt.xlim([1.e-3, 1.2e-1])
  plt.semilogx(ex_density, Lb, 'bs')
  plt.show()

  print("figure 1 done")
  _, ax1 = plt.subplots()
  ax1.set_xlim([1.e-3, 1.2e-1])
  ax1.set_xlabel('density(per unit cell)')
  ax1.semilogx(ex_density, seebeck, 'g^')
  ax1.set_ylabel('seebeck($\mu$V/K)')

  ax2 = ax1.twinx()
  ax2.set_xlim([1.e-3, 1.2e-1])
  ax2.semilogx(ex_density, conduct, 'r-', ex_density, thermal_conduct, 'b-')
  ax2.set_ylabel('conduct')
  plt.show()

  fig4, ax1 = plt.subplots()
  ax1.set_xlim([1.e-3, 1.2e-1])
  ax1.set_xlabel('density(per unit cell)')
  ax1.semilogx(ex_density, zT, 'ro')
  ax1.set_ylabel('zT')

  ax2 = ax1.twinx()
  ax2.set_xlim([1.e-3, 1.2e-1])
  ax2.semilogx(ex_density, PF, 'b-')
  ax2.set_ylabel('PF')
  plt.show()


if __name__ == "__main__":
  main()
  