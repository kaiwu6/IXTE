'''main'''
#import os
import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt
import BSEquation as BSE
from tight_binding import band_fitting
from transport_module import get_band_data, GetStats, correlation
from transport_module import thermal_result
from const import EVALUATION_POINTS, TEMP, NT_TRANS
from const import KAPPA_PH
#matplotlib.use('TkAgg')

def main():
  '''main program: calculate the transportation'''
  print('link start')
  dir_mos2 = '/Users/wk/Dropbox/2015rp/TMDC_Exciton/Tight-binding Model/MoS2/'
  seedname = 'MoS2'
  ix_band = BSE.compute_ixbands(dir_mos2, seedname)
  t_1 = band_fitting(ix_band[0])
  t_2 = band_fitting(ix_band[1])
  np.savetxt("t_1.dat", t_1)
  np.savetxt("t_2.dat", t_2)
  print("t_1: ", t_1)
  print("t_2: ", t_2)
  print("tight binding model completed.")
  print("start to calculate the transportations")
  ix_band1, ix_velocity1 = get_band_data(t_1)
  ix_band2, ix_velocity2 = get_band_data(t_2)

# Prepare evaluation points for the thermal coefficient:
  mu_points = list(map(lambda x: -TEMP * np.log(2.0*(1.0+x)),
                       range(EVALUATION_POINTS)))  # -np.log(fugacity)*temp
  print('mu_points:', mu_points)
  ex_L0 = []
  ex_L1 = []
  ex_L2 = []
  ex_density = []
  for mu_ in mu_points:
    print("mu: ", mu_)
    ex_n_1, ex_dn_1 = GetStats().bose_(ix_band1, mu_, TEMP)
    ex_n_2, ex_dn_2 = GetStats().bose_(ix_band2, mu_, TEMP)
    density = (ex_n_1.sum() + ex_n_2.sum()) / NT_TRANS
    print("density: ", density)

    cor_result1 = correlation(ix_band1-mu_, ix_velocity1, ex_dn_1)
    cor_result2 = correlation(ix_band2-mu_, ix_velocity2, ex_dn_2)
    ex_L0.append(cor_result1[0] + cor_result2[0])
    ex_L1.append(cor_result1[1] + cor_result2[1])
    ex_L2.append(cor_result1[2] + cor_result2[2])
    ex_density.append(density)

  # results = list(map(thermal_result, ex_L0, ex_L1, ex_L2))
  results = map(lambda x, y, z: thermal_result(x, y, z, KAPPA_PH),
                ex_L0, ex_L1, ex_L2)

  Lorenz, conduct, seebeck, thermal_conduct, zT, power_factor = list(zip(*results))
  print("Lorenz: ", Lorenz)
  print("zT: ", zT)
  print("seebeck: ", seebeck)
  print("conduct: ", conduct)
  print("thermal: ", thermal_conduct)
  print("compute finished")

if __name__ == "__main__":
  main()
  