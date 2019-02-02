import numpy as np
import BSEquation as BSE
from tight_binding import band_fitting
from transport_module import TightBindingModel, GetStats, correlation
from transport_module import thermal_result
from const import pi, KNUMB, EVALUATION_POINTS, TEMP, NT, KNUMB_TRANS, NT_TRANS

grid_number = KNUMB_TRANS
def generate_k_grid(grid_number):
  '''Generate the meshgrid'''
  lines = np.linspace(0, 2 * pi, grid_number)
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
  ix_band_1, _ = BSE.compute_ixbands(dir_mos2, seedname)
  t_1 = band_fitting(ix_band_1)
  np.savetxt("t_1.dat", t_1)
  print("t_1: ", t_1)
  print("tight binding model completed.")
  print("start to calculate the transportations")

  k_grid = generate_k_grid(grid_number)
  band_1 = TightBindingModel(k_grid, t_1)
  ex_band_1, ex_velocity_1 = band_1.get_results()
  print(ex_band_1.shape)
# Prepare evaluation points for the thermal coefficient:
  fugacity = list(map(lambda x: 0.2 * (1.0 + x), range(EVALUATION_POINTS)))
  mu_points = list(map(lambda x: -TEMP * np.log(x), fugacity))
  ex_L0 = []
  ex_L1 = []
  ex_L2 = []
  ex_density = []
  for mu_ in mu_points:
    print("mu: ", mu_)
    ex_n_1, ex_dn_1 = GetStats().fermi_(ex_band_1, mu_, TEMP)
    density = ex_n_1.sum() / NT_TRANS
    print("density: ", density)
    eng_mu_1 = ex_band_1 - mu_

    ex_l0_1, ex_l1_1, ex_l2_1 = correlation(eng_mu_1, ex_velocity_1, ex_dn_1)
    print('ex_l0_1:', ex_l0_1)
    print('ex_l1_1:', ex_l1_1)
    print('ex_l2_1:', ex_l2_1)
    ex_L0.append(ex_l0_1)
    ex_L1.append(ex_l1_1)
    ex_L2.append(ex_l2_1)
    ex_density.append(density)

  # results = list(map(thermal_result, ex_L0, ex_L1, ex_L2))
  results = map(thermal_result, ex_L0, ex_L1, ex_L2)

  [Lb, conduct, seebeck, thermal_conduct, zT, PF] = list(zip(*results))[:6]
  print("Lb: ", Lb)
  print("zT: ", zT)
  print("seebeck: ", seebeck)
  print("conduct: ", conduct)
  print("thermal: ", thermal_conduct)
  print("compute finished")
'''
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
'''

if __name__ == "__main__":
  main()
  