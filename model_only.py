'''
This is the script for band calculation only
'''
import numpy as np
import BSEquation as BSE
from tight_binding import band_fitting, band_expansion
from const import NT, KNUMB, K_ARRAY_NOPI
from basic_tools import plot_line_long
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

def dump_energy_data(band1, band2):
  '''Dump energy data'''
  k_x, k_y = list(zip(*K_ARRAY_NOPI))
  result = np.c_[k_x, k_y, band1, band2]
  np.savetxt("bands.dat", result)

def main():
  '''calculate the IX band and hopping parameter'''
  print("start to calculate the IX bands")
  dir_mos2 = '/Users/wk/Dropbox/2015rp/TMDC_Exciton/Tight-binding Model/MoS2/'
  seedname = 'MoS2'
  ix_band_1, ix_band_2 = BSE.compute_ixbands(dir_mos2, seedname)
  dump_energy_data(ix_band_1, ix_band_2)

  print('start plot')
  fig1 = plt.figure(figsize=(10, 6))
  a_x = fig1.add_subplot(111)
  new_plot = BSE.HighSymmetryPlot(NT, 1)
  new_plot.plot_one_band_long(a_x, ix_band_1, ix_band_2)
  plt.savefig('band.png')
  plt.show()
  
  #Start to band fitting
  t_1 = band_fitting(ix_band_1)
  t_2 = band_fitting(ix_band_2)
  print("save hopping to files: t_1.dat, t_2.dat")
  np.savetxt("t_1.dat", t_1)
  np.savetxt("t_2.dat", t_2)

  plot_points = plot_line_long()
  band1_fit = band_expansion(plot_points, *t_1)
  band2_fit = band_expansion(plot_points, *t_2)
  
  
  print("completed")

if __name__ == "__main__":
  main()
