'''main'''
#import os
import numpy as np
import BSEquation as BSE
# from const import *
from tight_binding import *
from transport_ex import cal_transport

def main():
  ''' main program: calculate the transportation'''
  print('link start')
  dir_mos2 = '/Users/wk/Dropbox/2015rp/TMDC_Exciton/Tight-binding Model/MoS2/'
  seedname = 'MoS2'

  # compute IX bands
  # todo: improve for the folder to dump the results.
  # question: where is the best place to output the plot? 
  ix_band_1, ix_band_2 = BSE.compute_ixbands(dir_mos2, seedname)
  t_1 = band_fitting(ix_band_1)
  t_2 = band_fitting(ix_band_2)
  #cal_transport(t_1, t_2)

if __name__ == "__main__":
  main()
  