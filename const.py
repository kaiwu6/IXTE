'''Global pnysical and lattica constant for MoS2 model '''
import math
import numpy as np

GS = 4.0
KAPPA_PH = 19.  # 34.0
TEMP = 0.026  # temperature in EV
KNUMB = 6
NT = KNUMB * KNUMB

pi = math.pi
V = 0.1 / NT

HBAR = 1.05e-34     #Planck constant
EV = 1.6e-19        #electron unit
KB = 1.38e-23       #Boltzmann constant

A_LAT = 3.16e-10    #a-axis of MoS2
C_LAT = 12.30e-10   #c-axis of MoS2
UNIT_CELL = A_LAT ** 2.0 * C_LAT * np.sqrt(3.0) / 2.0
R_TAU = 1.05e-14 * 43.6 * 5.4  # electron*racial
VELOCITY_UNIT = EV * A_LAT / HBAR / (np.sqrt(3.0) / 2.0)
    # eV is for switch EV to J ?

CONDUCT_UNIT = (VELOCITY_UNIT ** 2 * R_TAU /
                UNIT_CELL / EV * EV ** 2) / NT / 1.e6
SEEBECK_UNIT = KB / EV / TEMP * 1.e6
THERMAL_CONDUCT_UNIT = (VELOCITY_UNIT ** 2 * R_TAU / UNIT_CELL * KB) / TEMP / NT
# K-Grid for calculation
K_GRID = [[float(k_x) / KNUMB, float(k_y) / KNUMB] for k_x in range(KNUMB)
          for k_y in range(KNUMB)]

EVALUATION_POINTS = 4  #for transport calculation
