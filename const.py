'''Global pnysical and lattica constant for MoS2 model '''
import math
import numpy as np

GS = 4.0
KAPPA_PH = 34.0
TEMP = 0.026  # temperature in EV
KNUMB = 6
KNUMB_TRANS = 144
NT = KNUMB * KNUMB
NT_TRANS = KNUMB_TRANS * KNUMB_TRANS

pi = math.pi
V = 0.17 / NT
# 7 layers of hBN: 2.3nm, dielectric const 3.7

HBAR = 1.05e-34     #Planck constant
EV = 1.6e-19        #electron unit
KB = 1.38e-23       #Boltzmann constant

A_LAT = 3.16e-10    #a-axis of MoS2, unit m
C_LAT = 12.30e-10   #c-axis of MoS2, unit m
UNIT_CELL = A_LAT ** 2 * C_LAT * np.sqrt(3.0) / 2.0
R_TAU = 1.05e-14 * 43.6  # electron*racial
VELOCITY_UNIT = EV * A_LAT / HBAR / (np.sqrt(3.0) / 2.0)
    # eV is for switch EV to J ?

CONDUCT_UNIT = (VELOCITY_UNIT ** 2 * R_TAU /
                UNIT_CELL / EV * EV ** 2) / NT_TRANS
SEEBECK_UNIT = KB / EV / TEMP
THERMAL_CONDUCT_UNIT = (VELOCITY_UNIT ** 2 * R_TAU
                        / UNIT_CELL * KB) / TEMP / NT_TRANS

#K_ARRAY with 2pi for calculation
#K_ARRAY_NOPI without 2pi for indexing
K_ARRAY = [[2.0 * pi * k_x / KNUMB, 2.0 * pi * k_y / KNUMB]
           for k_x in range(KNUMB) for k_y in range(KNUMB)]

K_ARRAY_NOPI = [[float(k_x) / KNUMB, float(k_y) / KNUMB]
                for k_x in range(KNUMB) for k_y in range(KNUMB)]

EVALUATION_POINTS = 4  #for transport calculation
