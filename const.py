'''Global pnysical and lattica constant for MoS2 model '''
import numpy as np
import math

HBAR = 1.05e-34     #Planck constant
EV = 1.6e-19        #electron unit
KB = 1.38e-23       #Boltzmann constant

A_LAT = 3.16e-10    #a-axis of MoS2
C_LAT = 12.30e-10   #c-axis of MoS2
UNIT_CELL = A_LAT ** 2.0 * C_LAT * np.sqrt(3.0) / 2.0

KNUMB = 6
NT = KNUMB * KNUMB

pi = math.pi
cos = math.cos
sin = math.sin
sqrt = math.sqrt

V = 0.1 / NT