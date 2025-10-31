# Ilseung Han
#   03.04.2025
#   21.07.2025

import numpy as np
import astropy.units as u
import astropy.constants as co

    # Note that
    #   the stellar position was calculated from the corner in RAMSES;
    #   and this position was calculated from the center in POLARIS/RADMC-3D.

x = 7.9414097559E-02        # the x position of the central star in pc  (sink_00100.csv)
y = 8.7052647954E-02        # the y position of the central star in pc  (sink_00100.csv)
z = 7.9216680912E-02        # the z position of the central star in pc  (sink_00100.csv)
b = 0.169154432522779E+00   # the size of the simulation in pc          (info_00100.txt)
l = 0.308000000000000E+19   # the 1-pc length of the simulation in cm   (info_00100.txt)
l = l * 1e-02               # the 1-pc length of the simulation in m

xpos = (x - b/2) * l        # the x poistion of the central star in m   (POLARIS)
ypos = (y - b/2) * l        # the y position of the central star in m   (POLRAIS)
zpos = (z - b/2) * l        # the z position of the central star in m   (POLARIS)
print('%10.2e (m); %8.2f (au)' % (xpos, xpos / co.au.value))
print('%10.2e (m); %8.2f (au)' % (ypos, ypos / co.au.value))
print('%10.2e (m); %8.2f (au)' % (zpos, zpos / co.au.value))
