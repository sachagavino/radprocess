# Ilseung Han
#   16.07.2025

import os

def dustkappa(amin, amax, apow, model = 'astrosil', wave = [0.27, 4000], nwave = 200):
    os.system('rm -rf dustkappa_%s.inp' % model)
    os.system('optool %s 1.0 -mie -a %f %f %f -l %f %f %d -radmc'
              % (model, amin, amax, apow, wave[0], wave[1], nwave))
    os.system('mv dustkappa.inp dustkappa_%s.inp' % model)
