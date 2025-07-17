# Ilseung Han
#   16.07.2025

import os

def dustopac(ndust = 1, model = 'astrosil'):
    os.system('rm -rf dustopac.inp')
    with open('dustopac.inp', 'w') as f:
        f.write('2\n')                  # format number (always 2!)
        f.write('%d' % ndust)           # the number of dust species
        for i in range(ndust):
            f.write('\n-')
            f.write('\n1')              # input file style (1: dustkappa_*.inp)
            f.write('\n0')              # normal thermal grains (always 0!)
            f.write('\n%s' % model[i])  # the name of dust species
