from os import system

import sys
if sys.version_info.major > 2:
    from subprocess import STDOUT, run
else:
    from subprocess32 import STDOUT, run

def thermal(noscat=None, nphot_therm=None, nphot_scat=None, setthreads=1, \
        inclfreefree=None, nofreefree=None, inclgascont=None, nogascont=None, \
        verbose=True, timelimit=7200, nice=None):

    if nice != None:
        command="nice -{0:d} radmc3d mctherm".format(nice)
    else:
        command="radmc3d mctherm"

    if (noscat == True):
        command += " noscat"
    if (nphot_therm != None):
        command += " nphot_therm {0:d}".format(nphot_therm)
    if (setthreads != 1):
        command += " setthreads {0:d}".format(setthreads)
    if (inclfreefree == True):
        command += " inclfreefree"
    if (nofreefree == True):
        command += " nofreefree"
    if (inclgascont == True):
        command += " inclgascont"
    if (nogascont == True):
        command += " nogascont"

    if not verbose:
        f = open("thermal/radmc3d.out","w")
        output = run(command.split(" "), cwd = 'thermal/', stdout=f, stderr=f, timeout=timelimit)
        f.close()
    else:
        output = run(command.split(" "), cwd = 'thermal/', stderr=STDOUT, timeout=timelimit)

def localfield(nphot_mono=None, verbose=True, timelimit=7200):

    command="radmc3d mcmono"

    if (nphot_mono != None):
        command += " nphot_mono {0:d}".format(nphot_mono)

    if not verbose:
        f = open("thermal/radmc3d.out","w")
        output = run(command.split(" "), cwd = 'thermal/', stdout=f, stderr=f, timeout=timelimit)
        f.close()
    else:
        output = run(command.split(" "), cwd = 'thermal/', stderr=STDOUT, timeout=timelimit)

def image():
    print('making image here')
