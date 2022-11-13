import os, sys

import numpy as np
from glob import glob
import logging
from datetime import datetime as dtime
from gamryTools import Solution, CalStdFit
from gamryPlots import PlotY, PlotZ, PlotZvsf, PlotPhasevsf, PlotZfit

# Assign logger
log = logging.getLogger('HIPPOS')
stream = logging.StreamHandler(sys.stdout)
stream.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
log.setLevel(logging.DEBUG)
log.addHandler(stream)

date = '20220915'
circType = 'CPE'  # Options are 'CPE', 'RC', and 'RC-R'. If desired, a circuit string can be entered here instead.
initial_guess = None  # Required when circType is not in the above list. Ignored otherwise.
cmapName = 'viridis'
outFigName = 'GamryCal'
figSize = (6,6)
xtn = 'pdf'
PLOT_AIR = False

# Import data
add = None
gamryFiles = glob(os.path.join('calData', date, '*Temp*.txt'))

nSweeps = np.size(gamryFiles)
cals = np.empty(nSweeps, dtype=object)
calStd = CalStdFit(interpMethod='cubic')

for i, file in enumerate(gamryFiles):
    cals[i] = Solution(cmapName)
    cals[i].loadFile(file)

    if not np.isnan(cals[i].sigmaStd_Sm):
        cals[i].sigmaStdCalc_Sm = calStd(cals[i].T_K, lbl_uScm=cals[i].lbl_uScm)
    else:
        cals[i].sigmaStdCalc_Sm = 1e-8  # Default air conductivity
    cals[i].FitCircuit(circType=circType, initial_guess=initial_guess)
    cals[i].Kcell_pm = cals[i].sigmaStdCalc_Sm * cals[i].Rcalc_ohm

    if not PLOT_AIR and cals[i].comp == 'Air':
        gamryFiles.remove(file)
        cals.delete(cals[i])
        add = '_noAir'

# Sort list based on conductivity
listStandards = np.array([cal.sigmaStd_Sm for cal in cals])
iSort = np.argsort(listStandards)
cals = cals[iSort]

#Tcal_C = np.array([22.8, 22.3, 21.7, 21.5, 23.215, 21.8, 21.9])
#sigmaStdBottle_Sm = np.array([calStd(Tcal_C + 273.15)])
#sigmaStdlist_Sm = np.array([sol.sigmaStd_Sm for sol in cals])  # value on each bottle at 25 deg C, not using fits for temp dependence

PlotZfit(cals, figSize, xtn)

Rtick_ohm = np.array([cal.Rcalc_ohm for cal in cals])
PlotZ(cals, figSize, outFigName, xtn, Rtick_ohm)
PlotY(cals, figSize, outFigName, xtn, 1/Rtick_ohm)
PlotZvsf(cals, figSize, outFigName, xtn, Rtick_ohm)
PlotPhasevsf(cals, figSize, outFigName, xtn)

