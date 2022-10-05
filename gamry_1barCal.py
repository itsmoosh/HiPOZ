import os

import matplotlib.cm
import numpy as np
from glob import glob
import logging as log
from datetime import datetime as dtime
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import get_cmap
from gamryTools import Solution, readFloat, PlotY, PlotZ, PlotZvsf, PlotPhasevsf, AddTicksX, AddTicksY, readFloat

date = '20220915'
cmapName = 'viridis'

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}\usepackage{siunitx}\usepackage{upgreek}\sisetup{round-mode=places,scientific-notation=true,round-precision=2}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'STIXGeneral'
cmap = get_cmap(cmapName)
outFigName = 'GamryCal'
figSize = (6,6)
xtn = 'pdf'
PLOT_AIR = False

gamryDTstr = r'%m/%d/%Y-%I:%M %p'



# Import data
add = None
gamryFiles = glob(os.path.join('calData', date, '*Temp*.txt'))

nSweeps = np.size(gamryFiles)
cals = np.empty(nSweeps, dtype=object)

for i, file in enumerate(gamryFiles):
    cals[i] = Solution()
    with open(file) as f:
        f.readline()  # Skip intro line
        cals[i].time = dtime.strptime(f.readline()[:-1], gamryDTstr)  # Measurement time
        cals[i].T_K = readFloat(f.readline())  # Temp
        cals[i].P_MPa = readFloat(f.readline())  # Pressure
        cals[i].descrip = f.readline()  # Text description
        cals[i].Vdrive_V = readFloat(f.readline())  # Driving voltage
        cals[i].fStart_Hz = readFloat(f.readline())  # Driving voltage
        cals[i].fStop_Hz = readFloat(f.readline())  # Driving voltage
        cals[i].nfSteps = int(readFloat(f.readline()))  # Number of f steps

    if 'DIwater' in cals[i].descrip:
        cals[i].comp = 'Pure H2O'
        cals[i].sigmaStd_Sm = 0
        cals[i].legLabel = r'$\approx0$'
    elif 'Air' in cals[i].descrip:
        cals[i].comp = 'Air'
        cals[i].sigmaStd_Sm = np.nan
        cals[i].legLabel = 'Air'
    else:
        cals[i].comp = 'KCl'
        if 'uScm' in cals[i].descrip:
            cals[i].sigmaStd_Sm = float(cals[i].descrip.split(':')[-1].split('uScm')[0]) / 1e4
        elif 'mScm' in cals[i].descrip:
            cals[i].sigmaStd_Sm = float(cals[i].descrip.split(':')[-1].split('mScm')[0]) / 10
        elif 'Sm' in cals[i].descrip:
            cals[i].sigmaStd_Sm = float(cals[i].descrip.split(':')[-1].split('Sm')[0])
        else:
            cals[i].sigmaStd_Sm = np.nan
        cals[i].legLabel = f'{cals[i].sigmaStd_Sm:.4f}'

    _, cals[i].f_Hz, Zabs_ohm, Phi_ohm = np.loadtxt(file, skiprows=10, unpack=True)
    # Phase negated for convenience looking at plots in 1st quadrant
    cals[i].Z_ohm = Zabs_ohm * np.exp(-1j*np.deg2rad(Phi_ohm))

    if not PLOT_AIR and cals[i].comp == 'Air':
        gamryFiles.remove(file)
        cals.delete(cals[i])
        add = '_noAir'

listStandards = np.array([cal.sigmaStd_Sm for cal in cals])
uniqueListStandards = np.unique(listStandards)
standardsSort = np.argsort(uniqueListStandards)
for cal in cals:
    colorLabel = standardsSort[cal.sigmaStd_Sm == uniqueListStandards][0]
    cal.color = cmap(colorLabel / np.max(standardsSort))


# Sort list based on conductivity
iSort = np.argsort(listStandards)
cals = cals[iSort]

sigmaStdBottle_Sm = np.array([0.002179, 0.0081, 0.0423, 0.1926, 0.26, 1.401733, 7.661103]) # using fits based on actual T soln was measured at
# ('sigma_standards_T_table.xlsx')
sigmaStdlist_Sm = np.array([sol.sigmaStd_Sm for sol in cals])  # value on each bottle at 25 deg C, not using fits for temp dependence
iMaxSigma = np.where([sol.sigmaStd_Sm == np.max(sigmaStdlist_Sm[np.isfinite(sigmaStdlist_Sm)]) for sol in cals])[0][0]

Kcell_pm = np.real(cals[iMaxSigma].Z_ohm[0]) * np.max(sigmaStdBottle_Sm)

Rtick_ohm = Kcell_pm / sigmaStdBottle_Sm

PlotZ(cals, figSize, outFigName, xtn, Rtick_ohm)
PlotY(cals, figSize, outFigName, xtn, 1/Rtick_ohm)
PlotZvsf(cals, figSize, outFigName, xtn, Rtick_ohm)
PlotPhasevsf(cals, figSize, outFigName, xtn)


# 1/Z_cell = 1/R + i*omega*C -- Pan et al. (2021): https://doi.org/10.1029/2021GL094020
