import os

import matplotlib.cm
import numpy as np
from glob import glob
import logging as log
from datetime import datetime as dtime
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import get_cmap
from gamryTools import Solution, readFloat, PlotCondvsP

dates = ['20220922','20220923','20220924']
colorstr = 'gbryk'
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

fig = plt.figure(figsize=figSize)
grid = GridSpec(1, 1)
ax = fig.add_subplot(grid[0, 0])
ax.set_xlabel(r'Pressure $P$ ($\si{MPa}$)')
ax.set_ylabel(r'Conductivity $\sigma$ ($\si{S/m}$)')
ax.set_title(r'Conductivity vs Pressure')

# Import data
add = None
iColor = 0
for thisDate in dates:
    gamryFiles = glob(os.path.join('calData', thisDate, '*Temp*.txt'))

    nSweeps = np.size(gamryFiles)
    meas = np.empty(nSweeps, dtype=object)

    for i, file in enumerate(gamryFiles):
        meas[i] = Solution()
        with open(file) as f:
            f.readline()  # Skip intro line
            meas[i].time = dtime.strptime(f.readline()[:-1], gamryDTstr)  # Measurement time
            meas[i].T_K = readFloat(f.readline())  # Temp
            meas[i].P_MPa = readFloat(f.readline())  # Pressure
            meas[i].descrip = f.readline()  # Text description
            meas[i].Vdrive_V = readFloat(f.readline())  # Driving voltage
            meas[i].fStart_Hz = readFloat(f.readline())  # Driving voltage
            meas[i].fStop_Hz = readFloat(f.readline())  # Driving voltage
            meas[i].nfSteps = int(readFloat(f.readline()))  # Number of f steps

        if 'DIwater' in meas[i].descrip:
            meas[i].comp = 'Pure H2O'
            meas[i].sigmaStd_Sm = 0
            meas[i].legLabel = r'$\approx0$'
        elif 'Air' in meas[i].descrip:
            meas[i].comp = 'Air'
            meas[i].sigmaStd_Sm = np.nan
            meas[i].legLabel = 'Air'
        else:
            meas[i].comp = 'KCl'
            if '447uScm' in meas[i].descrip:
                # meas[i].sigmaStd_Sm = float(meas[i].descrip.split(':')[-1].split('uScm')[0]) / 1e4
                meas[i].sigmaStd_Sm = 447/1842 # 10kHz value from after-removal measurement on 20220923, phase at 10kHz is basically zero.
                meas[i].legLabel = '447uScm KCl std'

            elif '2764uScm' in meas[i].descrip:
                meas[i].sigmaStd_Sm = 2764/614 # 10kHz value from after-removal measurement on 20220923, phase at 10kHz is basically zero.
                meas[i].legLabel = '2764uScm KCl std'

        _, meas[i].f_Hz, Zabs_ohm, Phi_ohm = np.loadtxt(file, skiprows=10, unpack=True)
        # Phase negated for convenience looking at plots in 1st quadrant
        meas[i].Z_ohm = Zabs_ohm * np.exp(-1j*np.deg2rad(Phi_ohm))
        # meas[i].sigma_Sm = np.real(np.interp(10e3,meas[i].f_Hz,meas[i].Z_ohm))*meas[i].sigmaStd_Sm)
        # print(i)
        meas[i].sigma_Sm = np.real(np.interp(10e3,meas[i].f_Hz, Zabs_ohm)) * meas[i].sigmaStd_Sm * 1e-4

    listStandards = np.array([cal.sigmaStd_Sm for cal in meas])
    uniqueListStandards = np.unique(listStandards)
    standardsSort = np.argsort(uniqueListStandards)
    for cal in meas:
        colorLabel = standardsSort[cal.sigmaStd_Sm == uniqueListStandards][0]
        cal.color = cmap(colorLabel / np.max(standardsSort))


    # Sort list based on conductivity
    iSort = np.argsort(listStandards)
    meas = meas[iSort]

    Plist = [thisMeas.P_MPa for thisMeas in meas]
    Siglist = [thisMeas.sigma_Sm*1e-4 for thisMeas in meas]
    for thisMeas in meas:
        ax.plot(thisMeas.P_MPa,thisMeas.sigma_Sm, marker='o', markerfacecolor = colorstr[iColor], markeredgecolor  = 'k')
    iColor += 1

ax.grid()
plt.show()
plt.tight_layout()
if add is None:
    addBit = ''
else:
    addBit = add
outfName = f'{outFigName}CondvsP{addBit}.{xtn}'
fig.savefig(outfName, format=xtn, dpi=200)
print(f'Cond vs P plot saved to file: {outfName}')
plt.close()

# 1/Z_cell = 1/R + i*omega*C -- Pan et al. (2021): https://doi.org/10.1029/2021GL094020
