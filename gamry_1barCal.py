import os
import numpy as np
from glob import glob
import logging as log
from datetime import datetime as dtime
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}\usepackage{siunitx}\usepackage{upgreek}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'STIXGeneral'
outFigName = 'GamryCal'
figSize = (6,6)
xtn = 'pdf'
PLOT_AIR = True

gamryDTstr = r'%m/%d/%Y-%I:%M %p'

class Solution: 
    def __init__(self):
        comp = None  # Solute composition
        w_ppt = None  # Solute mass concentration in g/kg
        P_MPa = None  # Chamber pressure of measurement in MPa
        T_K = None  # Temperature of measurement in K
        Vdrive_V = None  # Driving voltage in V
        fStart_Hz = None  # Start frequency in Hz
        fStop_Hz = None  # Stop frequency in Hz
        nfSteps = None  # Number of frequency steps
        time = None  # Measurement start time
        sigmaStd_Sm = None  # Known conductivity of standard solutions (as applicable)
        descrip = None  # Text description (as applicable)

        # Outputs
        f_Hz = None  # Frequency values of Gamry sweep measurements in Hz
        Z_ohm = None  # Complex impedance values of Gamry sweep measurements in ohm 
        sigma_Sm = None  # DC electrical conductivity in S/m


def AddTicks(Rticks, lineList, ax):

    defaultTicks = plt.xticks()[0]
    nDefaultTicks = np.size(defaultTicks)
    allTicks = np.concatenate((defaultTicks, Rticks))
    plt.xticks(allTicks)
    Colors = [line.get_color() for line in lineList]
    [t.set_color(color) for t, color in zip(ax.xaxis.get_ticklabels()[nDefaultTicks:], Colors)]
    [t.set_color(color) for t, color in zip(ax.xaxis.get_ticklines()[nDefaultTicks:], Colors)]


def PlotZ(cals, figSize, outFigName, xtn, Rticks, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'$\mathrm{Re}\{Z\}$ ($\Omega$)')
    ax.set_ylabel(r'$-\mathrm{Im}\{Z\}$ ($\Omega$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- impedance')
    ax.set_xscale('log')
    ax.set_yscale('log')

    for thisCal in cals:
        if np.isnan(thisCal.sigmaStd_Sm):
            legLabel = 'Air'
        else:
            legLabel = f'{thisCal.sigmaStd_Sm:.4f}'

    lineList = [ax.plot(np.real(thisCal.Z_ohm), np.imag(thisCal.Z_ohm), label=legLabel)[0] for thisCal in cals]

    AddTicks(Rticks, lineList, ax)
    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Z{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.debug(f'Gary calibration plot saved to file: {outfName}')
    plt.close()


def PlotY(cals, figSize, outFigName, xtn, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'$\mathrm{Re}\{Y\}$ ($\mho$)')
    ax.set_ylabel(r'$-\mathrm{Im}\{Y\}$ ($\mho$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- conductance')
    ax.set_xscale('log')
    ax.set_yscale('log')

    for thisCal in cals:
        if np.isnan(thisCal.sigmaStd_Sm):
            legLabel = 'Air'
        else:
            legLabel = f'{thisCal.sigmaStd_Sm:.4f}'
        ax.plot(1/np.real(thisCal.Z_ohm), 1/np.imag(thisCal.Z_ohm), label=legLabel)

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Y{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.debug(f'Gary calibration plot saved to file: {outfName}')
    plt.close()


def PlotZvsf(cals, figSize, outFigName, xtn, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Frequency $f$ ($\si{Hz}$)')
    ax.set_ylabel(r'Impedance $|Z|$ ($\Omega$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- impedance spectrum')
    ax.set_xscale('log')
    ax.set_yscale('log')

    for thisCal in cals:
        if np.isnan(thisCal.sigmaStd_Sm):
            legLabel = 'Air'
        else:
            legLabel = f'{thisCal.sigmaStd_Sm:.4f}'
        ax.plot(thisCal.f_Hz, np.abs(thisCal.Z_ohm), label=legLabel)

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Zvsf{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.debug(f'Gary calibration plot saved to file: {outfName}')
    plt.close()


def readFloat(line):
    return float(line.split(':')[-1])


# Import data
airFname = os.path.join('calData', '20220621-1457_Temp301K.txt')
add = None
gamryFiles = glob(os.path.join('calData', '*.txt'))
if not PLOT_AIR and airFname in gamryFiles:
    gamryFiles.remove(airFname)
    add = '_noAir'

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
    elif 'Air' in cals[i].descrip:
        cals[i].comp = 'Air'
        cals[i].sigmaStd_Sm = np.nan
    else:
        cals[i].comp = 'KCl'
        cals[i].sigmaStd_Sm = float(cals[i].descrip.split(':')[-1].split('uScm')[0]) / 1e4

    _, cals[i].f_Hz, Zabs_ohm, Phi_ohm = np.loadtxt(file, skiprows=10, unpack=True)
    # Phase negated for convenience looking at plots in 1st quadrant
    cals[i].Z_ohm = Zabs_ohm * np.exp(-1j*np.deg2rad(Phi_ohm))

# Sort list based on conductivity
iSort = np.argsort(np.array([cal.sigmaStd_Sm for cal in cals]))
cals = cals[iSort]

sigmaStdBottle_Sm = np.array([0.002179, 0.0081, 0.0423, 0.1926, 0.26, 1.401733, 7.661103]) # using fits based on actual T soln was measured at
# ('sigma_standards_T_table.xlsx')
sigmaStdlist_Sm = np.array([sol.sigmaStd_Sm for sol in cals])  # value on each bottle at 25 deg C, not using fits for temp dependence
iMaxSigma = np.where([sol.sigmaStd_Sm == np.max(sigmaStdlist_Sm[np.isfinite(sigmaStdlist_Sm)]) for sol in cals])[0][0]

Kcell_pm = np.real(cals[iMaxSigma].Z_ohm[-1]) * np.max(sigmaStdBottle_Sm)

Rtick_ohm = Kcell_pm / sigmaStdBottle_Sm

PlotZ(cals, figSize, outFigName, xtn, Rtick_ohm)
PlotY(cals, figSize, outFigName, xtn, 1 / Rtick_ohm)
PlotZvsf(cals, figSize, outFigName, xtn, Rtick_ohm)


# 1/Z_cell = 1/R + i*omega*C -- Pan et al. (2021): https://doi.org/10.1029/2021GL094020
