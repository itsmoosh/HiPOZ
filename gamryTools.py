import os

import matplotlib.cm
import numpy as np
from glob import glob
import logging as log
from datetime import datetime as dtime
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import get_cmap

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
        legLabel = None  # Legend label
        color = None  # Color of lines

        # Outputs
        f_Hz = None  # Frequency values of Gamry sweep measurements in Hz
        Z_ohm = None  # Complex impedance values of Gamry sweep measurements in ohm 
        sigma_Sm = None  # DC electrical conductivity in S/m


def AddTicksX(Rticks, lineList, ax):
        defaultTicks = ax.get_xticks()
        defaultTickLabels = [f'$10^{{{np.log10(tick):.0f}}}$' for tick in defaultTicks]
        RtickLabels = [f'$\\num{{{Rtick}}}$' for Rtick in Rticks]
        nDefaultTicks = np.size(defaultTicks)
        allTicks = np.concatenate((defaultTicks, Rticks))
        allTickLabels = np.concatenate((defaultTickLabels, RtickLabels))

        ax.set_xticks(allTicks)
        ax.set_xticklabels(allTickLabels)
        lineColors = [line.get_color() for line in lineList[::5]]
        lineStyles = [line.set_linestyle('--') for line in lineList[::2]]
        plt.setp(ax.xaxis.get_ticklabels()[nDefaultTicks:], rotation='vertical')
        [plt.setp(tick, color=color) for tick, color in zip(ax.xaxis.get_ticklabels()[nDefaultTicks:], lineColors[1:])]



def AddTicksY(Rticks, lineList, ax):
        defaultTicks = ax.get_yticks()
        defaultTickLabels = [f'$10^{{{np.log10(tick):.0f}}}$' for tick in defaultTicks]
        RtickLabels = [f'$\\num{{{Rtick}}}$' for Rtick in Rticks]
        nDefaultTicks = np.size(defaultTicks)
        allTicks = np.concatenate((defaultTicks, Rticks))
        allTickLabels = np.concatenate((defaultTickLabels, RtickLabels))

        ax.set_yticks(allTicks)
        ax.set_yticklabels(allTickLabels)
        lineColors = [line.get_color() for line in lineList[::5]]
        lineStyles = [line.set_linestyle('--') for line in lineList[::2]]
        plt.setp(ax.yaxis.get_ticklabels()[nDefaultTicks:])
        [plt.setp(tick, color=color) for tick, color in zip(ax.yaxis.get_ticklabels()[nDefaultTicks:], lineColors[1:])]



def PlotZ(cals, figSize, outFigName, xtn, Rticks, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'$\mathrm{Re}\{Z\}$ ($\Omega$)')
    ax.set_ylabel(r'$-\mathrm{Im}\{Z\}$ ($\Omega$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- impedance')
    ax.set_xscale('log')
    ax.set_yscale('log')

    lineList = [ax.plot(np.real(thisCal.Z_ohm), np.imag(thisCal.Z_ohm), label=thisCal.legLabel, color=thisCal.color)[0] for thisCal in cals]
    AddTicksX(Rticks, lineList, ax)

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Z{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    print(f'Gamry calibration plot saved to file: {outfName}')
    plt.show()
    plt.close()


def PlotY(cals, figSize, outFigName, xtn, Rticks, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'$\mathrm{Re}\{Y\}$ ($\mho$)')
    ax.set_ylabel(r'$-\mathrm{Im}\{Y\}$ ($\mho$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- conductance')
    ax.set_xscale('log')
    ax.set_yscale('log')

    lineList = [ax.plot(1/np.real(thisCal.Z_ohm), 1/np.imag(thisCal.Z_ohm), label=thisCal.legLabel, color=thisCal.color)[0] for thisCal in cals]
    AddTicksX(Rticks, lineList, ax)

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Y{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    print(f'Gamry calibration plot saved to file: {outfName}')
    plt.close()


def PlotZvsf(cals, figSize, outFigName, xtn, Rticks, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Frequency $f$ ($\si{Hz}$)')
    ax.set_ylabel(r'Impedance $|Z|$ ($\Omega$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- impedance spectrum')
    ax.set_xscale('log')
    ax.set_yscale('log')

    lineList = [ax.plot(thisCal.f_Hz, np.abs(thisCal.Z_ohm), label=thisCal.legLabel, color=thisCal.color)[0] for thisCal in cals]
    AddTicksY(Rticks, lineList, ax)

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Zvsf{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    print(f'Gamry calibration plot saved to file: {outfName}')
    plt.close()

def PlotPhasevsf(cals, figSize, outFigName, xtn, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Frequency $f$ ($\si{Hz}$)')
    ax.set_ylabel(r'Phase ($^\circ$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- Phase vs Frequency')
    ax.set_xscale('log')

    lineList = [ax.plot(thisCal.f_Hz, -np.angle(thisCal.Z_ohm, deg=True), label=thisCal.legLabel)[0] for thisCal in cals]

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Phasevsf{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    print(f'Gamry calibration plot saved to file: {outfName}')
    plt.close()

def readFloat(line):
    return float(line.split(':')[-1])

def PlotCondvsP(meas, figSize, outFigName, xtn, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Pressure $P$ ($\si{MPa}$)')
    ax.set_ylabel(r'Conductivity $\sigma$ ($\si{S/m}$)')
    ax.set_title(r'Conductivity vs Pressure' + date)
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    Plist = [thisMeas.P_MPa for thisMeas in meas]
    Siglist = [thisMeas.sigma_Sm*1e-4 for thisMeas in meas]
    for thisMeas in meas:
        ax.plot(thisMeas.P_MPa,thisMeas.sigma_Sm, marker='o', markerfacecolor = 'g', markeredgecolor  = 'k')
    # lineList = [ax.plot(thisMeas.P_MPa, thisMeas.sigma_Sm, label=thisMeas.legLabel, color=thisMeas.color, marker='o')[0]  for thisMeas in meas]
    # AddTicksY(Rticks, lineList, ax)

    # ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
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
