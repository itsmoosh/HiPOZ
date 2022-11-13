import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import logging

# Assign logger
log = logging.getLogger('HIPPOS')

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{stix}\usepackage{siunitx}\usepackage{upgreek}\sisetup{round-mode=places,scientific-notation=true,round-precision=2}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'STIXGeneral'
PLOT_AIR = False
tfmt = '%b%d%H%M%S'
uScm = r'\mu{S\,cm^{-1}}'
MS_data = 'o'
LS_fit = '-'


def AddTicksX(Rticks, lineList, ax):
    defaultTicks = ax.get_xticks()
    defaultTickLabels = [f'$10^{{{np.log10(tick):.0f}}}$' for tick in defaultTicks]
    RtickLabels = [f'$\\num{{{Rtick}}}$' for Rtick in Rticks]
    nDefaultTicks = np.size(defaultTicks)
    allTicks = np.concatenate((defaultTicks, Rticks))
    allTickLabels = np.concatenate((defaultTickLabels, RtickLabels))

    ax.set_xticks(allTicks)
    ax.set_xticklabels(allTickLabels)
    lineColors = [line.get_color() for line in lineList]
    #lineStyles = [line.set_linestyle('--') for line in lineList[::2]]
    plt.setp(ax.xaxis.get_ticklabels()[nDefaultTicks:], rotation='vertical')
    [plt.setp(tick, color=color) for tick, color in zip(ax.xaxis.get_ticklabels()[nDefaultTicks:], lineColors)]


def AddTicksY(Rticks, lineList, ax):
    defaultTicks = ax.get_yticks()
    defaultTickLabels = [f'$10^{{{np.log10(tick):.0f}}}$' for tick in defaultTicks]
    RtickLabels = [f'$\\num{{{Rtick}}}$' for Rtick in Rticks]
    nDefaultTicks = np.size(defaultTicks)
    allTicks = np.concatenate((defaultTicks, Rticks))
    allTickLabels = np.concatenate((defaultTickLabels, RtickLabels))

    ax.set_yticks(allTicks)
    ax.set_yticklabels(allTickLabels)
    lineColors = [line.get_color() for line in lineList]
    #lineStyles = [line.set_linestyle('--') for line in lineList[::2]]
    plt.setp(ax.yaxis.get_ticklabels()[nDefaultTicks:], rotation='horizontal')
    [plt.setp(tick, color=color) for tick, color in zip(ax.yaxis.get_ticklabels()[nDefaultTicks:], lineColors)]


def PlotZ(sols, figSize, outFigName, xtn, Rticks, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'$\mathrm{Re}\{Z\}$ ($\Omega$)')
    ax.set_ylabel(r'$-\mathrm{Im}\{Z\}$ ($\Omega$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- impedance')
    ax.set_xscale('log')
    ax.set_yscale('log')

    dotList = [ax.scatter(np.real(sol.Z_ohm), -np.imag(sol.Z_ohm), marker=MS_data, label=sol.legLabel, color=sol.color) for sol in sols]
    lineList = [ax.plot(np.real(sol.Z_ohm), -np.imag(sol.Z_ohm), ls=LS_fit, label=f'{sol.legLabel} fit', color=sol.fitColor)[0] for sol in sols]
    AddTicksX(Rticks, lineList, ax)

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Z{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.info(f'Gamry calibration plot saved to file: {outfName}')
    plt.close()


def PlotY(sols, figSize, outFigName, xtn, Rticks, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'$\mathrm{Re}\{Y\}$ ($\mho$)')
    ax.set_ylabel(r'$-\mathrm{Im}\{Y\}$ ($\mho$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- conductance')
    ax.set_xscale('log')
    ax.set_yscale('log')

    dotList = [ax.scatter(1/np.real(sol.Z_ohm), -1/np.imag(sol.Z_ohm), marker=MS_data, label=sol.legLabel, color=sol.color) for sol in sols]
    lineList = [ax.plot(1/np.real(sol.Zfit_ohm), -1/np.imag(sol.Zfit_ohm), ls=LS_fit, label=f'{sol.legLabel} fit', color=sol.fitColor)[0] for sol in sols]
    AddTicksX(Rticks, lineList, ax)

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Y{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.info(f'Gamry calibration plot saved to file: {outfName}')
    plt.close()

    return


def PlotZvsf(sols, figSize, outFigName, xtn, Rticks, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Frequency $f$ ($\si{Hz}$)')
    ax.set_ylabel(r'Impedance $|Z|$ ($\Omega$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- impedance spectrum')
    ax.set_xscale('log')
    ax.set_yscale('log')

    dotList = [ax.scatter(sol.f_Hz, np.abs(sol.Z_ohm), marker=MS_data, label=sol.legLabel, color=sol.color) for sol in sols]
    lineList = [ax.plot(sol.f_Hz, np.abs(sol.Zfit_ohm), ls=LS_fit, label=f'{sol.legLabel} fit', color=sol.fitColor)[0] for sol in sols]
    AddTicksY(Rticks, lineList, ax)

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Zvsf{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.info(f'Gamry calibration plot saved to file: {outfName}')
    plt.close()

    return


def PlotPhasevsf(sols, figSize, outFigName, xtn, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Frequency $f$ ($\si{Hz}$)')
    ax.set_ylabel(r'Phase ($^\circ$)')
    ax.set_title(r'Calibration solution Gamry sweeps --- Phase vs Frequency')
    ax.set_xscale('log')

    dotList = [ax.scatter(sol.f_Hz, np.angle(sol.Z_ohm, deg=True), marker=MS_data, color=sol.color, label=sol.legLabel) for sol in sols]
    lineList = [ax.plot(sol.f_Hz, np.angle(sol.Zfit_ohm, deg=True), ls=LS_fit, color=sol.fitColor, label=f'{sol.legLabel} fit')[0] for sol in sols]

    ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Phasevsf{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.info(f'Gamry calibration plot saved to file: {outfName}')
    plt.close()

    return


def PlotCondvsP(sols, figSize, outFigName, xtn, add=None):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Pressure $P$ ($\si{MPa}$)')
    ax.set_ylabel(r'Conductivity $\sigma$ ($\si{S/m}$)')
    ax.set_title(r'Conductivity vs Pressure' + date)
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    Plist = [sol.P_MPa for sol in sols]
    Siglist = [sol.sigma_Sm * 1e-4 for sol in sols]
    for sol in sols:
        ax.plot(sol.P_MPa, sol.sigma_Sm, marker='o', markerfacecolor='g', markeredgecolor='k')
    # lineList = [ax.plot(sol.P_MPa, sol.sigma_Sm, label=sol.legLabel, color=sol.color, marker='o')[0]  for sol in sols]
    # AddTicksY(Rticks, lineList, ax)

    # ax.legend(title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}CondvsP{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.info(f'Cond vs P plot saved to file: {outfName}')
    plt.close()

    return


def PlotZfit(sols, figSize, xtn, outFigName=None):
    for sol in sols:
        fig = plt.figure(figsize=figSize)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        ax.grid()
        ax.set_axisbelow(True)
        ax.set_xlabel('Re($Z$)')
        ax.set_ylabel('$-$Im($Z$)')
        ax.scatter(np.real(sol.Z_ohm), -np.imag(sol.Z_ohm), marker=MS_data, label=sol.legLabel, color=sol.color)
        ax.plot(np.real(sol.Zfit_ohm), -np.imag(sol.Zfit_ohm), ls=LS_fit, label=f'{sol.legLabel} fit', color=sol.fitColor)
        ax.set_xlim(left=0)
        plt.legend()
        tstr = sol.time.strftime(tfmt)
        if outFigName is None:
            thisOutFigName = f'{sol.lbl_uScm}uScm_{tstr}'
        else:
            thisOutFigName = outFigName
        ax.set_title(f'Nyquist plot for ${sol.lbl_uScm}\,\si{{{uScm}}}$ at {tstr}, $K_\mathrm{{cell}}={sol.Kcell_pm:.2f}\,\si{{m^{{-1}}}}$')

        outfName = f'{thisOutFigName}Nyquist.{xtn}'
        fig.savefig(outfName, format=xtn, dpi=200)
        log.info(f'Nyquist plot saved to file: {outfName}')
        plt.close()

    return
