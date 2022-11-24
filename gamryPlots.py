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

def GetUnique(handles, labels, PAIRED=True):
    outLabels = np.unique(labels)
    if not PAIRED:
        dataLabels = np.where([not 'fit' in label for label in outLabels])[0]
        fitLabels = np.where(['fit' in label for label in outLabels])[0]
        outLabels = np.concatenate((outLabels[dataLabels], outLabels[fitLabels]))
    iUnique = np.array([next(i for i,label in enumerate(labels) if label == uLabel) for uLabel in outLabels])
    outHandles = [handles[i] for i in iUnique]
    outLabels = [labels[i] for i in iUnique]

    return outHandles, outLabels


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


def PlotZ(sols, figSize, outFigName, xtn, Rticks, add=None, LEG_PAIRS=True):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'$\mathrm{Re}\{Z\}$ ($\Omega$)')
    ax.set_ylabel(r'$-\mathrm{Im}\{Z\}$ ($\Omega$)')
    ax.set_title(r'Calibration solution complex impedance')
    ax.set_xscale('log')
    ax.set_yscale('log')

    dotList = [ax.scatter(np.real(sol.Z_ohm), -np.imag(sol.Z_ohm), marker=MS_data, label=f'{sol.legLabel} data', color=sol.color) for sol in sols]
    lineList = [ax.plot(np.real(sol.Z_ohm), -np.imag(sol.Z_ohm), ls=LS_fit, label=f'{sol.legLabel} fit', color=sol.fitColor)[0] for sol in sols]
    if Rticks is not None:
        AddTicksX(Rticks, lineList, ax)

    allHandles, allLabels = ax.get_legend_handles_labels()
    handles, labels = GetUnique(allHandles, allLabels, PAIRED=LEG_PAIRS)
    ax.legend(handles, labels, title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
    plt.tight_layout()
    if add is None:
        addBit = ''
    else:
        addBit = add
    outfName = f'{outFigName}Z{addBit}.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=200)
    log.info(f'Gamry calibration plot saved to file: {outfName}')
    plt.close()


def PlotY(sols, figSize, outFigName, xtn, Rticks, add=None, LEG_PAIRS=True):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'$\mathrm{Re}\{Y\}$ ($\mho$)')
    ax.set_ylabel(r'$-\mathrm{Im}\{Y\}$ ($\mho$)')
    ax.set_title(r'Calibration solution complex conductance')
    ax.set_xscale('log')
    ax.set_yscale('log')

    dotList = [ax.scatter(1/np.real(sol.Z_ohm), -1/np.imag(sol.Z_ohm), marker=MS_data, label=f'{sol.legLabel} data', color=sol.color) for sol in sols]
    lineList = [ax.plot(1/np.real(sol.Zfit_ohm), -1/np.imag(sol.Zfit_ohm), ls=LS_fit, label=f'{sol.legLabel} fit', color=sol.fitColor)[0] for sol in sols]
    if Rticks is not None:
        AddTicksX(Rticks, lineList, ax)

    allHandles, allLabels = ax.get_legend_handles_labels()
    handles, labels = GetUnique(allHandles, allLabels, PAIRED=LEG_PAIRS)
    ax.legend(handles, labels, title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
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


def PlotZvsf(sols, figSize, outFigName, xtn, Rticks, add=None, LEG_PAIRS=True):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Frequency $f$ ($\si{Hz}$)')
    ax.set_ylabel(r'Impedance $|Z|$ ($\Omega$)')
    ax.set_title(r'Calibration solution impedance spectrum')
    ax.set_xscale('log')
    ax.set_yscale('log')

    dotList = [ax.scatter(sol.f_Hz, np.abs(sol.Z_ohm), marker=MS_data, label=f'{sol.legLabel} data', color=sol.color) for sol in sols]
    lineList = [ax.plot(sol.f_Hz, np.abs(sol.Zfit_ohm), ls=LS_fit, label=f'{sol.legLabel} fit', color=sol.fitColor)[0] for sol in sols]
    if Rticks is not None:
        AddTicksY(Rticks, lineList, ax)

    allHandles, allLabels = ax.get_legend_handles_labels()
    handles, labels = GetUnique(allHandles, allLabels, PAIRED=LEG_PAIRS)
    ax.legend(handles, labels, title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
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


def PlotPhasevsf(sols, figSize, outFigName, xtn, add=None, LEG_PAIRS=True):
    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    ax.set_xlabel(r'Frequency $f$ ($\si{Hz}$)')
    ax.set_ylabel(r'Phase ($^\circ$)')
    ax.set_title(r'Calibration solution phase spectrum')
    ax.set_xscale('log')

    dotList = [ax.scatter(sol.f_Hz, np.angle(sol.Z_ohm, deg=True), marker=MS_data, color=sol.color, label=f'{sol.legLabel} data') for sol in sols]
    lineList = [ax.plot(sol.f_Hz, np.angle(sol.Zfit_ohm, deg=True), ls=LS_fit, color=sol.fitColor, label=f'{sol.legLabel} fit')[0] for sol in sols]

    allHandles, allLabels = ax.get_legend_handles_labels()
    handles, labels = GetUnique(allHandles, allLabels, PAIRED=LEG_PAIRS)
    ax.legend(handles, labels, title=r'$\sigma_\mathrm{std}$ ($\si{S/m}$)')
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
    # lineList = [ax.plot(sol.P_MPa, sol.sigma_Sm, label=f'{sol.legLabel} data', color=sol.color, marker='o')[0]  for sol in sols]
    # if Rticks is not None:
    #     AddTicksY(Rticks, lineList, ax)

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


def PlotZfit(sols, figSize, xtn, outFigName=None, LEG_PAIRS=True):
    for sol in sols:
        fig = plt.figure(figsize=figSize)
        grid = GridSpec(1, 1)
        ax = fig.add_subplot(grid[0, 0])
        ax.grid()
        ax.set_axisbelow(True)
        ax.set_xlabel('Re($Z$)')
        ax.set_ylabel('$-$Im($Z$)')
        ax.scatter(np.real(sol.Z_ohm), -np.imag(sol.Z_ohm), marker=MS_data, label=f'{sol.legLabel} data', color=sol.color)
        ax.plot(np.real(sol.Zfit_ohm), -np.imag(sol.Zfit_ohm), ls=LS_fit, label=f'{sol.legLabel} fit', color=sol.fitColor)
        ax.set_xlim(left=0)

        allHandles, allLabels = ax.get_legend_handles_labels()
        handles, labels = GetUnique(allHandles, allLabels, PAIRED=LEG_PAIRS)
        ax.legend(handles, labels)
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

def PlotSigma(allMeas, figSize, outFigName, xtn):
    colorstr = 'gbryk'
    iColor = 0

    fig = plt.figure(figsize=figSize)
    grid = GridSpec(1, 1)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    for ax in ax1, ax2:
        ax.grid()
        ax.set_axisbelow(True)
        ax.set_ylabel('$\sigma$(S/m)')
        ax.set_yscale('log')
    ax1.set_xlabel('$T$(K)')
    ax2.set_xlabel('$P$(MPa)')

    allPs = np.zeros_like(allMeas)
    allTs = np.zeros_like(allMeas)
    allSigs = np.zeros_like(allMeas)
    for i,meas in enumerate(allMeas):
        allPs[i] = [thisMeas.P_MPa for thisMeas in meas]
        allTs[i] = [thisMeas.T_K for thisMeas in meas]
        allSigs[i] = [thisMeas.sigma_Sm for thisMeas in meas]

    for meas in allMeas:
        Plist = np.array([thisMeas.P_MPa for thisMeas in meas])
        # klist = [thisMeas. for thisMeas in meas]
        iSort  = np.argsort(Plist)
        Tlist = np.array([thisMeas.T_K for thisMeas in meas])
        Siglist = np.array([thisMeas.sigma_Sm for thisMeas in meas])
        # for thisMeas in meas:
        #     ax.plot(thisMeas.P_MPa,thisMeas.sigma_Sm, marker='o', markerfacecolor = colorstr[iColor], markeredgecolor  = 'k')
        for thisMeas in meas:
            # ax1.plot(thisMeas.T_K,thisMeas.sigma_Sm, marker='o', markerfacecolor = colorstr[iColor], markeredgecolor  = 'k')
            ax1.plot(thisMeas.T_K,thisMeas.sigma_Sm, marker='o', markerfacecolor = thisMeas.color, markeredgecolor  = 'k')
            ax2.plot(thisMeas.P_MPa,thisMeas.sigma_Sm, marker='o', markerfacecolor = thisMeas.color, markeredgecolor  = 'k')
        iColor += 1

    # ax.set_yscale('log')
    plt.show()
    plt.tight_layout()
    outfName = f'{outFigName}CondvsTP.{xtn}'
    fig.savefig(outfName, format=xtn, dpi=500)
    print(f'Cond vs T plot saved to file: {outfName}')
    plt.close()
