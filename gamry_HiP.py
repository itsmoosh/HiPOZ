import os, sys, re
import numpy as np
from glob import glob
import logging
from datetime import datetime as dtime

from gamryTools import Solution, CalStdFit
from gamryPlots import PlotY, PlotZ, PlotZvsf, PlotPhasevsf, PlotZfit, PlotSigma

from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Assign logger
log = logging.getLogger('HIPPOS')
stream = logging.StreamHandler(sys.stdout)
stream.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
log.setLevel(logging.DEBUG)
log.addHandler(stream)

dates = ['20220922','20220923','20220924','20221010','20221011']
# dates = ['20221016']; # values read at 5e-4S/m
# dates = ['20221014']; # using the 1bar std yield around 5 S/m not 8. THere's a zero conductivity point at high pressure
# dates = ['20221010','20221011']
# dates = ['20220922','20220923']#,'20220924','20221010','20221011']
# dates = ['20220924']#,'20221010','20221011']
circType = 'CPE'  # Options are 'CPE', 'RC', and 'RC-R'. If desired, a circuit string can be entered here instead.
initial_guess = None  # Required when circType is not in the above list. Ignored otherwise.
cmapName = 'viridis'
outFigName = 'GamryMeas'
figSize = (6,6)
xtn = 'png'
PLOT_AIR = False

# Import data
add = None
nDates = np.size(dates)
allMeas = np.empty(nDates,dtype=object)
for d_ind,thisDate in enumerate(dates):
    fList = glob(os.path.join('calData', thisDate,thisDate+'*.txt'))
    gamryFiles = [f for f in fList if re.search(thisDate+'-'+'[0-9][0-9][0-9][0-9]_', f)]

    nSweeps = np.size(gamryFiles)
    meas = np.empty(nSweeps, dtype=object)
    calStd = CalStdFit(interpMethod='cubic')

    for i, file in enumerate(gamryFiles):
        meas[i] = Solution(cmapName)
        meas[i].loadFile(file)

        if not np.isnan(meas[i].sigmaStd_Sm):
            meas[i].sigmaStdCalc_Sm = calStd(meas[i].T_K, lbl_uScm=meas[i].lbl_uScm)
        else:
            meas[i].sigmaStdCalc_Sm = 1e-8  # Default air conductivity
        meas[i].FitCircuit(circType=circType, initial_guess=initial_guess)
        meas[i].Kcell_pm = meas[i].sigmaStdCalc_Sm * meas[i].Rcalc_ohm

        if not PLOT_AIR and meas[i].comp == 'Air':
            gamryFiles.remove(file)
            cals.delete(cals[i])
            add = '_noAir'
    #     with open(file) as f:
    #         f.readline()  # Skip intro line
    #         meas[i].time = dtime.strptime(f.readline()[:-1], gamryDTstr)  # Measurement time
    #         meas[i].T_K = readFloat(f.readline())  # Temp
    #         meas[i].P_MPa = readFloat(f.readline())  # Pressure
    #         meas[i].descrip = f.readline()  # Text description
    #         meas[i].Vdrive_V = readFloat(f.readline())  # Driving voltage
    #         meas[i].fStart_Hz = readFloat(f.readline())  # Driving voltage
    #         meas[i].fStop_Hz = readFloat(f.readline())  # Driving voltage
    #         meas[i].nfSteps = int(readFloat(f.readline()))  # Number of f steps
    #
    #     if 'DIwater' in meas[i].descrip:
    #         meas[i].comp = 'Pure H2O'
    #         meas[i].sigmaStd_Sm = 0
    #         meas[i].legLabel = r'$\approx0$'
    #     elif 'Air' in meas[i].descrip:
    #         meas[i].comp = 'Air'
    #         meas[i].sigmaStd_Sm = np.nan
    #         meas[i].legLabel = 'Air'
    #     else:
    #         meas[i].comp = 'KCl'
    #         if '447uScm' in meas[i].descrip:
    #             # meas[i].sigmaStd_Sm = float(meas[i].descrip.split(':')[-1].split('uScm')[0]) / 1e4
    #             meas[i].sigmaStd_Sm = 447/1842 # 10kHz value from after-removal measurement on 20220923, phase at 10kHz is basically zero.
    #             meas[i].legLabel = '447uScm KCl std'
    #
    #         elif '2764uScm' in meas[i].descrip:
    #             meas[i].sigmaStd_Sm = 2764/614 # 10kHz value from after-removal measurement on 20220923, phase at 10kHz is basically zero.
    #             meas[i].legLabel = '2764uScm KCl std'
    #
    #         elif '8Sm' in meas[i].descrip:
    #             meas[i].sigmaStd_Sm = 80000 / 28  # 10kHz value from after-removal measurement on 20220923, phase at 10kHz is basically zero.
    #             meas[i].legLabel = '8Sm KCl std'

        # _, meas[i].f_Hz, Zabs_ohm, Phi_ohm = np.loadtxt(file, skiprows=10, unpack=True)
        # # Phase negated for convenience looking at plots in 1st quadrant
        # meas[i].Z_ohm = Zabs_ohm * np.exp(-1j*np.deg2rad(Phi_ohm))
        # # meas[i].sigma_Sm = np.real(np.interp(10e3,meas[i].f_Hz,meas[i].Z_ohm))*meas[i].sigmaStd_Sm)
        # # print(i)
        # meas[i].sigma_Sm = np.real(np.interp(10e3,meas[i].f_Hz, Zabs_ohm)) * meas[i].sigmaStd_Sm * 1e-4
        # circuit.fit(meas[i].f_Hz,meas[i].Z_ohm)
        # meas[i].R0 = circuit.parameters_[0]
        # meas[i].Z_fit = circuit.predict(meas[i].f_Hz)
        meas[i].sigma_Sm = meas[i].Kcell_pm/meas[i].Rcalc_ohm


    listStandards = np.array([cal.sigmaStd_Sm for cal in meas])
    uniqueListStandards = np.unique(listStandards)
    # standardsSort = np.argsort(uniqueListStandards)
    # for cal in meas:
    #     colorLabel = standardsSort[cal.sigmaStd_Sm == uniqueListStandards][0]
    #     cal.color = cmap(colorLabel / np.max(standardsSort))
    #
    #
    # # Sort list based on conductivity
    # iSort = np.argsort(listStandards)
    # meas = meas[iSort]

    allMeas[d_ind] = meas

PlotSigma(allMeas,figSize,outFigName,xtn)

# 1/Z_cell = 1/R + i*omega*C -- Pan et al. (2021): https://doi.org/10.1029/2021GL094020
