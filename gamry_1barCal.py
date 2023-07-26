import os, sys

import numpy as np
from glob import glob
import logging
from gamryTools import Solution, CalStdFit
from gamryPlots import PlotY, PlotZ, PlotZvsf, PlotPhasevsf, PlotZfit

# Assign logger
log = logging.getLogger('HiPOZ')
stream = logging.StreamHandler(sys.stdout)
stream.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
log.setLevel(logging.DEBUG)
log.addHandler(stream)

conc = '1e-3 molal'
Vrecipe_mL = 5000
Vbeaker_mL = 200
DeltaVbeaker_mL = 0.05*Vbeaker_mL
DeltamSolute_g = 0.001
Ttap_C = 25

Sol = Solution(comp='NaCl')
mSolute_g, Vwater_mL = Sol.Recipe(float(conc.split(' ')[0]), units=conc.split(' ')[-1], vol_mL=Vrecipe_mL, TH2O_C=Ttap_C)
mSolute_g = float(f'{mSolute_g:.3f}')  # Truncate to precision of scale as in DeltamSolute_g
Sol.wMeas_ppt, Sol.Deltaw_ppt, Sol.wMeas_molal, Sol.Deltaw_molal = Sol.CalcConc(mSolute_g, Vbeaker_mL, Vwater_mL,
                                           DeltamSolute_g=DeltamSolute_g, DeltaVbeaker_mL=DeltaVbeaker_mL, TH2O_C=Ttap_C)
print(
f"""Recipe for {Sol.w_ppt:.5f} ppt = {Sol.w_molal:.5f} molal {Sol.comp} (aq):
    mSalt (g): {mSolute_g:.3f}
    VH2O (mL): {Vwater_mL:.0f}
""")
# Get power of 10 in uncertainties for format string
fmtw_ppt = f'.{np.maximum(0,-np.floor(np.log10(Sol.Deltaw_ppt)) + 1):g}f'
fmtw_molal = f'.{np.maximum(0,-np.floor(np.log10(Sol.Deltaw_molal)) + 1):g}f'
print(
f"""Actual concentration and uncertainty:
    wMeas  (g/kg): {Sol.wMeas_ppt:{fmtw_ppt}} +/- {Sol.Deltaw_ppt:{fmtw_ppt}}
    wMeas (molal): {Sol.wMeas_molal:{fmtw_molal}} +/- {Sol.Deltaw_molal:{fmtw_molal}}
""")

date = '20221123'
circType = 'RC'  # Options are 'CPE', 'RC', and 'RC-R'. If desired, a circuit string can be entered here instead.
initial_guess = None  # Required when circType is not in the above list. Ignored otherwise.
cmapName = 'viridis'
outFigName = 'GamryCal'
figSize = (6,5)
xtn = 'pdf'
PLOT_AIR = False  # Whether to include dry run spectra in analysis
ADD_TICKS = False  # Whether to add ticks marking the fitted resistance/reactance values on plots
PAIR_DATA_W_FIT = True  # Whether to pair data and fit together for each solution in legend labels
LOW_SIG_CUTOFF = True  # Whether to cut off plotting of data at high frequencies (> 100 kHz) for low-conductivity solutions (< 2000 uS/cm), which get bypassed by transmission line effects

# Import data
add = None
gamryFiles = glob(os.path.join('data', date, '*', '*Temp*.txt'))

nSweeps = np.size(gamryFiles)
cals = np.empty(nSweeps, dtype=object)
calStd = CalStdFit(interpMethod='cubic')

for i, file in enumerate(gamryFiles):
    cals[i] = Solution(cmapName=cmapName)
    cals[i].loadFile(file)

    if not np.isnan(cals[i].sigmaStd_Sm):
        cals[i].sigmaStdCalc_Sm = calStd(cals[i].T_K, lbl_uScm=cals[i].lbl_uScm)
    else:
        cals[i].sigmaStdCalc_Sm = 1e-8  # Default air conductivity
    cals[i].FitCircuit(circType=circType, initial_guess=initial_guess, PRINT=(i == 0))
    cals[i].Kcell_pm = cals[i].sigmaStdCalc_Sm * cals[i].Rcalc_ohm

    if not PLOT_AIR and cals[i].comp == 'Air':
        gamryFiles.remove(file)
        cals.delete(cals[i])
        add = '_noAir'

# Sort list based on conductivity
listStandards = np.array([cal.sigmaStd_Sm for cal in cals])
iSort = np.argsort(listStandards)
cals = cals[iSort]

# Tcal_C = np.array([22.8, 22.3, 21.7, 21.5, 23.215, 21.8, 21.9])
# sigmaStdBottle_Sm = np.array([calStd(Tcal_C + 273.15)])
# sigmaStdlist_Sm = np.array([sol.sigmaStd_Sm for sol in cals])  # value on each bottle at 25 deg C, not using fits for temp dependence

if LOW_SIG_CUTOFF:
    for sol in cals:
        if sol.lbl_uScm < 2000:
            sol.Z_ohm = sol.Z_ohm[sol.f_Hz < 1e5]
            sol.Zfit_ohm = sol.Zfit_ohm[sol.f_Hz < 1e5]
            sol.f_Hz = sol.f_Hz[sol.f_Hz < 1e5]

for sol in cals:
    sol.Z_ohm = sol.Z_ohm[sol.f_Hz > 1e1]
    sol.Zfit_ohm = sol.Zfit_ohm[sol.f_Hz > 1e1]
    sol.f_Hz = sol.f_Hz[sol.f_Hz > 1e1]


if ADD_TICKS:
    Rtick_ohm = np.array([cal.Rcalc_ohm for cal in cals])
    Ytick_mho = 1/Rtick_ohm
else:
    Rtick_ohm = None
    Ytick_mho = None


PlotZvsf(cals, figSize, outFigName, xtn, Rtick_ohm, LEG_PAIRS=PAIR_DATA_W_FIT)
PlotPhasevsf(cals, figSize, outFigName, xtn, LEG_PAIRS=PAIR_DATA_W_FIT)
PlotZ(cals, figSize, outFigName, xtn, Rtick_ohm, LEG_PAIRS=PAIR_DATA_W_FIT)
PlotY(cals, figSize, outFigName, xtn, Ytick_mho, LEG_PAIRS=PAIR_DATA_W_FIT)
PlotZfit(cals, figSize, xtn, LEG_PAIRS=PAIR_DATA_W_FIT)
