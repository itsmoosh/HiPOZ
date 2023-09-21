import os
import sys
import numpy as np
from glob import glob
import logging
from datetime import datetime as dtime
from gamryTools import ResistorData, CalStdFit, Solution
from gamryPlots import PlotZvsf, PlotPhasevsf, PlotZ, PlotY, PlotZfit

# Assign logger
log = logging.getLogger('HiPOZ')
stream = logging.StreamHandler(sys.stdout)
stream.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
log.setLevel(logging.DEBUG)
log.addHandler(stream)

R_nominal_ohm = 1000  # Nominal resistance value of the resistors

date = '20230720'
circType = 'RC'  # Options are 'CPE', 'RC', and 'RC-R'. If desired, a circuit string can be entered here instead.
initial_guess = None  # Required when circType is not in the above list. Ignored otherwise.
cmapName = 'viridis'
outFigName = 'ResistorCal'  # Output figure name
figSize = (6, 5)
xtn = 'pdf'
PAIR_DATA_W_FIT = True  # Whether to pair data and fit together for each solution in legend labels
LOW_SIG_CUTOFF = True  # Whether to cut off plotting of data at high frequencies (> 100 kHz) for low-conductivity solutions (< 2000 uS/cm), which get bypassed by transmission line effects
ADD_TICKS = False  # Whether to add ticks marking the fitted resistance/reactance values on plots

# Import data
add = None
gamryFiles = glob(os.path.join('data', date, '*', '*Temp*.txt'))

nSweeps = np.size(gamryFiles)
resistors = []  # Initialize as an empty list
calStd = CalStdFit(interpMethod='cubic')

# ... (Previous code remains unchanged)

for i, file in enumerate(gamryFiles):
    resistor = ResistorData(cmapName=cmapName)
    try:
        resistor.loadFile(file, gamryFiles)  # Pass gamryFiles as an argument
    except ValueError as ve:
        print(f"Error loading file {file}: {ve}")
        continue


    #print(f"Resistor {i + 1} - Impedance data length: {len(resistor.f_Hz)}")

    if len(resistor.f_Hz) == 0:
        print(f"Skipping Resistor {i + 1} - Empty impedance data")
        continue

    resistor.R_nominal_ohm = R_nominal_ohm

    # Debugging print statements
    print(f"Resistor {i + 1} - Frequency (Hz):", resistor.f_Hz)
    print(f"Resistor {i + 1} - Impedance (ohm):", resistor.Z_ohm)

    resistor.FitCircuit(circType=circType, initial_guess=initial_guess, PRINT=(i == 0))

    resistors.append(resistor)  # Append the resistor to the list

# Convert 'resistors' list to a NumPy array
resistors = np.array(resistors)

# Sort the array based on resistance values
listResistances = np.array([resistor.Rcalc_ohm for resistor in resistors])
filtered_resistances = [resistance for resistance in listResistances if resistance is not None]
iSort = np.argsort(filtered_resistances)
resistors = resistors[iSort]

# ... (Previous code remains unchanged)

if LOW_SIG_CUTOFF:
    for resistor in resistors:
        if resistor.Rcalc_ohm > 1e3:
            cutoff_index = resistor.f_Hz < 1e5  # Create the boolean index
            resistor.f_Hz = resistor.f_Hz[cutoff_index]
            resistor.Z_ohm = resistor.Z_ohm[cutoff_index]
            resistor.Zfit_ohm = resistor.Zfit_ohm[cutoff_index]

        print(f"Resistor after LOW_SIG_CUTOFF - Frequency (Hz):", resistor.f_Hz)
        print(f"Resistor after LOW_SIG_CUTOFF - Impedance (ohm):", resistor.Z_ohm)

for resistor in resistors:
    cutoff_index = resistor.f_Hz > 1e1  # Create the boolean index
    resistor.f_Hz = resistor.f_Hz[cutoff_index]
    resistor.Z_ohm = resistor.Z_ohm[cutoff_index]
    resistor.Zfit_ohm = resistor.Zfit_ohm[cutoff_index]

    print(f"Resistor after removing low-frequency points - Frequency (Hz):", resistor.f_Hz)
    print(f"Resistor after removing low-frequency points - Impedance (ohm):", resistor.Z_ohm)


if ADD_TICKS:
    Rtick_ohm = np.array([resistor.Rcalc_ohm for resistor in resistors])
    Ytick_mho = 1 / Rtick_ohm
else:
    Rtick_ohm = None
    Ytick_mho = None

PlotZvsf(resistors, figSize, outFigName, xtn, Rtick_ohm, LEG_PAIRS=PAIR_DATA_W_FIT)
PlotPhasevsf(resistors, figSize, outFigName, xtn, LEG_PAIRS=PAIR_DATA_W_FIT)
PlotZ(resistors, figSize, outFigName, xtn, Rtick_ohm, LEG_PAIRS=PAIR_DATA_W_FIT)
PlotY(resistors, figSize, outFigName, xtn, Ytick_mho, LEG_PAIRS=PAIR_DATA_W_FIT)
PlotZfit(resistors, figSize, xtn, LEG_PAIRS=PAIR_DATA_W_FIT)
