import numpy as np
import os
from glob import glob
import logging
from gamryTools import Resistor
from gamryPlots import PlotY, PlotZ, PlotZvsf, PlotPhasevsf, PlotZfit
from datetime import datetime as dtime

# Enable logging
logging.basicConfig(level=logging.DEBUG)

# Define the directory containing the data files
data_dir = '/Users/cal2019/Dropbox/ResearchVoyage/Models_code/HiPOZ/data/20230711/ConductivityData_100ohm_outsideInTube'

# Define the circuit type and other parameters
date = '20230711'
circType = 'RL'
cmapName = 'viridis'
outFigName = 'GamryCal'
figSize = (6, 5)
xtn = 'pdf'
PLOT_AIR = False
PAIR_DATA_W_FIT = True

# Import resistor data
data_files = glob(os.path.join(data_dir, '*.txt'))
resistors = []

for file in data_files:
    resistor = Resistor(cmapName=cmapName)
    with open(file, 'r') as f:
        lines = f.readlines()
        logging.debug(f"Processing file: {file}")
        header_index = None
        for i, line in enumerate(lines):
            logging.debug(f"Line: {line.strip()}")
            if line.strip() == 'Index #\tFrequency (Hz)\tImpedance (ohm)\tPhase (degrees)':
                header_index = i + 1
                break
        if header_index is None:
            continue  # Skip this file if the header line is not found

        for line in lines[header_index:]:
            columns = line.strip().split('\t')
            if len(columns) >= 4:
                timestamp = columns[0]
                impedance_real = columns[1]
                impedance_imag = columns[2]
                try:
                    color = 'red'  # Replace 'red' with a valid color value
                    resistor.addPoint(impedance_real, impedance_imag, timestamp, color)

                except ValueError:
                    logging.warning(f"Invalid data format. Skipping line: {line.strip()}")

        if resistor.numPoints > 0:
            color = 'red'  # Replace 'red' with a valid color value
            resistor.fitImpedanceVsf()  # Fit the impedance data and assign the fitting coefficients

            # Additional changes - Start
            resistor.lbl_uScm = 'label_value'  # Replace 'label_value' with the appropriate label in microsiemens per centimeter
            resistor.fitColor = 'red'  # Replace 'fit_color' with the appropriate color value for the fit
            resistor.Kcell_pm = 123.45  # Replace '123.45' with the appropriate value for Kcell_pm
            # Additional changes - End

            resistors.append(resistor)

logging.debug(f"len(resistors): {len(resistors)}")

print()
# Plotting
PlotZvsf(resistors, figSize, outFigName, xtn, Rticks=[100, 1000], LEG_PAIRS=PAIR_DATA_W_FIT)
# PlotPhasevsf(resistors, figSize, outFigName, xtn, LEG_PAIRS=PAIR_DATA_W_FIT)
# PlotZ(resistors, figSize, outFigName, xtn, Rticks=[100, 1000], LEG_PAIRS=PAIR_DATA_W_FIT)
# PlotY(resistors, figSize, outFigName, xtn, Rticks=[100, 1000], LEG_PAIRS=PAIR_DATA_W_FIT)
# PlotZfit(resistors, figSize, xtn, LEG_PAIRS=PAIR_DATA_W_FIT)
