import os, sys, re
import numpy as np
from glob import glob
import logging
from datetime import datetime as dtime

from gamryTools import Solution, CalStdFit, TimeSeries
from gamryPlots import PlotY, PlotZ, PlotZvsf, PlotPhasevsf, PlotZfit, PlotTimeseries

from PyQt5.QtWidgets import QApplication
from DataSelector import DataSelector

from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# import sys
# from PyQt5.QtWidgets import QApplication
# from FileSelector import FileSelector
#
# app = QApplication(sys.argv)  # Create an application if one doesn't exist
# ex = FileSelector()
# ex.show()
# app.exec_()  # Start the event loop

# After the window is closed, you can access the selected files
# calibration_files = ex.calibration_files
# data_files = ex.data_files

# Assign logger
log = logging.getLogger('HiPOZ')
stream = logging.StreamHandler(sys.stdout)
stream.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
log.setLevel(logging.DEBUG)
log.addHandler(stream)

PAN_DATA = False
if PAN_DATA:
    dates = ['nan']
else:
    # dates = ['20220922','20220923','20220924','20221010','20221011','20230223']
    dates = ['20231214'] # values read at 5e-4S/m
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

def main():
    app = QApplication(sys.argv)  # Ensure QApplication is initialized
    # Import data
    add = None
    nDates = np.size(dates)
    allMeas = np.empty(nDates,dtype=object)
    for d_ind,thisDate in enumerate(dates):
        print(f'Processing {thisDate}/')
        pathHead = os.path.join(os.path.join('data', thisDate))
        gamryFiles = glob(os.path.join(pathHead, 'Conductivity*', '*.txt'))
        # gamryFiles = [f for f in fList if re.search(thisDate+'-'+'[0-9][0-9][0-9][0-9]_', f)]
        nSweeps = np.size(gamryFiles)
        meas = np.empty(nSweeps, dtype=object)
        calStd = CalStdFit(interpMethod='cubic')

        lf = len(gamryFiles)
        for i, file in enumerate(gamryFiles):
            meas[i] = Solution(cmapName=cmapName)
            print(f'measurement {i} of {lf}')
            meas[i].loadFile(file, PAN=PAN_DATA)

            if not np.isnan(meas[i].sigmaStd_Sm):
                meas[i].sigmaStdCalc_Sm = calStd(meas[i].T_K, lbl_uScm=meas[i].lbl_uScm)
            else:
                meas[i].sigmaStdCalc_Sm = 1e-8  # Default air conductivity
            meas[i].FitCircuit(circType=circType, initial_guess=initial_guess, PRINT=False,BASIN_HOPPING=False,MULTIPROC=True)

        allMeas[d_ind] = meas

    timeseries = TimeSeries(allMeas)
    timeseries.organizeData()
    ds = DataSelector(timeseries)
    ds.show()

    # Check if the main window is visible
    if not ds.isVisible():
        print("The window is not visible")

    # The geometry of the window can be printed to ensure it's within the visible area of your screen
    print("Window geometry:", ds.geometry())

    # Start the event loop
    return_code = app.exec_()
    print("Event loop exited with code:", return_code)
    sys.exit(return_code)

# PlotTimeseries(allMeas,figSize,outFigName,xtn)
if __name__ == "__main__":
    main()
# 1/Z_cell = 1/R + i*omega*C -- Pan et al. (2021): https://doi.org/10.1029/2021GL094020
