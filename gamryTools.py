from collections.abc import Iterable
import numpy as np
import scipy.interpolate as sci
from datetime import datetime as dtime
from glob import glob
import logging
from matplotlib.cm import get_cmap
import matplotlib.colors as mc
import colorsys
from impedance.models.circuits import CustomCircuit
import schemdraw
import schemdraw.elements as elm

# Assign logger
log = logging.getLogger('HIPPOS')

gamryDTstr = r'%m/%d/%Y-%I:%M %p'
Lleads = 2

def readFloat(line):
    return float(line.split(':')[-1])

def LightenColor(color, lightnessMult=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    if type(color) == str:
        c = mc.cnames[color]
    else:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - lightnessMult * (1 - c[1]), c[2])

class Solution: 
    def __init__(self, cmapName):
        self.comp = None  # Solute composition
        self.w_ppt = None  # Solute mass concentration in g/kg
        self.P_MPa = None  # Chamber pressure of measurement in MPa
        self.T_K = None  # Temperature of measurement in K
        self.Vdrive_V = None  # Driving voltage in V
        self.fStart_Hz = None  # Start frequency in Hz
        self.fStop_Hz = None  # Stop frequency in Hz
        self.nfSteps = None  # Number of frequency steps
        self.time = None  # Measurement start time
        self.sigmaStd_Sm = None  # Known conductivity of standard solutions (as applicable)
        self.lbl_uScm = None  # Integer of the above in uS/cm
        self.descrip = None  # Text description (as applicable)
        self.legLabel = None  # Legend label
        self.color = None  # Color of lines
        self.cmap = get_cmap(cmapName)
        self.file = None  # File from which data has been loaded
        self.circFile = None  # File to print circuit diagram to
        self.xtn = 'pdf'

        # Outputs
        self.f_Hz = None  # Frequency values of Gamry sweep measurements in Hz
        self.Z_ohm = None  # Complex impedance values of Gamry sweep measurements in ohm
        self.Rcalc_ohm = None  # Fit from electrical impedance spectroscopy equivalent circuit modeling in ohm
        self.Runc_ohm = None  # Uncertainty in fit for Rcalc in ohm
        self.sigmaStdCalc_Sm = None  # Conductivity from interpolation of standard solution bottle label in S/m
        self.Kcell_pm = None  # Cell constant K in 1/m
        self.sigma_Sm = None  # DC electrical conductivity in S/m

    def loadFile(self, file):
        self.file = file
        with open(self.file) as f:
            f.readline()  # Skip intro line
            self.time = dtime.strptime(f.readline()[:-1], gamryDTstr)  # Measurement time
            self.T_K = readFloat(f.readline())  # Temp
            self.P_MPa = readFloat(f.readline())  # Pressure
            self.descrip = f.readline()  # Text description
            self.Vdrive_V = readFloat(f.readline())  # Driving voltage
            self.fStart_Hz = readFloat(f.readline())  # Driving voltage
            self.fStop_Hz = readFloat(f.readline())  # Driving voltage
            self.nfSteps = int(readFloat(f.readline()))  # Number of f steps

        if 'DIwater' in self.descrip:
            self.comp = 'Pure H2O'
            self.sigmaStd_Sm = 0
            self.legLabel = r'$\approx0$'
            self.lbl_uScm = 1
        elif 'Air' in self.descrip:
            self.comp = 'Air'
            self.sigmaStd_Sm = np.nan
            self.legLabel = 'Air'
            self.lbl_uScm = 1
        else:
            self.comp = 'KCl'
            if 'uScm' in self.descrip:
                self.sigmaStd_Sm = float(self.descrip.split(':')[-1].split('uScm')[0]) / 1e4
            elif 'mScm' in self.descrip:
                self.sigmaStd_Sm = float(self.descrip.split(':')[-1].split('mScm')[0]) / 10
            elif 'Sm' in self.descrip:
                self.sigmaStd_Sm = float(self.descrip.split(':')[-1].split('Sm')[0])
            else:
                self.sigmaStd_Sm = np.nan
            self.legLabel = f'{self.sigmaStd_Sm:.4f}'
            self.lbl_uScm = round(self.sigmaStd_Sm*1e4)

        self.color = self.cmap(np.log(self.lbl_uScm)/np.log(80000))
        self.fitColor = LightenColor(self.color, lightnessMult=0.4)

        _, self.f_Hz, Zabs_ohm, Phi_ohm = np.loadtxt(self.file, skiprows=10, unpack=True)
        self.Z_ohm = Zabs_ohm * np.exp(1j * np.deg2rad(Phi_ohm))

        return

    def FitCircuit(self, circType=None, initial_guess=None, BASIN_HOPPING=False, Kest_pm=None, PRINT=True, circFile=None):
        if Kest_pm is None:
            Kest_pm = 50
        if circType is None:
            circType = 'CPE'
        if circFile is None:
            self.circFile = f'{circType}circuit.{self.xtn}'
        else:
            self.circFile = circFile
        if circType == 'CPE':
            # Z_cell = R_0 + (R_0 + Z_CPE)/(1 + i*omega*C*(R_1 + Z_CPE)) -- Chin et al. (2018): https://doi.org/10.1063/1.5020076
            initial_guess = [Kest_pm/self.sigmaStdCalc_Sm, 8e-7, 0.85, 146.2e-12, 50]
            R0 = r'R_0'
            R1 = r'R_1'
            CPE1 = r'CPE_1'
            C1 = r'C_1'
            circStr = f'p({R1}-{CPE1},{C1})-{R0}'
            if PRINT:
                with schemdraw.Drawing(file=self.circFile, show=False) as circ:
                    circ.config(unit=Lleads)
                    circ += elm.Resistor().label(f'${R0}$').dot()
                    circ += (j1 := elm.Line().length(circ.unit/2).up())
                    circ += elm.Line().at(j1.start).length(circ.unit/2).down()
                    circ += elm.Resistor().right().label(f'${R1}$')
                    circ += elm.CPE().label(f'$Z_\mathrm{{{CPE1[:-2]}}}$')
                    circ += elm.Line().length(circ.unit/2).up().dot()
                    circ += (j2 := elm.Line().length(circ.unit/2).up())
                    circ += elm.Capacitor().endpoints(j1.end, j2.end).label(f'${C1}$').right()
                    circ += elm.Line().at(j2.start).length(circ.unit/2).right()
                log.info(f'Equivalent circuit diagram saved to file: {self.circFile}')
        elif circType == 'RC':
            # 1/Z_cell = 1/R + i*omega*C -- Pan et al. (2021): https://doi.org/10.1029/2021GL094020
            initial_guess = [Kest_pm/self.sigmaStdCalc_Sm, 146.2e-12]
            R1 = r'R_1'
            C1 = r'C_1'
            circStr = f'p({R1},{C1})'
            if PRINT:
                with schemdraw.Drawing(file=self.circFile, show=False) as circ:
                    circ.config(unit=Lleads)
                    circ += elm.Line().length(circ.unit/2).dot()
                    circ += (j1 := elm.Line().length(circ.unit/2).up())
                    circ += elm.Line().at(j1.start).length(circ.unit/2).down()
                    circ += elm.Resistor().right().label(f'${R1}$')
                    circ += elm.Line().length(circ.unit/2).up().dot()
                    circ += (j2 := elm.Line().length(circ.unit/2).up())
                    circ += elm.Capacitor().endpoints(j1.end, j2.end).label(f'${C1}$').right()
                    circ += elm.Line().at(j2.start).length(circ.unit/2).right()
                log.info(f'Equivalent circuit diagram saved to file: {self.circFile}')
        elif circType == 'RC-R':
            initial_guess = [Kest_pm/self.sigmaStdCalc_Sm, 146.2e-12, 50]
            R0 = r'R_0'
            R1 = r'R_1'
            C1 = r'C_1'
            circStr = f'p({R1},{C1})-{R0}'
            if PRINT:
                with schemdraw.Drawing(file=self.circFile, show=False) as circ:
                    circ.config(unit=Lleads)
                    circ += elm.Resistor().label(f'${R0}$').dot()
                    circ += (j1 := elm.Line().length(circ.unit/2).up())
                    circ += elm.Line().at(j1.start).length(circ.unit/2).down()
                    circ += elm.Resistor().right().label(f'${R1}$')
                    circ += elm.Line().length(circ.unit/2).up().dot()
                    circ += (j2 := elm.Line().length(circ.unit/2).up())
                    circ += elm.Capacitor().endpoints(j1.end, j2.end).label(f'${C1}$').right()
                    circ += elm.Line().at(j2.start).length(circ.unit/2).right()
                log.info(f'Equivalent circuit diagram saved to file: {self.circFile}')
        else:
            if initial_guess is None:
                raise ValueError(f'circuit type "{circType}" not recognized.')
            else:
                log.info(f'circuit type "{circType}" not recognized. Interpreting as circuit string.')
                circStr = circType
            
        log.debug(f'Fitting {circType} circuit to input file {self.file}')

        self.circuit = CustomCircuit(circStr, initial_guess=initial_guess)
        self.circuit.fit(self.f_Hz, self.Z_ohm, global_opt=BASIN_HOPPING)
        self.Zfit_ohm = self.circuit.predict(self.f_Hz)
        self.Rcalc_ohm = self.circuit.parameters_[0]
        self.Runc_ohm = self.circuit.conf_[0]
        log.debug(f'{self.circuit}' +
                  f'Fractional uncertainty in R: {self.Runc_ohm/self.Rcalc_ohm*100:.2f}%')

        return 


class expFit:
    def __init__(self, sigma0_Sm, lambd):
        self.sigma0_Sm = sigma0_Sm
        self.lambd = lambd

    def __call__(self, T_K):
        return self.sigma0_Sm * np.exp(self.lambd * T_K)


class CalStdFit:
    def __init__(self, interpMethod=None):
        self.calTemps_K = 273.15 + np.concatenate((np.array([5, 10]), np.arange(15, 31.5, 1)))
        self.calTable_Sm = {  # Conductivity measurements as a function of temperature as listed on the side of solution bottles
            23: 1e-4*np.array([15.32, 17.11, 18.32, 18.65, 18.97, 19.41, 19.85, 20.3, 20.92, 21.51, 22.1, 22.55, 23, 23.43, 23.9, 24.56, 25.2, 25.61, 26.4]),
            84: 1e-4*np.array([65, 67, 68, 70, 71, 73, 74, 76, 78, 79, 81, 82, 84, 86, 87, 89, 90, 92, 94]),
            447: 1e-4*np.array([278, 318, 361, 368, 376, 385, 394, 402, 411, 419, 430, 438, 447, 457, 467, 477, 487, 496, 505]),
            2070: 1e-4*np.array([1230, 1418, 1613, 1652, 1689, 1723, 1777, 1830, 1870, 1920, 1970, 2020, 2070, 2110, 2150, 2200, 2250, 2290, 2340]),
            2764: 1e-4*np.array([1752, 1998, 2244, 2296, 2346, 2400, 2448, 2502, 2554, 2606, 2658, 2712, 2764, 2816, 2872, 2922, 2952, 3030, 3082]),
            15000: 1e-4*np.array([9430, 10720, 12050, 12280, 12590, 12900, 13190, 13510, 13810, 14090, 14350, 14690, 15000, 15280, 15510, 15820, 16130, 16450, 16770]),
            80000: 1e-4*np.array([53800, 58900, 65700, 67200, 68700, 70100, 71500, 72900, 74300, 75800, 77200, 78600, 80000, 81600, 83100, 84700, 86200, 87700, 89400])
        }
        if interpMethod is None:
            fitsigma0_Sm = {23: 4.04e-6, 84: 6.882e-5, 447: 5.463e-5, 2070: 1.35e-4, 2764: 5.034e-4, 15000: 2.2572e-3, 80000: 2.3134e-2}  # Fit parameter for conductivity of each solution at 0 K
            fitlambd = {23: 0.0213, 84: 0.0161, 447: 0.0225, 2070: 0.0246, 2764: 0.0211, 15000: 0.0218, 80000: 0.0196}  # Fit parameter exponent for temp in K
            # Create function for evaluating fit conductivity
            self.sigmaCalc_Sm = {std: expFit(fitsigma0_Sm[std], fitlambd[std]) for std in self.calTable_Sm.keys()}

        else:
            # Interpolate label table to get conductivity
            self.sigmaCalc_Sm = {std: sci.interp1d(self.calTemps_K, self.calTable_Sm[std], kind=interpMethod, fill_value='extrapolate', bounds_error=False) for std in self.calTable_Sm.keys()}
            
    def __call__(self, T_K, lbl_uScm=None):
        if lbl_uScm is None:
            if isinstance(T_K, Iterable) and np.size(T_K) > 1:
                sigma_Sm = np.array([self.sigmaCalc_Sm[std](Ti_K) for std, Ti_K in zip(self.calTable_Sm.keys(), T_K)])
            else:
                sigma_Sm = np.array([self.sigmaCalc_Sm[std](T_K) for std in self.calTable_Sm.keys()])
        else:
            sigma_Sm = float(self.sigmaCalc_Sm[lbl_uScm](T_K))
        return sigma_Sm
