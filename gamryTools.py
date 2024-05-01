from collections.abc import Iterable
import os
import numpy as np
import scipy.interpolate as sci
from scipy.optimize import basinhopping
from datetime import datetime as dtime
from glob import glob
import logging
from matplotlib.cm import get_cmap
import matplotlib.colors as mc
import colorsys
from impedance.models.circuits import CustomCircuit
import schemdraw
import schemdraw.elements as elm
from PlanetProfile.Thermodynamics.MgSO4.MgSO4Props import Ppt2molal, Molal2ppt
from PlanetProfile.Utilities.defineStructs import Constants
from seafreeze.seafreeze import seafreeze
import multiprocessing as mp
from tqdm import tqdm


# Assign logger
log = logging.getLogger('HiPOZ')

gamryDTstr = r'%m/%d/%Y-%I:%M %p'
PanDTstr = r'%m/%d/%Y %I:%M:%S %p'
PSItoMPa = 6.89476e-3
soluteOptions = ['DIwater', 'KCl', 'NaCl', 'MgSO4']  # Solutes available in PlanetProfile's Constants
wKCl_ppt = {23:    0.0116,
            84:    0.04038,
            447:   0.2256,
            2070:  1.045,
            2764:  1.382,
            15000: 8.759,
            80000: 52.168}  # Concentration of KCl standards as reported on the bottles
pptUnits = ['ppt', 'gkg', 'g/kg']  # Available options for indicating pptm units
molalUnits = ['molal', 'molkg', 'mol/kg','mol_kg']  # Available options for indicating molal units
SiemensUnits = ['us_cm','ms_cm']
Saturated = ['Saturated']
allUnits = pptUnits + molalUnits + SiemensUnits + Saturated

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


def GetSigmaFromDescrip(descrip):
    """
    Read Gamry file description for key terms to get conductivity standard info
    :param descrip (string): String description from Gamry file.
    :return sigmaStd_Sm (float): Conductivity standard label value in S/m.
    """
    if 'uScm' in descrip:
        sigmaStd_Sm = float(descrip.split(':')[-1].split('uScm')[0]) / 1e4
    elif 'mScm' in descrip:
        sigmaStd_Sm = float(descrip.split(':')[-1].split('mScm')[0]) / 10
    elif 'Sm' in descrip:
        sigmaStd_Sm = float(descrip.split(':')[-1].split('Sm')[0])
    else:
        sigmaStd_Sm = np.nan
    return sigmaStd_Sm


def GetwFromDescrip(descrip, lbl_uScm=None):
    """
    Read Gamry file description for key terms to get concentration.
    :param descrip (string): String description from Gamry file.
    :param lbl_uScm (int): Conductivity label for KCl standards
    :return comp (string): Composition of the solution.
    :return w_ppt (float): Concentration value in ppt.
    :return w_molal (float): Concentration value in molal.
    """
    soluteCompare = [comp in descrip for comp in soluteOptions]  # Boolean list of comparisons to available options
    if np.sum(soluteCompare) > 1:
        raise ValueError(f'Multiple solutes found in descrip: "{descrip}".')

    if lbl_uScm is not None and not np.isnan(lbl_uScm):
        comp = 'KCl'
        w_ppt = wKCl_ppt[lbl_uScm]
        w_molal = Ppt2molal(w_ppt, m_gmol[comp])
    elif any(soluteCompare):
        if 'DIwater' in descrip:
            comp = 'PureH2O'
            w_ppt = 0
            w_molal = np.nan
        else:
            comp = list(soluteOptions)[np.where(soluteCompare)[0][0]]
            unitCompare = [units in descrip for units in allUnits]
            if not any(unitCompare):
                # raise ValueError(f'Units not found in descrip "{descrip}".')
                log.info(f'No units found in descrip: {descrip}')
            # unit = list(allUnits)[np.where(unitCompare)[0][0]]
            try:
                unit = list(allUnits)[np.where(unitCompare)[0][0]]
                # Strip whitespace characters like tabs and newlines
                cleaned_string = descrip.strip()

                # Split the string on underscore
                parts = cleaned_string.split('_')

                # Initialize a variable to hold the extracted number
                number_as_float = None

                # Iterate over the parts to extract digits
                for part in parts:
                    # Extract digits from the part
                    if comp not in part:
                        # Extract digits from the part
                        numeric_part = ''.join([char for char in part if char.isdigit()])

                        if numeric_part:  # Check if we found any digits
                            number_as_float = float(numeric_part)
                            break
                w = number_as_float
            except ValueError:
                raise ValueError('Concentration and units must appear at the beginning of descrip.')
            if unit in pptUnits:
                w_ppt = w + 0
                w_molal = Ppt2molal(w_ppt, Constants.m_gmol[comp])
            elif unit in Saturated:
                w_molal = 1000
                w_ppt = 1000
            else:
                w_molal = w + 0
                w_ppt = Molal2ppt(w_molal, Constants.m_gmol[comp])
    else:
        log.warning(f'Unable to interpret descrip "{descrip}".')
        comp, w_ppt, w_molal = (None, None, None)
    return comp, w_ppt, w_molal


class Solution: 
    def __init__(self, comp=None, cmapName='viridis'):
        self.comp = comp  # Solute composition
        self.circStr = None
        self.w_ppt = None  # Solute mass concentration in g solute/kg solution
        self.w_molal = None  # Solute mass concentration in mol solute/kg solvent
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
        self.ramp = None  # Pan data warming or cooling ramp direction
        self.wMeas_ppt = None  # Final solute mass concentration in g solute/kg solution based on measured-out quantities
        self.wMeas_molal = None  # Final solute mass concentration in mol solute/kg solvent based on measured-out quantities
        self.Deltaw_ppt = None  # Uncertainty in solute mass concentration ppt value
        self.Deltaw_molal = None  # Uncertainty in solute mass concentration molal value

        # Outputs
        self.f_Hz = None  # Frequency values of Gamry sweep measurements in Hz
        self.Z_ohm = None  # Complex impedance values of Gamry sweep measurements in ohm
        self.Rcalc_ohm = None  # Fit from electrical impedance spectroscopy equivalent circuit modeling in ohm
        self.Runc_ohm = None  # Uncertainty in fit for Rcalc in ohm
        self.sigmaStdCalc_Sm = None  # Conductivity from interpolation of standard solution bottle label in S/m
        self.Kcell_pm = None  # Cell constant K in 1/m
        self.sigma_Sm = None  # DC electrical conductivity in S/m

    def loadFile(self, file, PAN=False):
        self.file = file
        with open(self.file) as f:
            if PAN:
                f.readline()
                f.readline()
                self.time = dtime.strptime(f.readline()[:-2].split(':  ')[-1], PanDTstr)
                f.readline()
                f.readline()
                Pstring, Tstring = f.readline()[:-1].split(os.path.sep)[-1].split('-')
                self.P_MPa = float(Pstring.split('psi')[0]) * PSItoMPa
                if Tstring[0] == 'U':
                    self.ramp = 'warming'
                    Tstring = Tstring[1:]
                else:
                    self.ramp = 'cooling'
                if Tstring[0] == 'n':
                    negC = -1
                    Tstring = Tstring[1:]
                else:
                    negC = 1
                self.T_K = float(Tstring.split('deg')[0]) * negC + 273.15

                # Pan data files do not contain the following information
                self.descrip = 'Pan ' + self.ramp  # Text description
                self.Vdrive_V = np.nan  # Driving voltage
                self.fStart_Hz = np.nan  # Spectrum begin frequency
                self.fStop_Hz = np.nan  # Spectrum end frequency
                self.nfSteps = np.size(self.f_Hz)  # Number of f steps
            else:
                f.readline()  # Skip intro line
                self.time = dtime.strptime(f.readline()[:-1], gamryDTstr)  # Measurement time
                self.T_K = readFloat(f.readline())  # Temp
                self.P_MPa = readFloat(f.readline())  # Pressure
                self.descrip = f.readline()  # Text description
                self.Vdrive_V = readFloat(f.readline())  # Driving voltage
                self.fStart_Hz = readFloat(f.readline())  # Spectrum begin frequency
                self.fStop_Hz = readFloat(f.readline())  # Spectrum end frequency
                self.nfSteps = int(readFloat(f.readline()))  # Number of f steps

        if any([comp in self.descrip for comp in soluteOptions]):
            self.comp, self.w_ppt, self.w_molal = GetwFromDescrip(self.descrip)
            if 'DIwater' in self.descrip:
                self.sigmaStd_Sm = 0
                self.legLabel = r'$\approx0$'
                self.lbl_uScm = 1
            else:
                self.sigmaStd_Sm = np.nan
                self.legLabel = f'{self.w_ppt:.1f}\,ppt\,\\ce{{{self.comp}}}'
                self.lbl_uScm = np.minimum(1, round(self.w_ppt))
        elif 'Air' in self.descrip:
            self.comp = 'Air'
            self.sigmaStd_Sm = np.nan
            self.legLabel = 'Air'
            self.lbl_uScm = 1
        elif 'Pan' in self.descrip:
            self.comp = 'NaCl'
            self.sigmaStd_Sm = np.nan
            self.legLabel = self.descrip + f'{self.T_K:.0f}'
            self.lbl_uScm = np.nan
        elif 'special' in self.descrip.lower():
            self.comp = self.descrip.split(':')[-1].split('pecial_')[-1]
            self.sigmaStd_Sm = np.nan
            self.legLabel = self.comp
            self.lbl_uScm = np.nan
        else:
            self.sigmaStd_Sm = GetSigmaFromDescrip(self.descrip)
            self.legLabel = f'{self.sigmaStd_Sm:.4f}'
            self.lbl_uScm = np.round(self.sigmaStd_Sm*1e4)
            self.comp, self.w_ppt, self.w_molal = GetwFromDescrip(self.descrip, lbl_uScm=self.lbl_uScm)


        self.color = self.cmap(np.log(self.lbl_uScm)/np.log(80000))
        self.fitColor = LightenColor(self.color, lightnessMult=0.4)

        if PAN:
            _, self.f_Hz, Zprime_ohm, ZdblPrime_ohm, _, _ = np.loadtxt(self.file, skiprows=7, unpack=True)
            self.Z_ohm = Zprime_ohm + 1j*ZdblPrime_ohm
        else:
            _, self.f_Hz, Zabs_ohm, Phi_ohm = np.loadtxt(self.file, skiprows=10, unpack=True)
            self.Z_ohm = Zabs_ohm * np.exp(1j * np.deg2rad(Phi_ohm))
        return

    def FitCircuit(self, circType=None, initial_guess=None, BASIN_HOPPING=False,MULTIPROC=False, Kest_pm=None, PRINT=True, circFile=None):
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
            # initial_guess = [Kest_pm/self.sigmaStdCalc_Sm, 8e-7, 0.85, 146.2e-12, 50]
            initial_guess = [np.real(self.Z_ohm[0]),0.02, 1e-5, 0.9, 1e-6]
            R0 = r'R_0'
            R1 = r'R_1'
            CPE1 = r'CPE_1'
            C1 = r'C_1'
            self.circStr = f'{R0}-p({R1}-{CPE1},{C1})'
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
            self.circStr = f'p({R1},{C1})'
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

        # Perform multiple fits with random initial guesses
        if BASIN_HOPPING:
            n_trials = 5
            results = []
            bounds = [(0.01, 5), (1e-4, 1), (1e-6, 1), (0.8, 1), (1e-7, 1e-5)]  # CPE
            if MULTIPROC:
                log.info(f"Starting trials in multiprocessing mode with {mp.cpu_count()} cores")
                with mp.Pool(processes=mp.cpu_count()) as pool:
                    args = [(bounds) for _ in range(n_trials)]
                    for result in tqdm(pool.imap(self.optimize_circuit, args), total=n_trials):
                        results.append(result)
            else:
                for n in tqdm(range(n_trials), desc="Optimizing"):
                    log.info(f"Starting trial {n+1}")
                    result = self.optimize_circuit(bounds)
                    results.append(result)
            # Find the best result based on lowest cost
            best_result = min(results, key=lambda x: x.fun)

        # Modify the optimizer settings
        minimizer_kwargs = {
            "method": 'L-BFGS-B',
            "bounds": [(1e-16, 5), (1e-16, 1), (1e-16, 1), (1e-16, 1), (1e-16, 1e-4)],
            "options": {'disp': True}  # Display convergence messages
        }
        self.circuit = CustomCircuit(self.circStr, initial_guess=initial_guess)
        self.circuit.fit(self.f_Hz, self.Z_ohm)
        self.Zfit_ohm = self.circuit.predict(self.f_Hz)
        self.Rcalc_ohm = self.circuit.parameters_[0]
        self.Runc_ohm = self.circuit.conf_[0]
        log.debug(f'{self.circuit}' +
                  f'Fractional uncertainty in R: {self.Runc_ohm/self.Rcalc_ohm*100:.2f}%')
        return

    def optimize_circuit(self,bounds):
        try:
            initial_guess = [np.random.uniform(low, high) for low, high in bounds]
            result = basinhopping(self.create_and_fit_circuit, initial_guess, niter=100, stepsize=0.01,
                                  minimizer_kwargs={"method": "BFGS"})
            return result
        except Exception as e:
            print("Error in optimize_circuit:", str(e))
            return None  # Explicitly return None if there is an error

    # Setup basin hopping
    def create_and_fit_circuit(self,bounds):
        try:
            for i, (low, high) in enumerate(bounds):
                params[i] = np.clip(params[i], low, high)

            self.circuit = CustomCircuit(self.circStr, initial_guess=params)
            self.circuit.fit(self.f_Hz, self.Z_ohm)

            # Calculate the residuals (difference between observed and modeled impedance)
            Z_fit = self.circuit.predict(self.f_Hz)
            residuals = self.Z_ohm - Z_fit
            # Calculate sum of squared errors (SSE)
            sse = np.sum(np.abs(residuals)**2)
            return sse  # Return SSE as the cost
        except Exception as e:
            print("Error during fitting:", str(e))
            return float('inf')

    def Recipe(self, w, units='ppt', vol_mL=500, TH2O_C=25):
        """
        Determine the recipe for mixing the desired concentration of solution.
        :param w: Concentration in units of ppt (g solute per kg of solution) or molal (moles solute per kg of solvent).
        :param units: One of the available units defined in pptUnits or molalUnits.
        :param vol_mL: Volume of solution to mix
        :param TH2O_C: Temperature of DI water to use in Celsius
        :return mSolute_g: Mass of solute to mix in grams
        :return volOut_mL: Volume of DI water to mix in milliliters
        """
        failMsg = 'unable to continue with Recipe.'

        # Make sure we can do calculations for the set composition
        if self.comp is None:
            log.warning(f'Composition is unset, {failMsg}')
            return
        elif not self.comp in Constants.m_gmol.keys():
            log.warning(f'Composition "{self.comp}" not yet covered in PlanetProfile Constants.m_gmol dict, {failMsg}')
            return

        # Make sure all the fields we need are set
        if units not in allUnits:
            log.warning(f'Units "{units}" not recognized, {failMsg}')
            return

        if self.w_ppt is None and self.w_molal is None:
            if units in pptUnits:
                self.w_ppt = w
            elif units in molalUnits:
                self.w_molal = w
            else:
                log.warning(f'Concentration is unset (Solution.w_ppt or Solution.w_molal), {failMsg}')
                return

        # Set either missing w field
        if self.w_molal is None:
            self.w_molal = Ppt2molal(self.w_ppt, Constants.m_gmol[self.comp])
            log.debug('Added w_molal to Solution.')
        elif self.w_ppt is None:
            self.w_ppt = Molal2ppt(self.w_molal, Constants.m_gmol[self.comp])
            log.debug('Added w_ppt to Solution.')

        rhoH2O_kgm3 = seafreeze(np.array(tuple((Constants.P0, TH2O_C+Constants.T0))), 'water1').rho[0,0]
        mH2O_kg = rhoH2O_kgm3/1e3 * vol_mL/1e3
        mSolute_g = self.w_molal * mH2O_kg * Constants.m_gmol[self.comp]
        volOut_mL = vol_mL + 0

        return mSolute_g, volOut_mL

    def CalcConc(self, mSolute_g, Vbeaker_mL, Vtotal_mL, DeltamSolute_g=0.001, DeltaVbeaker_mL=10, TH2O_C=25):
        """
        Calculate actual concentration of solute based on measured-out amounts, with uncertainty.
        :param mSolute_g: Mass of solute in grams.
        :param Vbeaker_mL: Volume of beaker used to measure out the desired amount of water in milliliters.
        :param Vtotal_mL: Total measured-out amount of water in milliliters.
        :param DeltamSolute_g: Uncertainty in solute mass in g. Typically 0.001 g based on available weighing scale.
        :param DeltaVbeaker_mL: Uncertainty in measured volume of beaker used to measure out desired amount of water, in milliliters.
        :param TH2O_C: Temperature of DI water in Celsius
        :return wMeas_ppt: Measured mass concentration in g solute/kg solution.
        :return wMeas_molal: Measured mass concentration in mol solute/kg solvent.
        :return Deltaw_ppt: Uncertainty in measured mass concentration in g solute/kg solution.
        :return Deltaw_molal: Uncertainty in measured mass concentration in mol solute/kg solvent.
        """
        # Get water density and approximate uncertainty based on assuming +/- 1 K uncertainty in tap water temp
        rhoH2O_kgm3 = seafreeze(np.array(tuple((Constants.P0, TH2O_C + Constants.T0))), 'water1').rho[0,0]
        DeltarhoH2O_kgm3 = abs(seafreeze(np.array(tuple((Constants.P0, TH2O_C+Constants.T0-1))), 'water1').rho[0,0]
                               - seafreeze(np.array(tuple((Constants.P0, TH2O_C+Constants.T0+1))), 'water1').rho[0,0])/2
        # Calculate fractional uncertainties
        eps_mSolute = DeltamSolute_g/mSolute_g
        eps_Vtotal = (DeltaVbeaker_mL * np.sqrt(Vtotal_mL/Vbeaker_mL)) / Vtotal_mL
        eps_rhoH2O = DeltarhoH2O_kgm3/rhoH2O_kgm3

        mH2O_kg = rhoH2O_kgm3 * Vtotal_mL/1e6
        eps_mH2O = np.sqrt(eps_rhoH2O**2 + eps_Vtotal**2)
        molSolute = mSolute_g / Constants.m_gmol[self.comp]
        eps_molSolute = eps_mSolute  # Assume much higher precision in molecular weight
        wMeas_molal = molSolute / mH2O_kg
        Deltaw_molal = wMeas_molal * np.sqrt(eps_molSolute**2 + eps_mH2O**2)
        Deltaw_molal = float(f'{Deltaw_molal:.2e}')  # Truncate to 2 sig figs in precision
        mSolute_kg = mSolute_g/1e3
        mTotal_kg = mH2O_kg + mSolute_kg
        eps_mTotal = np.sqrt((mH2O_kg*eps_mH2O)**2 + (mSolute_kg*eps_mSolute)**2) / mTotal_kg
        wMeas_ppt = mSolute_g / mTotal_kg
        Deltaw_ppt = wMeas_ppt * np.sqrt(eps_mSolute**2 + eps_mTotal**2)
        Deltaw_ppt = float(f'{Deltaw_ppt:.2e}')  # Truncate to 2 sig figs in precision

        return wMeas_ppt, Deltaw_ppt, wMeas_molal, Deltaw_molal
class ResistorData:
    def __init__(self, comp=None, cmapName='viridis'):

        self.P_MPa = None  # Chamber pressure of measurement in MPa
        self.T_K = None  # Temperature of measurement in K
        self.Vdrive_V = None  # Driving voltage in V
        self.fStart_Hz = None  # Start frequency in Hz
        self.fStop_Hz = None  # Stop frequency in Hz
        self.nfSteps = None  # Number of frequency steps
        self.time = None  # Measurement start time
        self.descrip = None  # Text description (as applicable)
        self.legLabel = None  # Legend label
        self.color = None  # Color of lines
        self.cmap = get_cmap(cmapName)
        self.file = None  # File from which data has been loaded
        self.circFile = None  # File to print circuit diagram to
        self.xtn = 'pdf'
        # self.frequency = None
        # self.impedance = None
        # self.phase = None

    # Outputs
        # Initialize as empty lists
        self.f_Hz = []  # Frequency values of Gamry sweep measurements in Hz
        self.Z_ohm = []  # Complex impedance values of Gamry sweep measurements in ohm
        self.phase = []  # Phase values of Gamry sweep measurements in degrees
        # self.f_Hz = None  # Frequency values of Gamry sweep measurements in Hz
        # self.Z_ohm = None  # Complex impedance values of Gamry sweep measurements in ohm
        self.Rcalc_ohm = None  # Fit from electrical impedance spectroscopy equivalent circuit modeling in ohm
        self.Runc_ohm = None  # Uncertainty in fit for Rcalc in ohm

    def convert_to_numpy(self):
        self.f_Hz = np.array(self.f_Hz)
        magnitude = np.array(self.Z_ohm)
        phase_rad = np.radians(np.array(self.phase))
        self.Z_ohm = magnitude * np.exp(1j * phase_rad)

    def loadFile(self, filename, all_files):
        resistor_data = self  # Create an instance of the class

        # Read data from the file
        with open(filename, 'r') as file:
            lines = file.readlines()

            # Extract relevant data starting from line 11 (skip the first 10 lines with non-data content)
            for line in lines[10:]:
                # Split the line by spaces and extract data fields
                fields = line.strip().split()
                frequency = float(fields[1])
                impedance = float(fields[2])
                phase = float(fields[3])

                # Add data to the impedance_data list
                resistor_data.f_Hz.append(frequency)
                resistor_data.Z_ohm.append(impedance)
                resistor_data.phase.append(phase)

                # Debug print statements
                print("Frequency (Hz):", frequency)
                print("Impedance (ohm):", impedance)
                print("Phase (degrees):", phase)

        # Convert the lists to NumPy arrays
        resistor_data.f_Hz = np.array(resistor_data.f_Hz)
        resistor_data.Z_ohm = np.array(resistor_data.Z_ohm)
        resistor_data.phase = np.array(resistor_data.phase)

        # Check if impedance data is not empty
        if len(resistor_data.f_Hz) == 0 or len(resistor_data.Z_ohm) == 0:
            raise ValueError('Impedance data is empty. Please load valid data.')

        # Calculate the color index based on some property (here, it's the length of f_Hz array)
        color_index = len(resistor_data.f_Hz)

        # Set the fit color using the cmap defined in the class constructor
        resistor_data.fitColor = resistor_data.cmap(color_index / len(all_files))

        return resistor_data
    def FitCircuit(self, circType=None, initial_guess=None, BASIN_HOPPING=False, Kest_pm=None, PRINT=True,
                   circFile=None):

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
#            initial_guess = [Kest_pm / self.sigmaStdCalc_Sm, 8e-7, 0.85, 146.2e-12, 50]
            initial_guess = [0.01,  8e-7, 0.85, 146.2e-12, 50] #CP for resistor testing
            R0 = r'R_0'
            R1 = r'R_1'
            CPE1 = r'CPE_1'
            C1 = r'C_1'
            circStr = f'p({R1}-{CPE1},{C1})-{R0}'
            if PRINT:
                with schemdraw.Drawing(file=self.circFile, show=False) as circ:
                    circ.config(unit=Lleads)
                    circ += elm.Resistor().label(f'${R0}$').dot()
                    circ += (j1 := elm.Line().length(circ.unit / 2).up())
                    circ += elm.Line().at(j1.start).length(circ.unit / 2).down()
                    circ += elm.Resistor().right().label(f'${R1}$')
                    circ += elm.CPE().label(f'$Z_\mathrm{{{CPE1[:-2]}}}$')
                    circ += elm.Line().length(circ.unit / 2).up().dot()
                    circ += (j2 := elm.Line().length(circ.unit / 2).up())
                    circ += elm.Capacitor().endpoints(j1.end, j2.end).label(f'${C1}$').right()
                    circ += elm.Line().at(j2.start).length(circ.unit / 2).right()
                log.info(f'Equivalent circuit diagram saved to file: {self.circFile}')
        elif circType == 'RC':
            # 1/Z_cell = 1/R + i*omega*C -- Pan et al. (2021): https://doi.org/10.1029/2021GL094020
            # initial_guess = [Kest_pm / self.sigmaStdCalc_Sm, 146.2e-12]
            initial_guess = [1, 146.2e-12]
            R1 = r'R_1'
            C1 = r'C_1'
            circStr = f'p({R1},{C1})'
            if PRINT:
                with schemdraw.Drawing(file=self.circFile, show=False) as circ:
                    circ.config(unit=Lleads)
                    circ += elm.Line().length(circ.unit / 2).dot()
                    circ += (j1 := elm.Line().length(circ.unit / 2).up())
                    circ += elm.Line().at(j1.start).length(circ.unit / 2).down()
                    circ += elm.Resistor().right().label(f'${R1}$')
                    circ += elm.Line().length(circ.unit / 2).up().dot()
                    circ += (j2 := elm.Line().length(circ.unit / 2).up())
                    circ += elm.Capacitor().endpoints(j1.end, j2.end).label(f'${C1}$').right()
                    circ += elm.Line().at(j2.start).length(circ.unit / 2).right()
                log.info(f'Equivalent circuit diagram saved to file: {self.circFile}')
        elif circType == 'RC-R':
#            initial_guess = [Kest_pm / self.sigmaStdCalc_Sm, 146.2e-12, 50]
            initial_guess = [100, 146.2e-12, 50] #cp resistors
            R0 = r'R_0'
            R1 = r'R_1'
            C1 = r'C_1'
            circStr = f'p({R1},{C1})-{R0}'
            if PRINT:
                with schemdraw.Drawing(file=self.circFile, show=False) as circ:
                    circ.config(unit=Lleads)
                    circ += elm.Resistor().label(f'${R0}$').dot()
                    circ += (j1 := elm.Line().length(circ.unit / 2).up())
                    circ += elm.Line().at(j1.start).length(circ.unit / 2).down()
                    circ += elm.Resistor().right().label(f'${R1}$')
                    circ += elm.Line().length(circ.unit / 2).up().dot()
                    circ += (j2 := elm.Line().length(circ.unit / 2).up())
                    circ += elm.Capacitor().endpoints(j1.end, j2.end).label(f'${C1}$').right()
                    circ += elm.Line().at(j2.start).length(circ.unit / 2).right()
                log.info(f'Equivalent circuit diagram saved to file: {self.circFile}')
        else:
            if initial_guess is None:
                raise ValueError(f'circuit type "{circType}" not recognized.')
            else:
                log.info(f'circuit type "{circType}" not recognized. Interpreting as circuit string.')
                circStr = circType

        log.debug(f'Fitting {circType} circuit to input file {self.file}')
        log.debug(f'Frequency (Hz) for fitting: {self.f_Hz}') #cp
        log.debug(f'Impedance (ohm) for fitting: {self.Z_ohm}') #cp

        # Convert impedance data to NumPy arrays if not done already
        self.convert_to_numpy()

        # Check if impedance data is not empty (CP)
        if len(self.f_Hz) == 0 or len(self.Z_ohm) == 0:
            raise ValueError('Impedance data is empty. Please load valid data before fitting.')

        self.circuit = CustomCircuit(circStr, initial_guess=initial_guess)
        self.circuit.fit(self.f_Hz, self.Z_ohm, global_opt=BASIN_HOPPING)
        self.Zfit_ohm = self.circuit.predict(self.f_Hz)
        self.Rcalc_ohm = self.circuit.parameters_[0]
        self.Runc_ohm = self.circuit.conf_[0]
        log.debug(f'{self.circuit}' +
                  f'Fractional uncertainty in R: {self.Runc_ohm / self.Rcalc_ohm * 100:.2f}%')

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


