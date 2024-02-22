from impedance import preprocessing
from impedance.models.circuits import Randles, CustomCircuit
import matplotlib.pyplot as plt
from impedance.visualization import plot_nyquist
import pandas as pd
import numpy as np


# Define the file path for NaCl saturated
NaCl_file_path = '/Users/cal2019/Dropbox/ResearchVoyage/Models_code/HiPOZ/data/20240215/ConductivityData_RC_1kohm_10nF/20240215-123154_Temp295K.txt'

# Define the file path for KCl standard
KCl_file_path = '/Users/cal2019/Dropbox/ResearchVoyage/Models_code/HiPOZ/data/20231214/ConductivityData_KClStd_80ms_cm/20231214-110301_Temp293K.txt'


# Read the data into a pandas DataFrame, skip the header lines, and specify the delimiter
data = pd.read_csv(NaCl_file_path, skiprows=9, delimiter='\t')

# Display the first few rows of the DataFrame to verify the data was loaded correctly
print(data.head())

# Convert the pandas Series to numpy arrays
frequencies = np.array(data['Frequency (Hz)'])
impedance = np.array(data['Impedance (ohm)'])
phase = np.array(data['Phase (degrees)'])

# Combine impedance and phase into a single complex column
Z = impedance * np.exp(1j * np.radians(phase))

# RC Circuit
R1 = r'R_1'
C1 = r'C_1'
RC_circuit = f'p({R1},{C1})'
# Kest_pm = 80
# sigmaStdCalc_Sm = 8
# RC_initial_guess = [Kest_pm/sigmaStdCalc_Sm, 146.2e-12]
RC_initial_guess = [1e3, 10e-9]

# R0-R1C1 (RCR) circuit

# R0-R1C1-R2C2 circuit


# fit the CustomCircuit to the data
circuit = CustomCircuit(RC_circuit, initial_guess=RC_initial_guess)
circuit.fit(frequencies, Z)
print(circuit)