from impedance import preprocessing
from impedance.models.circuits import Randles, CustomCircuit
import matplotlib.pyplot as plt
from impedance.visualization import plot_nyquist
import pandas as pd
import numpy as np


# # Load data from the example EIS data
# frequencies, Z = preprocessing.readCSV('./data/2023091/ConductivityData_80mSKCltest/20230921-105603_Temp294K.txt')

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

# Lima 2017 Circuit
# Define your circuit parameters
# Parameters
R0 = r'R_0'
R1 = r'R_1'
C1 = r'C_1'
C2 = r'C_2'
CPE1 = r'CPE_1'
L1 = r'L_1'
R2 = r'R_2'

# Define the custom circuit as a string
# Lima_2017_circuit = f'p({C1},{R1})-p({C2},{CPE1})-{R2}-{L1}'
Lima_2017_circuit = f'p(C_1, R_1)-p(C_2,CPE_1)-R_2-L_1'
Lima_initial_guess = [3e-6, 2e3, 3e-6, 4e-6, 0.75, 9, 4e-6]

# Randles Circuit (Rahmawati 2022 model a)
randles = Randles(initial_guess=[.01, .005, .001, 200, .1])
# randles.fit(frequencies,Z)
# print(randles)

# Rahmawati 2022 model b circuit
CPE = r'CPE'
R1 = r'R_1'
Rahmawati_model_b = f'R_1-CPE'
Rahmawati_model_b_initial_guess = [1, 4e-6, 9e-1]

#  Rahmawati 2022 model c circuit
Rahmawati_model_c = f'R_1-p(C_1, R_2)-CPE'
Rahmawati_model_c_initial_guess = [1, 1, 1, 1, 1]

# Marshall and Catherine Circuits (addendum to Rahmawati 2022 model b)

R1 = r'R_1'
CPE1 = r'CPE_1'
R2 = r'R_2'
C1 = r'C_1'
MarshallCircuitIdea = f'p({R1}-{CPE1},{R2}-{C1})'

R1 = r'R_1'
CPE1 = r'CPE_1'
C1 = r'C_1'
MarshallCircuitIdea2 = f'{R1}-p({CPE1},{C1})'

MarshallInitialGuess = [10, 5e-6, 8e-1, 1e-6]

R1 = r'R_1'
CPE1 = r'CPE_1'
C1 = r'C_1'
CatherineIdea = f'p({R1},{C1})-{CPE1}'
CatherineInitialGuess = [10, 1, 1, 1]

# fit the CustomCircuit to the data
circuit = CustomCircuit(RC_circuit, initial_guess=RC_initial_guess)
circuit.fit(frequencies, Z)
print(circuit)