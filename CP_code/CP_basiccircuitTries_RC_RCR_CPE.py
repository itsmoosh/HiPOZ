from impedance import preprocessing
from impedance.models.circuits import Randles, CustomCircuit
import matplotlib.pyplot as plt
from impedance.visualization import plot_nyquist
import pandas as pd
import numpy as np

# # Load data from the example EIS data
# frequencies, Z = preprocessing.readCSV('./data/2023091/ConductivityData_80mSKCltest/20230921-105603_Temp294K.txt')

# Define the file path for NaCl saturated
NaCl_file_path = 'data/20231214/ConductivityData_NaCl_4mol_kg/20231214-154247_Temp294K.txt'

# Define the file path for KCl standard
KCl_file_path = 'data/20231214/ConductivityData_KClStd_23us_cm/20231214-104516_Temp293K.txt'


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

# # keep only the impedance data in the first quandrant
# frequencies, Z = preprocessing.ignoreBelowX(frequencies, Z)

#-------------------------------------------------------
# example circuit given in impedance.py
DEFAULT_circuit = 'R0-p(R1,C1)-p(R2-Wo1,C2)'
DEFAULT_initial_guess = [.01, .01, 100, .01, .05, 100, 1]

#-----
# CPE circuit
# Z_cell = R_0 + (R_0 + Z_CPE)/(1 + i*omega*C*(R_1 + Z_CPE)) -- Chin et al. (2018): https://doi.org/10.1063/1.5020076
R0 = r'R_0'
R1 = r'R_1'
CPE1 = r'CPE_1'
C1 = r'C_1'
CPE_circuit = f'p({R1}-{CPE1},{C1})-{R0}'
# is the below for a non-standard solution using a standard?
# CPE_initial_guess = [Kest_pm / self.sigmaStdCalc_Sm, 8e-7, 0.85, 146.2e-12, 50]
CPE_initial_guess = [10, 8e-7, 0.85, 146.2e-12, 50]
#-----

# RC-R circuit
R0 = r'R_0'
R1 = r'R_1'
C1 = r'C_1'
#RC_R_circuit = f'p({R1},{C1})-{R0}'
RC_R_circuit = f'R_0-p(R_1,C_1)'
#RC_R_initial_guess = [10, 146.2e-12, 50]
RC_R_initial_guess = [.01, .005, .001, 200, .1]

#-----

# RC Circuit
R1 = r'R_1'
C1 = r'C_1'
RC_circuit = f'p({R1},{C1})'
Kest_pm = 80
sigmaStdCalc_Sm = 8
RC_initial_guess = [Kest_pm/sigmaStdCalc_Sm, 146.2e-12]
#-----

#fit the CustomCircuit to the data
circuit = CustomCircuit(DEFAULT_circuit, initial_guess=(DEFAULT_initial_guess))
circuit.fit(frequencies, Z)
print(circuit)

# Fit the Randles circuit to the data
randles_circuit = Randles(initial_guess=RC_R_initial_guess)
randles_circuit.fit(frequencies, Z)



#-------------------------------------------------------

# Plotting impedance vs. frequency
# fig, ax = plt.subplots()
# ax.plot(frequencies, np.abs(Z), 'o-', label='Data')  # Plotting magnitude of impedance
# ax.set_xscale('log')  # Set x-axis to logarithmic scale for frequency
# ax.set_yscale('log')  # Set y-axis to logarithmic scale for impedance magnitude
# ax.set_xlabel('Frequency (Hz)')
# ax.set_ylabel('Impedance Magnitude (ohm)')
# ax.set_title('Impedance vs. Frequency')
# ax.legend()

# Nyquist plot using matplotlib directly
# fig, ax = plt.subplots()
# ax.plot(np.real(impedance), -np.imag(impedance), ls='-', marker='o')
# ax.set_xlabel('Real(Z)')
# ax.set_ylabel('-Imag(Z)')
# plt.show()

# # Nyquist with impedance.py function
# fig, ax = plt.subplots(figsize=(5, 5))
# circuit.plot(f_data=frequencies, Z_data=Z, kind='nyquist', ax=ax)
# plt.show()

# ...

f_pred = np.logspace(5, -2)

randles_fit = randles_circuit.predict(f_pred)
customCircuit_fit = circuit.predict(f_pred)
print(randles_fit)
print(customCircuit_fit)

fig, ax = plt.subplots(figsize=(5, 5))
randles_circuit.plot(ax=ax)
circuit.plot(f_data=frequencies, Z_data=Z, kind='nyquist', ax=ax)



ax.legend(['Data', 'Randles', 'Custom Circuit'])
plt.show()
