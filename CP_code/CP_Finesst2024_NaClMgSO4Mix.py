from impedance import preprocessing
from impedance.models.circuits import Randles, CustomCircuit
import matplotlib.pyplot as plt
from impedance.visualization import plot_nyquist
import pandas as pd
import numpy as np


# # Load data from the example EIS data
# frequencies, Z = preprocessing.readCSV('./data/2023091/ConductivityData_80mSKCltest/20230921-105603_Temp294K.txt')

# Define the file path for NaCl saturated
NaCl_file_path_1 = '/Users/cal2019/Dropbox/ResearchVoyage/Models_code/HiPOZ/data/20240118/ConductivityData_MgSO4_sat/20240118-120442_Temp294K.txt'
NaCl_file_path_2 = '/Users/cal2019/Dropbox/ResearchVoyage/Models_code/HiPOZ/data/20240118/ConductivityData_MgSO4_NaCl_50-50sat/20240118-121643_Temp294K.txt'
NaCl_file_path_3 = '/Users/cal2019/Dropbox/ResearchVoyage/Models_code/HiPOZ/data/20231214/ConductivityData_NaClSaturatedSoln/20231214-113705_Temp293K.txt'



# Read the data into a pandas DataFrame, skip the header lines, and specify the delimiter
data_1 = pd.read_csv(NaCl_file_path_1, skiprows=9, delimiter='\t')
data_2 = pd.read_csv(NaCl_file_path_2, skiprows=9, delimiter='\t')
data_3 = pd.read_csv(NaCl_file_path_3, skiprows=9, delimiter='\t')


# Convert the pandas Series to numpy arrays
frequencies = np.array(data_1['Frequency (Hz)'])
impedance = np.array(data_1['Impedance (ohm)'])
phase = np.array(data_1['Phase (degrees)'])

# Convert the pandas Series to numpy arrays
frequencies_2 = np.array(data_2['Frequency (Hz)'])
impedance_2 = np.array(data_2['Impedance (ohm)'])
phase_2 = np.array(data_2['Phase (degrees)'])

# Convert the pandas Series to numpy arrays
frequencies_3 = np.array(data_3['Frequency (Hz)'])
impedance_3 = np.array(data_3['Impedance (ohm)'])
phase_3 = np.array(data_3['Phase (degrees)'])

# Combine impedance and phase into a single complex column
Z = impedance * np.exp(1j * np.radians(phase))
Z_2 = impedance_2 * np.exp(1j * np.radians(phase_2))
Z_3 = impedance_3 * np.exp(1j * np.radians(phase_3))

# Plotting impedance vs. frequency
fig, ax = plt.subplots()
ax.plot(frequencies, np.abs(Z), 'o-', label='MgSO4 Saturated')  # Plotting magnitude of impedance
ax.plot(frequencies_2, np.abs(Z_2), 'o-', label='MgSO4-NaCl saturated 50-50 mixture')  # Plotting magnitude of impedance
ax.plot(frequencies_3, np.abs(Z_3), 'o-', label='NaCl Saturated')  # Plotting magnitude of impedance
ax.set_xscale('log')  # Set x-axis to logarithmic scale for frequency
ax.set_yscale('log')  # Set y-axis to logarithmic scale for phase
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Impedance Magnitude (ohm)')
ax.set_title('Impedance vs. Frequency')
ax.legend()

# Plotting phase vs frequency
fig, ax = plt.subplots()
ax.plot(frequencies, phase, 'o-', label='MgSO4 Saturated')  # Plotting magnitude of impedance
ax.plot(frequencies_2, phase_2, 'o-', label='MgSO4-NaCl saturated 50-50 mixture')  # Plotting magnitude of impedance
ax.plot(frequencies_3, phase_3, 'o-', label='NaCl Saturated')  # Plotting magnitude of impedance
ax.set_xscale('log')  # Set x-axis to logarithmic scale for frequency
ax.set_yscale('log')  # Set y-axis to logarithmic scale for phase
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Phase')
ax.set_title('Phase vs. Frequency')
ax.legend()

# Rahmawati 2022 model b circuit
CPE = r'CPE'
R1 = r'R_1'
Rahmawati_model_b = f'R_1-CPE'
Rahmawati_model_b_initial_guess = [1, 4e-6, 9e-1]

CustomCircuit = CustomCircuit(initial_guess=Rahmawati_model_b_initial_guess,
                              circuit=Rahmawati_model_b)


CustomCircuit.fit(frequencies, Z)

CustomCircuit.fit(frequencies, Z_2)
print(CustomCircuit)
CustomCircuit.fit(frequencies, Z_3)


#ax.plot(frequencies, impedance)



#Create a single Nyquist plot for all data and fits
fig, ax = plt.subplots(figsize=(12, 6))  # Adjust the figsize as needed
CustomCircuit.plot(ax=ax, f_data=frequencies, Z_data=Z, kind='nyquist', label='MgSO4 Saturated', color='blue')
CustomCircuit.plot(ax=ax, f_data=frequencies_2, Z_data=Z_2, kind='nyquist', label='MgSO4-NaCl 50-50', color='green')
CustomCircuit.plot(ax=ax, ls='--', f_data=frequencies_3, Z_data=Z_3, kind='nyquist', label='NaCl Saturated', color='red')

#Add labels and legend
ax.set_xlabel('Real Imp (ohm)', fontsize=14)
ax.set_ylabel('-Imaginary Imp (ohm)', fontsize=14)
ax.set_title('Nyquist Plot for Different Conditions', fontsize=14)
ax.legend(fontsize=14)

# Set the aspect ratio to 'auto'
ax.set_aspect('0.2')

# Set x and y axis limits
ax.set_xlim(0, 140)
ax.set_ylim(0, 650)

plt.show()
# # CustomCircuit.plot(f_data=frequencies, Z_data=Z, kind='bode')
# plt.title('Bode Plot MgSO4 Saturated')
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Magnitude (ohm)')
# CustomCircuit.plot(f_data=frequencies_2, Z_data=Z_2, kind='bode')
# plt.title('Bode Plot MgSO4-NaCl 50-50')
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Magnitude (ohm)')
# CustomCircuit.plot(f_data=frequencies_3, Z_data=Z_3, kind='bode')
# plt.title('Bode Plot NaCl Saturated')
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Magnitude (ohm)')

#

