from impedance import preprocessing
from impedance.models.circuits import Randles, CustomCircuit
import matplotlib.pyplot as plt
from impedance.visualization import plot_nyquist
import pandas as pd
import numpy as np



# 1. Initialize equivalent circuits
randles = Randles(initial_guess=[.01, .005, .001, 200, .1])
randlesCPE = Randles(initial_guess=[.01, .005, .001, 200, .1, .9], CPE=True)
customCircuit = CustomCircuit(initial_guess=[.01, .005, .1, .005, .1, .001, 200],
                              circuit='R_0-p(R_1,C_1)-p(R_2,C_2)-Wo_1')
customConstantCircuit = CustomCircuit(initial_guess=[None, .005, .1, .005, .1, .001, None],
                                      constants={'R_0': 0.02, 'Wo_1_1': 200},
                                      circuit='R_0-p(R_1,C_1)-p(R_2,C_2)-Wo_1')

print(customConstantCircuit)

# 2. Formulate Data
#Several convenience functions for importing data exist in the
# impedance.preprocessing module, including one for reading
# simple .csv files where frequencies are stored in the first column,
# real parts of the impedance are in the second column,
# and imaginary parts of the impedance are in the third column.

# Define the file path for NaCl saturated
NaCl_file_path = '/Users/cal2019/Dropbox/ResearchVoyage/Models_code/HiPOZ/data/20231212/ConductivityData_1mol_kg_NaCl/20231212-150318_Temp294K.txt'

# Define the file path for KCl standard
KCl_file_path = '/Users/cal2019/Dropbox/ResearchVoyage/Models_code/HiPOZ/data/20231214/ConductivityData_KClStd_23us_cm/20231214-104516_Temp293K.txt'

# Read the data into a pandas DataFrame, skip the header lines, and specify the delimiter
data = pd.read_csv(NaCl_file_path, skiprows=9, delimiter='\t')

# # Display the first few rows of the DataFrame to verify the data was loaded correctly
print(data.head())

# Convert the pandas Series to numpy arrays
frequencies = np.array(data['Frequency (Hz)'])
impedance = np.array(data['Impedance (ohm)'])
phase = np.array(data['Phase (degrees)'])

# Combine impedance and phase into a single complex column
Z = impedance * np.exp(1j * np.radians(phase))


# 3. Fit the equivalent circuits to a spectrum
randles.fit(frequencies, Z)
randlesCPE.fit(frequencies, Z)
customCircuit.fit(frequencies, Z)
customConstantCircuit.fit(frequencies, Z)

print(customConstantCircuit)

# 4. Predict circuit model and visualize using plotting
# method included in the package
randles.plot(f_data=frequencies, Z_data=Z, kind='nyquist')
randlesCPE.plot(f_data=frequencies, Z_data=Z, kind='nyquist')
customCircuit.plot(f_data=frequencies, Z_data=Z, kind='nyquist')
customConstantCircuit.plot(f_data=frequencies, Z_data=Z, kind='nyquist')

randles.plot(f_data=frequencies, Z_data=Z, kind='bode')
randlesCPE.plot(f_data=frequencies, Z_data=Z, kind='bode')
customCircuit.plot(f_data=frequencies, Z_data=Z, kind='bode')
customConstantCircuit.plot(f_data=frequencies, Z_data=Z, kind='bode')


plt.show()



