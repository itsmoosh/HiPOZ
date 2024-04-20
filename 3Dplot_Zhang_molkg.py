import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

DO_MOLKG = False
model= 'SPCE'

def Molal2ppt(b_molkg, m_gmol):
    """ Convert dissolved salt concentration from molality to ppt

        Args:
            b_molkg (float, shape N): Concentration(s) in mol/kg (molal)
            m_gmol (float, shape 1 or shape N): Molecular weight of solute in g/mol
        Returns:
            w_ppt (float, shape N): Corresponding concentration(s) in ppt by mass
    """
    m_kgmol = m_gmol / 1e3
    w_ppt = b_molkg*m_kgmol / (1 + b_molkg*m_kgmol) * 1e3

    return w_ppt


def Ppt2molal(w_ppt, m_gmol):
    """ Convert dissolved salt concentration from ppt to molal

        Args:
            w_ppt (float, shape N): Corresponding concentration(s) in ppt by mass
            m_gmol (float, shape 1 or shape N): Molecular weight of solute in g/mol
        Returns:
            b_molkg (float, shape N): Concentration(s) in mol/kg (molal)
    """
    m_kgmol = m_gmol / 1e3
    w_frac = w_ppt / 1e3
    b_molkg = w_frac / (1 - w_frac) / m_kgmol

    return b_molkg

# Read the data from the Excel file
df = pd.read_excel('ForPythonPlots.xlsx', sheet_name='SummaryOurData')

# Filter the data for SPCE method
if model == 'both':
    df_filtered = pd.concat((df[df['Method'] == 'SPCE'], df[df['Method'] == 'TIP4P']))
    modelName = 'Combined SPCE+TIP4P'
else:
    df_filtered = df[df['Method'] == model]
    modelName = model

# Create a 3D scatter plot with colormap for wtpct
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Define colormap and normalize it
colormap = plt.cm.coolwarm
normalize = plt.Normalize(vmin=min(df_filtered['wtpct']), vmax=max(df_filtered['wtpct']))

# Plot the 3D scatter plot
sc = ax.scatter(df_filtered['T_K'], df_filtered['P_MPa'], df_filtered['Cond_S_m'], c=df_filtered['wtpct'], cmap=colormap, marker='o', label=f'{modelName} data')

# Define custom equation
def custom_equation(data, A_1, A_2, A_3, A_4, A_5, n):
    T, P, wtpct = data
    Z = (A_1*T + A_2) * wtpct**n * np.exp((-A_3*wtpct + A_5*P) / (T - A_4)) # Zhang 2020 with P modification
    return Z

# Fit the custom equation to the data
T = df_filtered['T_K']
P = df_filtered['P_MPa']
wtpct = df_filtered['wtpct']
m_molkg = Ppt2molal(wtpct*10, 58.44)
density = df_filtered['Dens_kg_m3']
Z_actual = df_filtered['Cond_S_m']
initial_guess = [1.818, -442.0, 160.4, -616.1, 1, 1.0]  # A_1, A_2, A_3, A_4, A_5, n; values from Zhang et al. (2020) for mostly-NaCl solution
# initial_guess = [3.095e-2, -7.259, -6.451e7, 1.382e9, 2.871e5, 1.022] # Catherine values
# initial_guess = [0.18, -41.7, -96.2, 63.9, -0.05, 1.2] # mol/kg values
# initial_guess = [0.02, -4.97, 18.3, 62.3, -0.05, 1.26] # g/kg values

if DO_MOLKG:
    fitVar = m_molkg
    fitType = 'mol/kg'
else:
    fitVar = wtpct
    fitType = 'wt%'
popt, pcov = curve_fit(custom_equation, (T, P, fitVar), Z_actual, p0=initial_guess)

A_1, A_2, A_3, A_4, A_5, n = popt
uncA_1, uncA_2, uncA_3, uncA_4, uncA_5, uncn = np.sqrt(np.diag(pcov))
print(f'{modelName} fit parameters with m in {fitType}:')
print(f"A_1: {A_1} +/- {uncA_1}")
print(f"A_2: {A_2} +/- {uncA_2}")
print(f"A_3: {A_3} +/- {uncA_3}")
print(f"A_4: {A_4} +/- {uncA_4}")
print(f"A_5: {A_5} +/- {uncA_5}")
print(f"n: {n} +/- {uncn}")

print(f'Covariance matrix:')
lines = [f' & '.join([f'\\num{{{val:.3e}}}' for val in line]) for line in pcov]
matrix = f'\\\\ \n        '.join(lines)
pcovStr = f"""
    $\\begin{{bmatrix}}
        {matrix}\\\\
    \\end{{bmatrix}}$ \\\\
"""
print(pcovStr)

# Compute R-squared value
Z_pred = custom_equation((T, P, fitVar), A_1, A_2, A_3, A_4, A_5, n)
meanZ = np.mean(Z_actual)
N = np.size(Z_actual)
RMSE = np.sqrt(np.sum((Z_pred - Z_actual)**2) / N)
RMSEmean = np.sqrt(np.sum((Z_actual - meanZ)**2) / N)
r_squared = 1 - (RMSE/RMSEmean)**2
print(f"R-squared value: {r_squared}")
print(f'RMS error: {RMSE}')

fit = ax.scatter(df_filtered['T_K'], df_filtered['P_MPa'], Z_pred, c='k', label='Fit')


# Set labels and title
ax.set_title('NaCl(aq) at Icy Moon Subsurface Ocean Conditions', x=0.6, y=0.9)  # Adjust the 'y' value to lower the title
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (MPa)')
ax.set_zlabel('Conductivity (S/m)')
ax.legend(loc='best')

# Add colorbar
sm_obj = plt.cm.ScalarMappable(cmap=colormap, norm=normalize)
sm_obj.set_array([])
cbar = plt.colorbar(sm_obj, label='Concentration(wt%)', ax=ax, shrink=0.5)  # Adjust the 'shrink' parameter to make colorbar shorter

# Show plot
plt.show()
