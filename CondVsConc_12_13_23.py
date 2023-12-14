import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Read the Excel file
df = pd.read_excel('NaClCondDataDefCompare.xlsx', sheet_name='CondVsConc')

# Group the data by source
grouped_data = df.groupby('Source')

# Create the plot
fig, ax = plt.subplots()

# Plot the data for each source
legend_handles = []
for i, (source, data) in enumerate(grouped_data):
    line = ax.plot(data['wtpct'], data['S_m'], marker='o', markersize=8, linestyle='-', label=source)
    legend_handles.append(Line2D([], [], marker='o', markersize=8, linestyle='-', label=source, color=line[0].get_color()))

# Set labels and title with color matching the plot
ax.set_xlabel('Concentration (wt %)', fontsize=12, color=legend_handles[0].get_color())
ax.set_ylabel('Conductivity (S/m)', fontsize=12, color=legend_handles[0].get_color())
ax.set_title('Conductivity vs Concentration for Each Source', fontsize=14)

# Set axis tick label font size
ax.tick_params(axis='both', labelsize=10)

# Add legend with no lines
ax.legend(handles=legend_handles, fontsize=10, loc='upper left')

# Show the plot
plt.show()
