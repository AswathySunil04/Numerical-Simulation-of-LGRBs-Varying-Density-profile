import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

# 1. Load and clean data (Targeting 'a' instead of 'Ljet')
log_data = []
log_path = '/home/krishna/GRBs_JET/spin_profile_1.out'
with open(log_path, 'r') as f:
    for line in f:
        match = re.search(r"Time=([\d.e+-]+) \| a=([\d.e+-]+) \| Ljet=([\d.e+-]+)", line)
        if match:
            log_data.append([float(match.group(1))/3e10, float(match.group(2))])

df_full = pd.DataFrame(log_data, columns=['Time', 'a_spin']).drop_duplicates().sort_values('Time')

# 2. Filter (0.5 to 7s) and Subsample (Every 15th point for clarity)
df_window = df_full[(df_full['Time'] >= 0.5) & (df_full['Time'] <= 7.0)].copy()
df_sub = df_window.iloc[::15].copy() 

# 3. Stats Prep 
# Spin is very stable; we assume a 0.5% uncertainty for the Chi2 test
df_sub['a_err'] = df_sub['a_spin'] * 0.005  
x_data, y_data, y_err = df_sub['Time'].values, df_sub['a_spin'].values, df_sub['a_err'].values

# 4. Polynomial Fit (Degree 2)
# Since Spin is linear, we fit directly on the values (no log needed)
poly_coeffs = np.polyfit(x_data, y_data, 2) 
poly_func = np.poly1d(poly_coeffs)
y_fit = poly_func(x_data)

# 5. Chi-Squared Calculation
chi_sq = np.sum(((y_data - y_fit) / y_err)**2)
dof = len(y_fit) - 3 
red_chi = chi_sq / dof

# 6. "PRETTY" SQUARE PLOT
fig, ax = plt.subplots(figsize=(8, 8))

# Error Bars (White-filled circles, no caps)
ax.errorbar(x_data, y_data, yerr=y_err, 
            fmt='s', color='black', markerfacecolor='white', markeredgecolor='black',
            markersize=5, elinewidth=0.8, capsize=0, alpha=0.8, 
            label='Simulated Spin ($a$)')

# Fit Line
ax.plot(x_data, y_fit, color='navy', linewidth=2.5, 
        label=f'Spin Evolution Fit (Deg 2)')

# Formatting (Linear Y-axis)
ax.set_ylim(min(y_data)*0.98, max(y_data)*1.02)
ax.set_xlim(0, 7.5)
ax.tick_params(axis='both', which='both', width=1.2, direction='in')
ax.tick_params(which='major', length=10)
ax.tick_params(labelsize=15)

ax.set_ylabel('Black Hole Spin ($a$)', fontsize=21)
ax.set_xlabel('Time (s)', fontsize=21)

# Stats Box
stats_text = (f"$\chi^2 = {chi_sq:.2f}$\n"
              f"$DOF = {dof}$\n"
              f"$\chi^2_{{red}} = {red_chi:.2f}$")

ax.text(0.05, 0.05, stats_text, transform=ax.transAxes, fontsize=14,
        bbox=dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.9))

ax.legend(fontsize=12, frameon=False, loc='upper right')
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.show()