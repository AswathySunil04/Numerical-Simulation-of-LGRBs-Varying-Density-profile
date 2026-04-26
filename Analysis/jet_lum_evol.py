import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

# 1. Load and clean data
log_data = []
log_path = '/home/krishna/GRBs_JET/spin_profile_1.out'
with open(log_path, 'r') as f:
    for line in f:
        match = re.search(r"Time=([\d.e+-]+) \| a=([\d.e+-]+) \| Ljet=([\d.e+-]+)", line)
        if match:
            log_data.append([float(match.group(1))/3e10, float(match.group(3))])

df_full = pd.DataFrame(log_data, columns=['Time', 'Ljet']).drop_duplicates().sort_values('Time')

# 2. INCREASE SUBSAMPLING (Every 15th point for "breathing room")
df_window = df_full[(df_full['Time'] >= 0.5) & (df_full['Time'] <= 7.0)].copy()
df_sub = df_window.iloc[::15].copy() 

# 3. Stats Prep (1.5% error)
df_sub['Ljet_err'] = df_sub['Ljet'] * 0.015  
x_data, y_data, y_err = df_sub['Time'].values, df_sub['Ljet'].values, df_sub['Ljet_err'].values

# 4. Degree 2 Fit
log_y = np.log10(y_data)
poly_coeffs = np.polyfit(x_data, log_y, 2) 
poly_func = np.poly1d(poly_coeffs)
y_fit = 10**poly_func(x_data)

# 5. Stats
chi_sq = np.sum(((y_data - y_fit) / y_err)**2)
dof = len(y_fit) - 3 
red_chi = chi_sq / dof

# 6. "PRETTY" SQUARE PLOT
fig, ax = plt.subplots(figsize=(8, 8))

# PRETTY SETTINGS: capsize=0 (no T-bars), thinner elinewidth, smaller markers
ax.errorbar(x_data, y_data, yerr=y_err, 
            fmt='o', color='black', markerfacecolor='black', markeredgecolor='black',
            markersize=4, elinewidth=0.8, capsize=0, alpha=0.8, 
            label='Simulated Data (Subsampled)')

# Smoother Fit Line
ax.plot(x_data, y_fit, color='darkgreen', linewidth=2.5, 
        label=f'Poly Fit (Deg 2)')

# Formatting (Clean Axis)
ax.set_yscale('log')
ax.set_xlim(0, 7.5)
ax.tick_params(axis='both', which='both', width=1.2, direction='in')
ax.tick_params(which='major', length=10)
ax.tick_params(which='minor', length=5)
ax.tick_params(labelsize=15)

ax.set_ylabel('Jet Luminosity (erg/s)', fontsize=21)
ax.set_xlabel('Time (s)', fontsize=21)

# Cleaner Stats Box
stats_text = (f"$\chi^2 = {chi_sq:.1f}$\n"
              f"$DOF = {dof}$\n"
              f"$\chi^2_{{red}} = {red_chi:.2f}$")

ax.text(0.05, 0.05, stats_text, transform=ax.transAxes, fontsize=14,
        bbox=dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.9))

ax.legend(fontsize=12, frameon=False)
ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.show()