import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.signal import savgol_filter  # Standard DS smoothing library

# --- 1. Data Extraction (Spin) ---
spin_path = '/home/krishna/GRBs_JET/spin_profile_1.out'
times_list = []
a_list = [] 
try:
    with open(spin_path, 'r') as f:
        for line in f: 
            match = re.search(r"Time=([\d.e+-]+) \| a=([\d.e+-]+)", line)
            if match:
                times_list.append(float(match.group(1)))
                a_list.append(float(match.group(2)))
    t_spin = np.array(times_list)
    a_values = np.array(a_list)
except Exception as e:
    print(f"Failed to parse file: {e}")

# --- 2. Alignment ---
# Assuming beta_72 is your extracted jet head velocity array
time_beta = np.linspace(t_spin.min(), t_spin.max(), len(beta_72))
matched_a = np.interp(time_beta, t_spin, a_values)

# --- 3. Filtering & Cleaning ---
active_idx = 53
a_plot = matched_a[:active_idx]
beta_plot = np.array(beta_72[:active_idx])
beta_clean = np.clip(beta_plot, 1e-3, 0.999) # Avoid 0 for log plot

# --- 4. SMOOTHENING ---
# window_length: must be odd and less than len(beta_clean)
# polyorder: usually 2 or 3 for physical trends
window = 7  
beta_smooth = savgol_filter(beta_clean, window_length=window, polyorder=2)

# --- 5. PLOTTING ---
plt.figure(figsize=(8, 7))

# Plot Raw Data (Faint)
plt.plot(a_plot, beta_clean, c='blue', label='Raw Data')

# Plot Smoothed Data (Bold)
plt.plot(a_plot, beta_smooth, c='red', linewidth= 1, label='Smoothed Trend (Savgol)')

# Formatting
plt.ylabel(r'$\beta_{Jet}$ (Velocity)', fontsize=18)
plt.xlabel(r'$a_{spin}$ (BH Spin)', fontsize=18)
plt.yscale('log')
plt.ylim(1e-2, 1.5) 
plt.grid(True, which='both', linestyle=':', alpha=0.6)

# Data Science Insight: Add an arrow to show time direction
# (Spin usually decreases over time)
plt.annotate('Evolution Over Time', xy=(a_plot[-1], beta_smooth[-1]), 
             xytext=(a_plot[0], beta_smooth[0]*0.5),
             arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=8),
             fontsize=12)

plt.legend(frameon=False, fontsize=12)
plt.tight_layout()

plt.savefig('Beta_vs_A_Smoothed.png', dpi=480)
plt.show()