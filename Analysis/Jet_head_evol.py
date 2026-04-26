import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
 
time_sec = df['TimeIndex'].values * 0.1
pos_cm = df['z'].values
C_LIGHT = 3e10 
 
mask_window = (time_sec >= 0) & (time_sec <= 5.0)
t_fit = time_sec[mask_window]
y_fit = pos_cm[mask_window]

slope, intercept, r_value, p_value, std_err = linregress(t_fit, y_fit)
beta_window = slope / C_LIGHT

plt.figure(figsize=(10, 6)) 
plt.plot(t_fit, y_fit, 'r-', linewidth=1.5, label='Breakout Window (0-5s)') 
plt.plot(t_fit, slope * t_fit + intercept, color='blue', linestyle='--', 
         linewidth=2, label=f'Linear Fit: $\\beta$ = {beta_window:.3f}') 
plt.axvline(x=4.3, color='black', linestyle=':', label='Physical Breakout (4.3s)') 
plt.title('Jet Head position vs Time', fontsize=14)
plt.xlabel('Time (seconds)', fontsize=12)
plt.ylabel('Jet Head Position [cm]', fontsize=12)
plt.xlim(0, 5.5)   
plt.ylim(0, 1.4e11)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

plt.tight_layout()
plt.show()  
print(f"Calculated Velocity (beta): {beta_window:.4f} c")
print(f"Fit Quality (R-squared): {r_value**2:.5f}")