import pandas as pd
import matplotlib.pyplot as plt

beta_clean = np.clip(beta_72, 0, 0.999)
beta_series = pd.Series(beta_72) 
window_size = 7
beta_smoothed = beta_series.rolling(window=window_size,min_periods=1, center=True).mean()
 
plt.figure(figsize=(8, 6.5)) 
plt.semilogy(Step_2/10, beta_clean, c='red', label='Raw Data') 
plt.semilogy(Step_2/10, beta_smoothed, c='blue', linewidth=1, 
             label=f'Moving Average (window={window_size})')
plt.ylabel(r'$\beta_{Jet}$', fontsize=22, color='black')
plt.xlabel(r'Time [sec]', fontsize=22, color='black') 
plt.xlim(0, 6)

plt.legend(fontsize=12)
plt.xticks(color='black')
plt.yticks(color='black')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()

plt.savefig('Temporal_evolution_beta_smoothed.png', dpi=480, bbox_inches='tight')
plt.show()