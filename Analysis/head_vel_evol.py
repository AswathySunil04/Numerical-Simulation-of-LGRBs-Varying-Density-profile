import numpy as np
import pandas as pd
import pyPLUTO.pload as pp
from pyPLUTO import *

# Configuration
wdir = '/home/krishna/GRBs_JET/Spin_profile_1/'
file_indices = range(1, 295)
search_radius = 1000
ix = 1  # fixed

# Results storage
results = []
prev_iy = None

for idx in file_indices:
    try:
        data = pp.pload(idx, w_dir=wdir, datatype='dbl')
    except Exception as e:
        print(f"[{idx:03}] Failed to load file: {e}")
        continue

    prs_line = data.prs[ix, :]

    if prev_iy is None:
        iy = np.argmax(prs_line)  # First frame: full search
    else:
        # Local search, but only allow iy >= prev_iy
        lower = max(prev_iy, 0)
        upper = min(prev_iy + search_radius + 1, len(prs_line))
        local_window = prs_line[lower:upper]
        if len(local_window) > 0:
            iy_local = np.argmax(local_window)
            iy = lower + iy_local
        else:
            iy = prev_iy  # fallback

    # Enforce monotonicity
    if prev_iy is not None and iy < prev_iy:
        iy = prev_iy

    # Collect values
    x_pos = data.x1[ix]
    y_pos = data.x2[iy]
    vx1 = data.vx1[ix, iy]
    vx2 = data.vx2[ix, iy]
    gamma_beta = np.sqrt(vx1**2 + vx2**2)
    lorentz_factor = np.sqrt(1 + gamma_beta**2)
    pressure = data.prs[ix, iy]
    density = data.rho[ix, iy]

    results.append((idx, ix, iy, x_pos, y_pos, vx1, vx2, lorentz_factor, pressure, density))

    print(f"[{idx:03}] ix = {ix}, iy = {iy}, r = {x_pos:.2e}, z = {y_pos:.2e}, "
          f"vx1 = {vx1:.3f}, vx2 = {vx2:.3f}, "
          f"Lorentz = {lorentz_factor:.3f}, P = {pressure:.3e}, ρ = {density:.3e}")

    prev_iy = iy

# Save output
df = pd.DataFrame(results, columns=[
    'TimeIndex', 'ix', 'iy', 'r', 'z', 'vx1', 'vx2', 'Lorentz', 'Pressure', 'Density'
])

iy_1= df['iy']
y_pos_11 = df['z']
df.to_csv("jet_head_monotonic_iy_0.01.csv", index=False)


import numpy as np


y_pos_all = [y_pos_11]
denominator = 0.1 * 3e10  
# Create dictionaries to store beta and gamma values independently for each j
beta_results = {}
gamma_results = {}

# Loop over each y_pos_j (for each j from 1 to 7)
for j, y_pos_j in enumerate(y_pos_all, start=1):
    beta_values_j = []
    gamma_values_j = []

    for i in range(1, len(y_pos_j)):  # i from 2 to 7
        beta = (y_pos_j[i] - y_pos_j[i - 1]) / denominator
        gamma = 1 / np.sqrt(1 - beta**2)
        beta_values_j.append(beta)
        gamma_values_j.append(gamma)

    # Store results for each j independently
    beta_results[f"beta_{j}"] = beta_values_j
    gamma_results[f"gamma_{j}"] = gamma_values_j

    # Optional: Print the results
    print(f"\nFor j = {j}:")
    for i_val, (b, g) in enumerate(zip(beta_values_j, gamma_values_j), start=2):
        print(f"  i = {i_val}: beta = {b:.6f}, gamma = {g:.6f}")

# Now you have the beta and gamma values stored in dictionaries
# Example access:
print("\nAll Beta Results for j=1:", beta_results["beta_1"])
print("\nAll Gamma Results for j=1:", gamma_results["gamma_1"])



beta_72 =  beta_results["beta_1"]
gamma_72 = gamma_results["gamma_1"]

Step = np.arange(1,295)
Step_2 = (Step[:-1] + Step[1:]) / 2





plt.figure(figsize=(8, 6.5))
# plt.plot(steps1/10, gamma_1 ,c='blue')
plt.semilogy(Step_2/10, beta_72, c='red')
plt.ylabel(r'$\beta_{Jet}$', fontsize=22, color='black')
plt.xlabel(r'Time [sec]', fontsize=22, color='black')
plt.legend()
plt.xticks(color='black')
plt.yticks(color='black')
plt.tight_layout()
plt.xlim(0, 6)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('Temporal_evolution_beta_1_semilogy.png', dpi=480, bbox_inches='tight')
plt.show()