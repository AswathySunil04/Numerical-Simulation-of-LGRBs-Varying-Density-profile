import pyPLUTO as pypl
import pyPLUTO.pload as pp
import pyPLUTO.Image as img
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
plt.rcParams['image.cmap'] = 'jet'
wdir = '/home/krishna/GRBs_JET/Spin_profile_1/'

# Choose the time step you want to visualize (replace 'your_desired_step' with the actual step)
desired_step = 53
# Loading the datainto a pload object D.
D = pp.pload(desired_step, w_dir=wdir, datatype='dbl')
I = img.Image()
print(D.vars)
display_object = I.pldisplay(D, np.log10(D.rho), x1=D.x1, x2=D.x2, label1='x-grid', label2='y-grid',
                             title=f'', fontsize = 20,
            figsize=[10,8])



cbar = plt.colorbar(display_object, ax=plt.gca(), pad=0.01)
cbar.set_label(r'$Log_{10}$ [$\rho$]  (gr/cm$^{3})$', fontsize=20)
plt.ylabel('z-grid',fontsize=20)
plt.xlabel('r-grid',fontsize=20)
plt.text(
    0.05, 0.9,
    f'Time Index: {desired_step}',
    transform=plt.gca().transAxes,
    bbox=dict(facecolor='white', alpha=0.85),
    fontsize=14
)
plt.tight_layout()
plt.show()