"""Plot broadband fluxes and all five picket-fence optical depths."""
import numpy as np
import matplotlib.pyplot as plt

flux = np.loadtxt('FMS_RC_flx.out')
pressure = flux[:, 1] / 1e5
fig, axes = plt.subplots(1, 2, figsize=(11, 6), sharey=True)
for column, label in enumerate(('SW up', 'SW down', 'SW net', 'LW up', 'LW down', 'LW net'), start=2):
    axes[0].plot(flux[:, column], pressure, label=label)
for column, label in enumerate(('Visible 1', 'Visible 2', 'Visible 3', 'IR 1', 'IR 2'), start=8):
    axes[1].plot(flux[:, column], pressure, label=label)
axes[0].set_xlabel(r'Flux [W m$^{-2}$]')
axes[0].set_ylabel('Pressure [bar]')
axes[1].set_xlabel('Optical depth')
axes[1].set_xscale('log')
for ax in axes:
    ax.set_yscale('log')
    ax.legend()
axes[0].invert_yaxis()
fig.tight_layout()
plt.show()
