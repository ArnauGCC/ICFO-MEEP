import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import math

theta_in = np.arange(0,85,5)
kxs = np.empty((50,theta_in.size))
thetas = np.empty((50,theta_in.size))
Rmeep = np.empty((50,theta_in.size))

for j in range(theta_in.size):
  f = np.genfromtxt("flux_t{}.dat".format(theta_in[j]), delimiter=",")
  kxs[:,j] = f[:,0]
  thetas[:,j] = f[:,2]
  Rmeep[:,j] = f[:,3]

wvl = f[:,1]
# create a 2d matrix for the wavelength by repeating the column vector for each angle
wvls = np.matlib.repmat(np.reshape(wvl, (wvl.size,1)),1,theta_in.size)

n1=1
n2=3.5

# compute angle of refracted planewave in medium n2
# for incident planewave in medium n1 at angle theta_in
theta_out = lambda theta_in: math.asin(n1*math.sin(theta_in)/n2)

# compute Fresnel reflectance for P-polarization in medium n2
# for incident planewave in medium n1 at angle theta_in
Rfresnel = lambda theta_in: math.fabs((n1*math.cos(theta_out(theta_in)) -
                                       n2*math.cos(theta_in)) /
                                      (n1*math.cos(theta_out(theta_in)) +
                                       n2*math.cos(theta_in)))**2

Ranalytic = np.empty((50, theta_in.size))
for m in range(wvl.size):
    for n in range(theta_in.size):
        Ranalytic[m,n] = Rfresnel(math.radians(thetas[m,n]))

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(9,6))

wvl_idx = 0
axs[0, 0].plot(theta_in,Rmeep[wvl_idx,:], 'bo-', label="meep")
axs[0, 0].plot(theta_in,Ranalytic[wvl_idx,:], 'ro-', label="analytic")
axs[0, 0].set_xlabel("angle of incident planewave (degrees)")
axs[0, 0].set_ylabel(f"reflectance at $\\lambda$ = {wvl[wvl_idx]:.1f} μm")
axs[0, 0].legend()

pc = axs[0, 1].pcolormesh(kxs, wvls, Rmeep, cmap='hot_r', shading='gouraud',
                          vmin=0, vmax=Rmeep.max())
axs[0, 1].axis([kxs[0,0], kxs[0,-1], wvl[-1], wvl[0]])
axs[0, 1].set_yticks([t for t in np.arange(wvl[-1], wvl[0] + 0.1, 0.1)])
axs[0, 1].set_xlabel("Bloch-periodic wavevector ($k_x/2π$)")
axs[0, 1].set_ylabel("wavelength (μm)")
axs[0, 1].set_title("reflectance (meep)")
cbar = fig.colorbar(pc, ax=axs[0, 1])
cbar.set_ticks([t for t in np.arange(0, 0.4, 0.1)])
cbar.set_ticklabels(["{:.1f}".format(t) for t in np.arange(0, 0.4, 0.1)])

pc = axs[1, 0].pcolormesh(thetas, wvls, Rmeep, cmap='hot_r', shading='gouraud',
                          vmin=0, vmax=Rmeep.max())
axs[1, 0].axis([thetas.min(), thetas.max(), wvl[-1], wvl[0]])
axs[1, 0].set_xticks([t for t in range(0,100,20)])
axs[1, 0].set_yticks([t for t in np.arange(wvl[-1], wvl[0] + 0.1, 0.1)])
axs[1, 0].set_xlabel("angle of incident planewave (degrees)")
axs[1, 0].set_ylabel("wavelength (μm)")
axs[1, 0].set_title("reflectance (meep)")
cbar = fig.colorbar(pc, ax=axs[1, 0])
cbar.set_ticks([t for t in np.arange(0, 0.4, 0.1)])
cbar.set_ticklabels(["{:.1f}".format(t) for t in np.arange(0, 0.4, 0.1)])

pc = axs[1, 1].pcolormesh(thetas, wvls, Ranalytic, cmap='hot_r', shading='gouraud',
                          vmin=0, vmax=Ranalytic.max())
axs[1, 1].axis([thetas.min(), thetas.max(), wvl[-1], wvl[0]])
axs[1, 1].set_xticks([t for t in range(0,100,20)])
axs[1, 1].set_yticks([t for t in np.arange(wvl[-1], wvl[0] + 0.1, 0.1)])
axs[1, 1].set_xlabel("angle of incident planewave (degrees)")
axs[1, 1].set_ylabel("wavelength (μm)")
axs[1, 1].set_title("reflectance (analytic)")
cbar = fig.colorbar(pc, ax=axs[1, 1])
cbar.set_ticks([t for t in np.arange(0, 0.4, 0.1)])
cbar.set_ticklabels(["{:.1f}".format(t) for t in np.arange(0, 0.4, 0.1)])

fig.tight_layout()
fig.savefig('reflectance_angular_spectrum.png', bbox_inches='tight', dpi=150)