import meep as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from packages.figures import *
from matplotlib.ticker import MultipleLocator


""" PARAMETERS """
alpha_deg = 15                          # triangle angle degrees
n = 2.00                                # index prism
dpml = 1                                # thickness of PML
pad = 0.5                               # pad between prism and waveguide
wg_width = 1                            # waveguide width
prism_length = 100                      # length of prism
offsx = -16                             # offset of prism from the center of the cell (== 0 --> prism in the center of the cell)
offsy = -8                              # offset y-axis
fmax = 1
fmin = 0.52
gamma = 0.65

fcen = (fmax+fmin)/2                    # pulse center frequency
df =  fmax - fmin

alpha = math.radians(alpha_deg)         # triangle angle
theta_inc = math.asin( math.sin(math.asin(1.0/n) - alpha) * n) + alpha
theta_inc *= gamma

k = mp.Vector3(1).rotate(mp.Vector3(0, 0, -1), theta_inc)


sx = prism_length + 2*dpml                              # cell size x-axis
sy = 1.5 * prism_length * math.atan(alpha)              # cell size y-axis (with 1.5 scale margin)
prism = create_prism(alpha_deg, n, prism_length, offsx, offsy)
wg_y = -pad - wg_width/2 + offsy
wg = create_h_waveguide(wg_y, wg_width, n)
source_size = sy/3                                      # source size


src = [mp.GaussianBeamSource(
        src=mp.GaussianSource(fcen, fwidth=df),
        center=mp.Vector3(-sx/2 + dpml, sy/2 - source_size/2),
        size=mp.Vector3(0, source_size),
        beam_x0=sx*k/4,                                 # relatiu al centre de la font
        beam_kdir=k,
        beam_w0=1,                                      # beam waist
        beam_E0=mp.Vector3(0, 0, 1),
        )]

resolution = n*6/fcen
sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
                    geometry=[prism, wg],
                    sources=src,
                    resolution=resolution,
                    boundary_layers=[mp.PML(dpml)]
                    )


nfreq = 1500
dft_pt = mp.Vector3(prism_length/2 + offsx, wg_y)
dft_region = sim.add_dft_fields([mp.Ez], fcen, df, nfreq, center=dft_pt, size=mp.Vector3(y=wg_width))

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)),
        until=260
        )


spectral_power = np.zeros(nfreq)

for i in range(nfreq):
        E = sim.get_dft_array(dft_region, mp.Ez, i)
        spectral_power[i] = np.sum(np.abs(E)**2)

freqs = np.linspace(fcen-df/2, fcen+df/2,nfreq)

[x,y,z,w] = sim.get_array_metadata(dft_cell=dft_region)


if mp.am_master():
        #:float.2f
        np.savez(
        f"dft-gam{gamma}_fmax{fmax}_w{wg_width}_al{alpha_deg}_n{n}_pad{pad}.npz",
        freqs=freqs,
        power=spectral_power,
        resolution=resolution,
        fcen=fcen,
        df=df
        )

        plt.plot(freqs,spectral_power)
        plt.xlabel(r'frequency $f (kHz)$')
        plt.ylabel("Fourier transform (PDF)")
        plt.grid(True, which="both", alpha=0.3)
        plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))
        plt.show()


sim.plot2D(fields=mp.Ez,
        eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
        field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
        boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
        output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)))

plt.show()