import meep as mp
import math
import matplotlib.pyplot as plt
from packages.figures import *


""" PARAMETERS """
alpha_deg = 15                          # triangle angle degrees
n = 2.00                                # index prism
dpml = 1                                # thickness of PML
pad = 0.5                               # pad between prism and waveguide
wg_width = 1                            # waveguide width
prism_length = 90                       # length of prism
offsx = -6                              # offset of prism from the center of the cell (== 0 --> prism in the center of the cell)
fcen = 0.75                             # pulse center frequency
source_size = prism_length/4            # source size
gamma = 0.65


sxy = prism_length + 2*dpml             # cell size
prism = create_ideal_prism(alpha_deg, n, prism_length, sxy, offsx)
wg = create_h_waveguide(-pad - wg_width/2, wg_width, n)


alpha = math.radians(alpha_deg)         # triangle angle
theta_inc = math.asin( math.sin(math.asin(1.0/n) - alpha) * n) + alpha
theta_inc *= gamma

k = mp.Vector3(1).rotate(mp.Vector3(0, 0, -1), theta_inc)

src = [mp.GaussianBeamSource(
        src=mp.GaussianSource(fcen, fwidth=1 * fcen),
        center=mp.Vector3(-sxy/2 + dpml, sxy/2 - source_size/2 - 12),
        size=mp.Vector3(0, source_size),
        beam_x0=sxy*k/4,                                # relatiu al centre de la font
        beam_kdir=k,
        beam_w0=1,                                      # beam waist
        beam_E0=mp.Vector3(0, 0, 1),
        )]


sim = mp.Simulation(cell_size=mp.Vector3(sxy, sxy),
                    geometry=[prism, wg],
                    sources=src,
                    resolution=n*6/fcen,
                    boundary_layers=[mp.PML(dpml)]
                    )

sim.run(mp.at_beginning(mp.output_epsilon),
        mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)),
        until=160
        )
        
sim.plot2D(fields=mp.Ez,
           eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
           field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
           boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
           output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sxy, sxy)))

plt.show()