import meep as mp
import math
import matplotlib.pyplot as plt
from packages.figures import create_ideal_prism

""" PARAMETERS """
alpha_deg = 20                          # triangle angle degrees
n = 2.00                                # index prism
dpml = 1                                # thickness of PML
prism_length = 35                       # length of prism
offs = -6                               # offset of prism from the center of the cell (== 0 --> prism in the center of the cell)
fcen = 0.75                             # pulse center frequency
source_size = prism_length/4            # source size


sxy = prism_length + 2*dpml             # cell size
alpha = math.radians(alpha_deg)         # triangle angle
prism = create_ideal_prism(alpha_deg, n, prism_length, sxy, offs)

theta_inc = math.asin( math.sin(math.asin(1.0/n) - alpha) * n) + alpha

print("theta_inc = ", theta_inc*180/math.pi, " deg")

theta = math.asin(math.sin(theta_inc)/n) + alpha
print("theta_pris = ", theta*180/math.pi, " deg")


"""
El tipus d'ona que s'obtindrà a la base del prisma depèn del factor gamma.

        gamma = 1:      Incidència amb angle crític
        gamma > 1:      Incidència subcrítica (transmissió) (provar vals 1.2 ~ 1.5)
        gamma < 1:      Incidència supercrítica (ona evanescent) (provar vals 0.5 ~ 0.8)
"""

gamma = 0.8
theta_inc *= gamma

k = mp.Vector3(1).rotate(mp.Vector3(0, 0, 1), -theta_inc)

src = [mp.GaussianBeamSource(
        src=mp.GaussianSource(fcen, fwidth=1 * fcen),
        center=mp.Vector3(-sxy/2 + dpml, sxy/2 - source_size/2),
        size=mp.Vector3(0, source_size),
        beam_x0=sxy*k/4,                                # relatiu al centre de la font
        beam_kdir=k,
        beam_w0=1,                                      # beam waist
        beam_E0=mp.Vector3(0, 0, 1),
        )]       

sim = mp.Simulation(cell_size=mp.Vector3(sxy, sxy),
                    geometry=[prism],
                    sources=src,
                    resolution=n*6/fcen,
                    boundary_layers=[mp.PML(dpml)]
                    )


sim.run(mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ez", mp.output_efield_z),
        until=15
        )
        
sim.plot2D(fields=mp.Ez,
           eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
           field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
           boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
           output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sxy, sxy)))

plt.show()


