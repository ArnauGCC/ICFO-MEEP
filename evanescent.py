import meep as mp
import matplotlib.pyplot as plt
import math

dpml = 1        # thickness of PML
sxy = 12        # cell size
freq = 0.8
n_default = 1.5
theta_inc = math.asin(1/n_default)


evanescent = True

if evanescent:
    theta_inc *=0.8
else:
    theta_inc *=1.2 


k_point = mp.Vector3(1, 0).rotate(mp.Vector3(0, 0, -1), theta_inc)

src = [mp.GaussianBeamSource(
                                src=mp.ContinuousSource(frequency=freq),
                                center=mp.Vector3(-sxy/2+dpml, sxy/4),
                                size=mp.Vector3(0, sxy/4),
                                beam_x0=mp.Vector3(0, (sxy/2+dpml)*math.sin(theta_inc)),
                                beam_kdir=k_point,
                                beam_w0=4,          # beam waist
                                beam_E0=mp.Vector3(0, 0, 1)
                            )]

sim = mp.Simulation(cell_size=mp.Vector3(sxy, sxy),
                     geometry=[mp.Block(#center=mp.Vector3(), 
                                        #size=mp.Vector3(mp.inf, 2), 
                                        center=mp.Vector3(0, -sxy/4), 
                                        size=mp.Vector3(sxy, sxy/2),
                                        material=mp.Medium(index=1))],
                     boundary_layers=[mp.PML(dpml)],
                     sources=src,
                     resolution=10,
                     default_material=mp.Medium(index=n_default)
                     )

sim.run(until=20)

sim.plot2D(fields=mp.Ez,
           eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
           field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
           boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
           output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sxy, sxy)))

plt.show()