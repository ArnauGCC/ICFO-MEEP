import meep as mp
import matplotlib.pyplot as plt



do_difference = False
prev_flux = None

sxy = 10
dpml = 1
sim = mp.Simulation(cell_size=mp.Vector3(sxy, sxy),
                            sources=[mp.Source(mp.ContinuousSource(frequency=1, slowness=0),
                                                component=mp.Ez,
                                                center=mp.Vector3(-4,0),
                                                size = mp.Vector3(y=4))],
                            resolution=8,
                            boundary_layers=[mp.PML(dpml)]
                            )

src_region = sim.add_energy(1, 0, 1, mp.FluxRegion(center=mp.Vector3(-3.5,0), size=mp.Vector3(y=5)))

def print_flux(sim):
    global prev_flux
    if do_difference:
        flux = sim.get_field_point(mp.Ez, mp.Vector3(-3.5,0))
        if prev_flux:   print("flux_diff =", flux - prev_flux)
        else:           print("flux_diff =", flux)

        prev_flux = flux

    else:   print("flux =", sim.get_field_point(mp.Ez, mp.Vector3(-3.5,0)))

sim.run(
    mp.at_every(1, lambda sim: print_flux(sim)),
    until=2500,
)


sim.plot2D(fields=mp.Ez,
            eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
            field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
            boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
            output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sxy, sxy)))

plt.show()