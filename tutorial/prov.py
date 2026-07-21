import meep as mp
import numpy as np
import matplotlib.pyplot as plt

cell = mp.Vector3(16,8,0)

geometry = [mp.Block(mp.Vector3(mp.inf,1,mp.inf),
                     center=mp.Vector3(),
                     material=mp.Medium(epsilon=12))]

sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
                     component=mp.Ez,
                     center=mp.Vector3(-7,0))]

pml_layers = [mp.PML(1.0)]
resolution = 10

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

sim.run(mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)), until=200)


sim.plot2D(fields=mp.Ez,
           eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none', 'alpha':0.3},
           field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
           boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
           output_plane=mp.Volume(center=mp.Vector3(), size=cell))

#eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
#ez_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Ez)
#plt.figure()
#plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
#plt.imshow(ez_data.transpose(), interpolation='spline36', cmap='RdBu', alpha=0.9)
#plt.axis('off')

plt.show()