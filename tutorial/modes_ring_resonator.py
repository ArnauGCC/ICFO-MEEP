import meep as mp
from pathlib import Path
import os
import subprocess

n = 3.4                 # index of waveguide
w = 1                   # width of waveguide
r = 1                   # inner radius of ring
pad = 4                 # padding between waveguide and edge of PML
dpml = 2                # thickness of PML
sxy = 2*(r+w+pad+dpml)  # cell size

c1 = mp.Cylinder(radius=r+w, material=mp.Medium(index=n))
c2 = mp.Cylinder(radius=r)

fcen = 0.15              # pulse center frequency
df = 0.1                 # pulse frequency width
src = mp.Source(mp.GaussianSource(fcen, fwidth=df), mp.Ez, mp.Vector3(r+0.1))

sim = mp.Simulation(cell_size=mp.Vector3(sxy, sxy),
                    geometry=[c1, c2],
                    sources=[src],
                    resolution=10,
                    symmetries=[mp.Mirror(mp.Y)],       # Simertric respecte les y
                                                        # o qualsevol recta
                    boundary_layers=[mp.PML(dpml)])


harminv = mp.Harminv(mp.Ez, mp.Vector3(r+0.1), fcen, df)
sim.run(mp.at_beginning(mp.output_epsilon),
        mp.after_sources(harminv),
        until_after_sources=300)


sim.use_output_directory()
sim.run(mp.at_every(1/fcen/20, mp.output_png(mp.Ez, "-Zc dkbluered")), until=1/fcen)

file_name = Path(__file__).stem     #Folder name for output files --> stem treu el .py
subprocess.run(["magick", "{}-out/{}-ez-*.png".format(file_name, file_name), "{}-out/ez.gif".format(file_name)])

for mode in harminv.modes:
    subprocess.run([f"rm {file_name}-out/*.png"], shell=True)

    sim.reset_meep()

    freq = mode.freq
    df = 0.01

    src = mp.Source(mp.GaussianSource(freq, fwidth=df), mp.Ez, mp.Vector3(r+0.1))
    
    sim = mp.Simulation(cell_size=mp.Vector3(sxy, sxy),
                    geometry=[c1, c2],
                    sources=[src],
                    symmetries=[mp.Mirror(mp.Y)],
                    resolution=10,
                    boundary_layers=[mp.PML(dpml)])
    
    sim.run(until_after_sources=200)

    sim.use_output_directory()
    sim.run(mp.at_every(1/freq/20, mp.output_png(mp.Ez, "-Zc dkbluered")), until=1/freq)
    
    subprocess.run(["magick", "{}-out/{}-ez-*.png".format(file_name, file_name),
                     "{}-out/ez{:.3f}.gif".format(file_name, freq)])

