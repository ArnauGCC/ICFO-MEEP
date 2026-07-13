import meep as mp
import matplotlib.pyplot as plt
import numpy as np

resolution = 50  # pixels/μm

cell_size = mp.Vector3(14, 14)

pml_layers = [mp.PML(thickness=2)]

# rotation angle (in degrees) of waveguide, counter clockwise (CCW) around z-axis
rot_angle = np.radians(20)

w = 1.0  # width of waveguide

geometry = [
    mp.Block(
        center=mp.Vector3(),
        size=mp.Vector3(mp.inf, w, mp.inf),
        e1=mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle),
        e2=mp.Vector3(y=1).rotate(mp.Vector3(z=1), rot_angle),
        material=mp.Medium(epsilon=12),
    )
]

fsrc = 0.15  # frequency of eigenmode or constant-amplitude source
bnum = 1  # band number of eigenmode

kpoint = mp.Vector3(x=1).rotate(mp.Vector3(z=1), rot_angle)

#compute_flux = False  # compute flux (True) or plot the field profile (False)

eig_src = True  # eigenmode (True) or constant-amplitude (False) source
"""
if eig_src:
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fsrc, fwidth=0.2 * fsrc),
            #if compute_flux
            #else mp.ContinuousSource(fsrc),
            center=mp.Vector3(),
            size=mp.Vector3(y=3 * w),
            direction=mp.NO_DIRECTION,
            eig_kpoint=kpoint,
            eig_band=bnum,
            eig_parity=mp.EVEN_Y + mp.ODD_Z if rot_angle == 0 else mp.ODD_Z,
            eig_match_freq=True,
        )
    ]
else:
    sources = [
        mp.Source(
            src=mp.GaussianSource(fsrc, fwidth=0.2 * fsrc),
            #if compute_flux
            #else mp.ContinuousSource(fsrc),
            center=mp.Vector3(),
            size=mp.Vector3(y=3 * w),
            component=mp.Ez,
        )
    ]
"""

sources = [mp.GaussianBeamSource(
            src=mp.ContinuousSource(frequency=fsrc),
            center=mp.Vector3(),
            beam_x0=mp.Vector3(-8,0),
            beam_kdir=kpoint,
            beam_w0=3.0,          # beam waist
            component=mp.Ez
        )]

sim = mp.Simulation(
    cell_size=cell_size,
    resolution=resolution,
    boundary_layers=pml_layers,
    sources=sources,
    #geometry=geometry,
    symmetries=[mp.Mirror(mp.Y)] if rot_angle == 0 else [],
)


sim.run(until=100)
sim.plot2D(
    output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(10, 10)),
    fields=mp.Ez,
    field_parameters={"alpha": 0.9},
)
plt.show()