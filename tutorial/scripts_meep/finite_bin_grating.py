"""
import meep as mp
import math
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

resolution = 25        # pixels/μm

dpml = 1.0             # PML thickness
dsub = 3.0             # substrate thickness
dpad = 3.0             # padding between grating and PML
gp = 10.0              # grating period
gh = 0.5               # grating height
gdc = 0.5              # grating duty cycle

nperiods = 10          # number of unit cells in finite periodic grating

ff_distance = 1e8      # far-field distance from near-field monitor
ff_angle = 20          # far-field cone angle
ff_npts = 500          # number of far-field points

ff_length = ff_distance*math.tan(math.radians(ff_angle))
ff_res = ff_npts/ff_length

sx = dpml+dsub+gh+dpad+dpml
cell_size = mp.Vector3(sx)

pml_layers = [mp.PML(thickness=dpml,direction=mp.X)]

symmetries = [mp.Mirror(mp.Y)]

wvl_min = 0.4           # min wavelength
wvl_max = 0.6           # max wavelength
fmin = 1/wvl_max        # min frequency
fmax = 1/wvl_min        # max frequency
fcen = 0.5*(fmin+fmax)  # center frequency
df = fmax-fmin          # frequency width

src_pt = mp.Vector3(-0.5*sx+dpml+0.5*dsub)
sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=src_pt)]

k_point = mp.Vector3()

glass = mp.Medium(index=1.5)

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    k_point=k_point,
                    default_material=glass,
                    sources=sources)

nfreq = 21
n2f_pt = mp.Vector3(0.5*sx-dpml-0.5*dpad)
n2f_obj = sim.add_near2far(fcen, df, nfreq, mp.Near2FarRegion(center=n2f_pt))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, n2f_pt, 1e-9))

    
### OJUUU!!! fa:
### center=mp.Vector3(ff_distance,0.5*ff_length)
### pq està explotant simetria en y axis
### per tant només calcula la meitat de les y --> 0.5*ff_length
    
ff_source = sim.get_farfields(n2f_obj, ff_res, center=mp.Vector3(ff_distance,0.5*ff_length), size=mp.Vector3(y=ff_length))

sim.reset_meep()

### unit cell with periodic boundaries

sy = gp
cell_size = mp.Vector3(sx,sy)

sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=src_pt, size=mp.Vector3(y=sy))]

geometry = [mp.Block(material=glass, size=mp.Vector3(dpml+dsub,mp.inf,mp.inf), center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub))),
            mp.Block(material=glass, size=mp.Vector3(gh,gdc*gp,mp.inf), center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh))]

sim = mp.Simulation(resolution=resolution,
                    split_chunks_evenly=True,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources,
                    symmetries=symmetries)

n2f_obj = sim.add_near2far(fcen, df, nfreq, mp.Near2FarRegion(center=n2f_pt, size=mp.Vector3(y=sy)), nperiods=nperiods)

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, n2f_pt, 1e-9))

ff_unitcell = sim.get_farfields(n2f_obj, ff_res, center=mp.Vector3(ff_distance,0.5*ff_length), size=mp.Vector3(y=ff_length))

sim.reset_meep()

### finite periodic grating with flat surface termination extending into PML

num_cells = 2*nperiods+1
sy = dpml+num_cells*gp+dpml
cell_size = mp.Vector3(sx,sy)

pml_layers = [mp.PML(thickness=dpml)]

sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df, is_integrated=True),
                     component=mp.Ez,
                     center=src_pt,
                     size=mp.Vector3(y=sy))]

geometry = [mp.Block(material=glass, size=mp.Vector3(dpml+dsub,mp.inf,mp.inf), center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub)))]

for j in range(num_cells):
    geometry.append(mp.Block(material=glass,
                             size=mp.Vector3(gh,gdc*gp,mp.inf),
                             center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,-0.5*sy+dpml+(j+0.5)*gp)))

sim = mp.Simulation(resolution=resolution,
                    split_chunks_evenly=True,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources,
                    symmetries=symmetries)

n2f_obj = sim.add_near2far(fcen, df, nfreq, mp.Near2FarRegion(center=n2f_pt, size=mp.Vector3(y=sy-2*dpml)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, n2f_pt, 1e-9))

ff_supercell = sim.get_farfields(n2f_obj, ff_res, center=mp.Vector3(ff_distance,0.5*ff_length), size=mp.Vector3(y=ff_length))

norm_err = LA.norm(ff_unitcell['Ez']-ff_supercell['Ez'])/nperiods
print("error:, {}, {}".format(nperiods,norm_err))


freqs = mp.get_near2far_freqs(n2f_obj)
wvl = np.divide(1,freqs)
ff_lengths = np.linspace(0,ff_length,ff_npts)
angles = [math.degrees(math.atan(f)) for f in ff_lengths/ff_distance]

wvl_slice = 0.5
idx_slice = np.where(np.asarray(freqs) == 1/wvl_slice)[0][0]

rel_enh = np.absolute(ff_unitcell['Ez'])**2/np.absolute(ff_source['Ez'])**2

plt.figure(dpi=150)

plt.subplot(1,2,1)
plt.pcolormesh(wvl,angles,rel_enh,cmap='Blues',shading='auto')
plt.axis([wvl_min,wvl_max,0,ff_angle])
plt.xlabel("wavelength (μm)")
plt.ylabel("angle (degrees)")
plt.grid(linewidth=0.5,linestyle='--')
plt.xticks([t for t in np.arange(wvl_min,wvl_max+0.1,0.1)])
plt.yticks([t for t in range(0,ff_angle+1,10)])
plt.title("far-field spectra")

plt.subplot(1,2,2)
plt.plot(angles,rel_enh[:,idx_slice],'bo-')
plt.xlim(0,ff_angle)
plt.ylim(0)
plt.xticks([t for t in range(0,ff_angle+1,10)])
plt.xlabel("angle (degrees)")
plt.ylabel("relative enhancement")
plt.grid(axis='x',linewidth=0.5,linestyle='--')
plt.title("f.-f. spectra @  λ = {:.1} μm".format(wvl_slice))

plt.tight_layout(pad=0.5)
plt.show()
"""


import meep as mp
import numpy as np
import math
import matplotlib.pyplot as plt

# True:  plot the scattered fields in the air region adjacent to the grating
# False: plot the diffraction spectra based on a 1d cross section of the scattered fields
field_profile = False

resolution = 50         # pixels/μm

dpml = 1.0              # PML thickness
dsub = 2.0              # substrate thickness
dpad = 1.0              # flat-surface padding
gp = 1.0                # grating periodicity
gh = 0.5                # grating height
gdc = 0.5               # grating duty cycle
num_cells = 5           # number of grating unit cells

# air region thickness adjacent to grating
dair = 10 if field_profile else dpad

wvl = 0.5               # center wavelength
fcen = 1/wvl            # center frequency

k_point = mp.Vector3()

glass = mp.Medium(index=1.5)

pml_layers = [mp.PML(thickness=dpml)]

symmetries=[mp.Mirror(mp.Y)]

sx = dpml+dsub+gh+dair+dpml
sy = dpml+dpad+num_cells*gp+dpad+dpml
cell_size = mp.Vector3(sx,sy)

src_pt = mp.Vector3(-0.5*sx+dpml+0.5*dsub)
sources = [mp.Source(mp.GaussianSource(fcen,fwidth=0.2*fcen,is_integrated=True),
                     component=mp.Ez,
                     center=src_pt,
                     size=mp.Vector3(y=sy))]

geometry = [mp.Block(material=glass,
                     size=mp.Vector3(dpml+dsub,mp.inf,mp.inf),
                     center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub)))]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources,
                    symmetries=symmetries)

mon_pt = mp.Vector3(0.5*sx-dpml-0.5*dair)
near_fields = sim.add_dft_fields([mp.Ez], fcen, 0, 1, center=mon_pt, size=mp.Vector3(dair if field_profile else 0,sy-2*dpml))

sim.run(until_after_sources=100)

flat_dft = sim.get_dft_array(near_fields, mp.Ez, 0)

sim.reset_meep()

for j in range(num_cells):
  geometry.append(mp.Block(material=glass,
                           size=mp.Vector3(gh,gdc*gp,mp.inf),
                           center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,-0.5*sy+dpml+dpad+(j+0.5)*gp)))

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources,
                    symmetries=symmetries)

near_fields = sim.add_dft_fields([mp.Ez], fcen, 0, 1, center=mon_pt, size=mp.Vector3(dair if field_profile else 0,sy-2*dpml))

sim.run(until_after_sources=100)

grating_dft = sim.get_dft_array(near_fields, mp.Ez, 0)

scattered_field = grating_dft-flat_dft
scattered_amplitude = np.abs(scattered_field)**2

[x,y,z,w] = sim.get_array_metadata(dft_cell=near_fields)

if field_profile:
  if mp.am_master():
    plt.figure(dpi=150)
    plt.pcolormesh(x,y,np.rot90(scattered_amplitude),cmap='inferno',shading='gouraud',vmin=0,vmax=scattered_amplitude.max())
    plt.gca().set_aspect('equal')
    plt.xlabel('x (μm)')
    plt.ylabel('y (μm)')

    # ensure that the height of the colobar matches that of the plot
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)
    plt.tight_layout()
    plt.show()
else:
  ky = np.fft.fftshift(np.fft.fftfreq(len(scattered_field), 1/resolution))
  FT_scattered_field = np.fft.fftshift(np.fft.fft(scattered_field))
  if mp.am_master():
    plt.figure(dpi=150)
    plt.subplots_adjust(hspace=0.3)

    plt.subplot(2,1,1)
    plt.plot(y,scattered_amplitude,'bo-')
    plt.xlabel("y (μm)")
    plt.ylabel("field amplitude")

    plt.subplot(2,1,2)
    plt.plot(ky,np.abs(FT_scattered_field)**2,'ro-')
    plt.gca().ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    plt.xlabel(r'wavevector k$_y$, 2π (μm)$^{-1}$')
    plt.ylabel("Fourier transform")
    plt.gca().set_xlim([-3, 3])

    plt.tight_layout(pad=1.0)
    plt.show()