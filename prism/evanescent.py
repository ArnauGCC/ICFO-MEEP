import meep as mp
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import math

dpml = 1        # thickness of PML
sxy = 52        # cell size
freq = 0.8
n_default = 1.5
theta_inc = math.asin(1/n_default)

"""
evanescent = True

if evanescent:
    theta_inc *=0.8
else:
    theta_inc *=1.2 
"""

k_point = mp.Vector3(1, 0).rotate(mp.Vector3(0, 0, -1), theta_inc)

df = freq*0.3
src = [mp.GaussianBeam2DSource(
                                src=mp.GaussianSource(frequency=freq, fwidth=df),
                                center=mp.Vector3(-sxy/2+dpml, sxy/4),
                                size=mp.Vector3(0, sxy/4),
                                beam_x0=mp.Vector3(0, (sxy/2+dpml)*math.sin(theta_inc)),
                                beam_kdir=k_point,
                                beam_w0=20,          # beam waist
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

refl_fr = [
            mp.FluxRegion(center=mp.Vector3(sxy/2 - 2*dpml, 3*sxy/10 + dpml), size=mp.Vector3(y=2*sxy/5 - dpml)),
            mp.FluxRegion(center=mp.Vector3(y=sxy/2 - 2*dpml), size=mp.Vector3(x=sxy))
          ]
wg_fr = mp.FluxRegion(center=mp.Vector3(sxy/2 - 2*dpml), size=mp.Vector3(y=3))

nfreq = 100

wg_region = sim.add_flux(freq, df, nfreq, wg_fr)
refl_region = sim.add_flux(freq, df, nfreq, *refl_fr)

sim.run(mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)), until=150)

wg_flux = mp.get_fluxes(wg_region)
refl_flux = mp.get_fluxes(refl_region)
freqs = mp.get_flux_freqs(refl_region)

plt.figure()
plt.plot(freqs,wg_flux)
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Flux")
plt.title("Waveguide SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))


plt.figure()
plt.plot(freqs,refl_flux)
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Flux")
plt.title("Reflected SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))
plt.show()

sim.plot2D(fields=mp.Ez,
           eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
           field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
           boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
           output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sxy, sxy)))

plt.show()