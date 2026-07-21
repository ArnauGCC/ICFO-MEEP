import meep as mp
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from packages.figures import *
import math
import numpy as np

run_meep = True    

dpml = 1        # thickness of PML
sxy = 18        # cell size
freq = 1
n_prsim = 1.8
n_wg = 1.5
pad = 0.525
wg_width = 1/(2*freq*n_wg)
n_eff =  n_prsim*math.sin(math.radians(45.6))
theta_inc = math.radians(90) - math.asin(n_eff/n_prsim)

prism = create_h_waveguide(sxy/4, sxy/2, n_prsim)
wg = create_h_waveguide(-wg_width/2 -pad/2, wg_width, n_wg)

k_point = mp.Vector3(1, 0).rotate(mp.Vector3(0, 0, -1), theta_inc)

df = freq*0.3
src = [mp.GaussianBeam2DSource(
                                src=mp.GaussianSource(frequency=freq, fwidth=df),
                                center=mp.Vector3(-sxy/2+dpml, sxy/4),
                                size=mp.Vector3(0, sxy/4),
                                #beam_x0=mp.Vector3(0, (sxy/2+dpml)*math.sin(theta_inc)),
                                beam_x0=sxy*k_point/4,
                                beam_kdir=k_point,
                                beam_w0=1,          # beam waist
                                beam_E0=mp.Vector3(0, 0, 1)
                            )]

resolution = n_prsim*16*freq
sim = mp.Simulation(cell_size=mp.Vector3(sxy, sxy),
                     geometry=[prism, wg],
                     boundary_layers=[mp.PML(dpml)],
                     sources=src,
                     resolution=resolution)

pad_fr = 0.2
refl_fr = [
            mp.FluxRegion(center=mp.Vector3(sxy/2 - dpml - pad_fr, sxy/2 - (sxy/6 - dpml - pad_fr)/2 - dpml - pad_fr), size=mp.Vector3(y=sxy/6 - dpml - pad_fr)),
            mp.FluxRegion(center=mp.Vector3(sxy/6 - dpml/2 - pad_fr, sxy/2 - dpml - pad_fr), size=mp.Vector3(x=2*sxy/3 - dpml))
          ]
wg_fr = mp.FluxRegion(center=mp.Vector3(sxy/2 - dpml - pad_fr, -(wg_width + pad)/2), size=mp.Vector3(y=12*wg_width))
#wg_fr = mp.FluxRegion(center=mp.Vector3(sxy/2 - dpml - pad_fr, sxy/4 -(wg_width + pad + dpml)), size=mp.Vector3(y=sxy/2 - pad -dpml))


src_fr = mp.FluxRegion(center=mp.Vector3(-sxy/2 + dpml + pad_fr, sxy/4), size=mp.Vector3(0, sxy/3),)

nfreq = 600

df *= 1.2
wg_region = sim.add_flux(freq, df, nfreq, wg_fr)
refl_region = sim.add_flux(freq, df, nfreq, *refl_fr)
src_region = sim.add_flux(freq, df, nfreq, src_fr)

if run_meep:
    sim.run(mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)), until=90)

    src_flux = mp.get_fluxes(src_region)
    wg_flux = mp.get_fluxes(wg_region)
    refl_flux = mp.get_fluxes(refl_region)
    freqs = mp.get_flux_freqs(refl_region)


    if mp.am_master():
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

        plt.figure()
        plt.plot(freqs,src_flux)
        plt.xlabel(r'frequency $f (kHz)$')
        plt.ylabel("Flux")
        plt.title("Source SPD")
        plt.grid(True, which="both", alpha=0.3)
        plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

        plt.figure()
        plt.plot(freqs,np.divide(refl_flux, src_flux))
        plt.xlabel(r'frequency $f (kHz)$')
        plt.ylabel("Flux")
        plt.title("Reflected/Source SPD")
        plt.grid(True, which="both", alpha=0.3)
        plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

        plt.figure()
        plt.plot(freqs,np.divide(wg_flux, src_flux))
        plt.xlabel(r'frequency $f (kHz)$')
        plt.ylabel("Flux")
        plt.title("Waveguide/Source SPD")
        plt.grid(True, which="both", alpha=0.3)
        plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

        plt.figure()
        plt.plot(freqs,np.add(np.divide(wg_flux, src_flux), np.divide(refl_flux, src_flux)))
        plt.xlabel(r'frequency $f (kHz)$')
        plt.ylabel("Flux")
        plt.title("Ratios sum SPD")
        plt.grid(True, which="both", alpha=0.3)
        plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

        plt.show()
        

sim.plot2D(fields=mp.Ez,
        eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
        field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
        boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
        output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sxy, sxy)))

plt.show()