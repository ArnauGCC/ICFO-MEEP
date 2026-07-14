import meep as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from packages.figures import *
from matplotlib.ticker import MultipleLocator


"""EXECUTION"""
run_meep = True
run_until_source = False
compute_dft_fields = True
set_manual_inc_angle = False


""" PARAMETERS """
alpha_deg = 45                          # triangle angle (base - side) degrees
n = 2.25                                # index prism
dpml = 5                                # thickness of PML
pad = 0.5                               # pad between prism and waveguide
prism_length = 120                      # length of prism
offsx = -1.2*prism_length/2             # offset of prism from the center of the cell (== 0 --> left vertex of the prism in the center of the cell)
offsy = -4                              # offset y-axis
offs_deg = +20                          # offset in degrees from critical angle

#wvlength = 1
#wg_width = 1                   # waveguide width

#fmax = 1.25                                
#fmin = 0.5
theta_inc = 00                          # if set_manual_inc_angle theta_inc is set


#fcen = (fmax+fmin)/2                    # pulse center frequency
#df =  fmax - fmin
fcen = 1/3
df = 0.75
wg_width = 1/(fcen*n)

if not set_manual_inc_angle:
    alpha = math.radians(alpha_deg)         # triangle angle
    theta_inc = math.asin( math.sin(math.asin(1.0/n) - alpha) * n) + alpha
    theta_inc += math.radians(offs_deg)

k = mp.Vector3(1).rotate(mp.Vector3(0, 0, -1), theta_inc)

sx = prism_length + 2*dpml                              # cell size x-axis
sy = (int)(1.5 * prism_length * math.atan(alpha))       # cell size y-axis (with 1.5 scale margin)
prism = create_ideal_prism(alpha_deg, n, mp.Vector3(offsx, offsy), sx)
wg_y = -pad - wg_width/2 + offsy
wg = create_h_waveguide(wg_y, wg_width, n)
src_size = mp.Vector3(y=sy/6 - offsy  if offsx > -sx/2 
                         else (sy/2 - offsy - dpml - (-offsx - sx/2 + dpml - offsy)*math.tan(alpha))*0.5)
src_center = mp.Vector3(-sx/2 + dpml, sy/2 - src_size.y/2 - dpml)

src = [mp.GaussianBeamSource(
        src=mp.GaussianSource(fcen, fwidth=df),
        center=src_center,
        size=src_size,
        beam_x0=src_center + sx*k/4,                                 # relatiu al centre de la font
        beam_kdir=k,
        beam_w0=0.8,                                      # beam waist
        beam_E0=mp.Vector3(0, 0, 1),
        )]

resolution = n*8*fcen
sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
                    geometry=[prism, wg],
                    sources=src,
                    resolution=resolution,
                    boundary_layers=[mp.PML(dpml)]
                    )

if run_meep:

    nfreq = 500
    dft_pt = mp.Vector3(sx/2 - dpml, wg_y)
    if compute_dft_fields: 
        dft_pad = 5
        dft_region_src = sim.add_dft_fields([mp.Ez], fcen, df, nfreq, 
                                            center=src_center, 
                                            size=src_size) 
        dft_region_wg = sim.add_dft_fields([mp.Ez], fcen, df, nfreq, 
                                           center=dft_pt - mp.Vector3(x=dft_pad), 
                                           size=mp.Vector3(y=wg_width))
        dft_region_refl_side = sim.add_dft_fields([mp.Ez], fcen, df, nfreq,
                                                  center=mp.Vector3(sx/2-dpml - dft_pad, sy/4 + offsy/2),
                                                  size=mp.Vector3(y=sy/2-offsy - 2*dft_pad))
        dft_region_refl_top = sim.add_dft_fields([mp.Ez], fcen, df, nfreq,
                                                 center=mp.Vector3(sx/4 - dpml - dft_pad, sy/2 - dft_pad - dpml),
                                                 size=mp.Vector3(x=sx/2))

    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.at_every(1, mp.to_appended("ez_ini", mp.output_efield_z)),
            until_after_sources=0)
    
    if compute_dft_fields:
                src_SPD = np.zeros(nfreq)              # Source Spectral Power Density

                for i in range(nfreq):
                        E = sim.get_dft_array(dft_region_src, mp.Ez, i)
                        src_SPD[i] = np.sum(np.abs(E)**2)

                freqs = np.linspace(fcen-df/2, fcen+df/2,nfreq)


    if not run_until_source:
        sim.run(mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)),
                until=350)
        
        if compute_dft_fields:
                wg_SPD = np.zeros(nfreq)
                refl_SPD = np.zeros(nfreq)

                for i in range(nfreq):
                        Ewg = sim.get_dft_array(dft_region_wg, mp.Ez, i)
                        wg_SPD[i] = np.sum(np.abs(Ewg)**2)

                        Erefl = sim.get_dft_array(dft_region_refl_side, mp.Ez, i)
                        Erefl = np.concatenate((
                                                Erefl,
                                                sim.get_dft_array(dft_region_refl_top, mp.Ez, i)
                                                ))
                        refl_SPD[i] = np.sum(np.abs(Erefl)**2)

                #[x,y,z,w] = sim.get_array_metadata(dft_cell=dft_region)


if compute_dft_fields and run_meep and mp.am_master():
    if not run_until_source:
        #:float.2f
        np.savez(
        f"TMP-dft-off{offs_deg}_fcen{fcen}_w{wg_width}_al{alpha_deg}_n{n}_pad{pad}.npz",
        freqs=freqs,
        src_power=src_SPD,
        wg_power=wg_SPD,
        refl_power=refl_SPD,
        resolution=resolution,
        fcen=fcen,
        df=df,
        prism_length=prism_length
        )

    plt.figure()
    plt.plot(freqs,src_SPD)
    plt.xlabel(r'frequency $f (kHz)$')
    plt.ylabel("Fourier transform (SPD)")
    plt.title("Source SPD")
    plt.grid(True, which="both", alpha=0.3)
    plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))
    
    if not run_until_source:
        plt.figure()
        plt.plot(freqs,wg_SPD)
        plt.xlabel(r'frequency $f (kHz)$')
        plt.ylabel("Fourier transform (SPD)")
        plt.title("Waveguide SPD")
        plt.grid(True, which="both", alpha=0.3)
        plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))
        
        plt.figure()
        plt.plot(freqs,refl_SPD)
        plt.xlabel(r'frequency $f (kHz)$')
        plt.ylabel("Fourier transform (SPD)")
        plt.title("Reflction SPD")
        plt.grid(True, which="both", alpha=0.3)
        plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))
        
    plt.show()


sim.plot2D(fields=mp.Ez,
        eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
        field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
        boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
        output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)))

plt.show()