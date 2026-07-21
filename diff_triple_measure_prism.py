import meep as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from packages.figures import *
from matplotlib.ticker import MultipleLocator
from scipy.optimize import brentq


def f(x, n, d):
    return (
        np.pi * d * np.sqrt(n**2 - x**2)
        - np.arctan(
            np.sqrt((x**2 - 1)/(n**2 - x**2))
        )
    )

"""EXECUTION"""
run_meep = True
compute_power = True
do_plots = True
compute_modes_coeff = True
set_manual_inc_angle = False


""" PARAMETERS """
alpha_deg = 45                          # triangle angle (base - side) degrees
n = 1.50                                # index waveguide
dpml = 2                                # thickness of PML
pad = 0.375                             # pad between prism and waveguide
prism_length = 30                       # length of prism
offsx = -1.2*prism_length/2             # offset of prism from the center of the cell (== 0 --> left vertex of the prism in the center of the cell)
offsy = -prism_length/5                 # offset y-axis
offs_deg = +2.25                        # offset in degrees from critical angle
df_factor = 0.3							# fwidth of source
theta_inc = 45                          # if set_manual_inc_angle theta_inc is set


def main():
	fcen = 1
	wg_width = 0.9/(2*fcen*n)
	n_p = n + 0.29							# index prism
	n_eff = brentq(f, 1.01, n-0.01, args=(n, wg_width))
	df = fcen*df_factor

	alpha = math.radians(alpha_deg)         # triangle angle
	if not set_manual_inc_angle:
		theta_inc = math.asin( math.sin(math.asin(n_eff/n_p) - alpha) * n_p) + alpha
		theta_inc += math.radians(offs_deg)

	k = mp.Vector3(1).rotate(mp.Vector3(0, 0, -1), theta_inc)

	sx = prism_length + 2*dpml                              # cell size x-axis
	sy = (int)(1.5 * prism_length * math.atan(alpha))       # cell size y-axis (with 1.5 scale margin)
	prism = create_ideal_prism(alpha_deg, n_p, mp.Vector3(offsx, offsy), sx, sy)

	wg_y = -pad - wg_width/2 + offsy
	wg = create_h_waveguide(wg_y, wg_width, n)

	src_size = mp.Vector3(y=sy/6 - offsy  if offsx > -sx/2 
				else (sy/2 - offsy - dpml - (-offsx - sx/2 + dpml - offsy)*math.tan(alpha))*0.65)
	src_center = mp.Vector3(-sx/2 + dpml, sy/3 - src_size.y/2 - dpml)

	src = [mp.GaussianBeam2DSource(
		src=mp.GaussianSource(fcen, fwidth=df),
		center=src_center,
		size=src_size,
		beam_x0=src_center + sx*k/4,                                 # relatiu al centre de la font
		beam_kdir=k,
		beam_w0=40,                                      # beam waist
		beam_E0=mp.Vector3(0, 0, 1),
		)]


	resolution = n_p*20*fcen
	fr_pad = 0.5										# frame region pad
	nfreq = 200
	src_fr = mp.FluxRegion(center=src_center + mp.Vector3(x = fr_pad), size=src_size*1.8)


	"""RESET to get Source Flux"""
	if (compute_power or compute_modes_coeff) and run_meep:
		sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
							sources=src,
							resolution=resolution,
							boundary_layers=[mp.PML(dpml)]
							)
			
		src_region = sim.add_flux(fcen, df, nfreq, src_fr)

		sim.run(mp.at_beginning(mp.output_epsilon),
				until_after_sources=2)
		
		src_flux = mp.get_fluxes(src_region)
		sim.reset_meep()


	sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
				geometry=[prism, wg],
				sources=src,
				resolution=resolution,
				boundary_layers=[mp.PML(dpml)]
				)

	wg_mult = 4

	wg_fr = mp.FluxRegion(center=mp.Vector3(sx/2 - dpml - fr_pad, wg_y), size=mp.Vector3(y=wg_mult*wg_width))
	refl_fr = [
				mp.FluxRegion(
					center=mp.Vector3(sx/2-dpml - fr_pad, sy/4 + (offsy - dpml)/2 + wg_mult*wg_width/2),
					size=mp.Vector3(y=sy/2-offsy - 2*fr_pad - dpml - wg_mult*wg_width)
					), 
				mp.FluxRegion(
					center=mp.Vector3(sx/4 - dpml - fr_pad, sy/2 - fr_pad - dpml),
					size=mp.Vector3(x=sx/2)
				)]

	if compute_power:
		wg_region = sim.add_flux(fcen, df, nfreq, wg_fr)
		refl_region = sim.add_flux(fcen, df, nfreq, *refl_fr)

	elif compute_modes_coeff:
		wg_region = sim.add_flux(fcen, df, nfreq, wg_fr)

	if run_meep:
	
		sim.run(mp.at_beginning(mp.output_epsilon),
				mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)), 
				until=mp.stop_when_energy_decayed(dt=int(1/df), decay_by=1e-6))

		
		if compute_power:
			wg_flux = mp.get_fluxes(wg_region)
			refl_flux = mp.get_fluxes(refl_region)
			freqs = mp.get_flux_freqs(refl_region)

			if compute_modes_coeff: 
				res = sim.get_eigenmode_coefficients(wg_region, bands=[1])
				coeff = res.alpha
				forward = np.abs(coeff[0, :, 0])**2

				#return np.sum(forward)/np.sum(src_flux)
			

			if mp.am_master():
				np.savez(
					f"DIFF_TMP-off{offs_deg}_fcen{fcen:.2f}_w{wg_width:.2f}_al{alpha_deg}_n{n}_pad{pad}_df{df_factor}.npz",
					freqs=freqs,
					src_flux=src_flux,
					wg_flux=wg_flux,
					refl_flux=refl_flux,
					resolution=resolution,
					fcen=fcen,
					df=df,
					prism_length=prism_length
				)

				if do_plots:
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

					if compute_modes_coeff:
						plt.figure()
						plt.plot(freqs,forward)
						plt.xlabel(r'frequency $f (kHz)$')
						plt.ylabel("Power band 1")
						plt.title("Energy for fundamental mode")
						plt.grid(True, which="both", alpha=0.3)
						plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

					plt.show()
		
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
		
					plt.show()

		elif compute_modes_coeff: 
				res = sim.get_eigenmode_coefficients(wg_region, bands=[1])
				coeff = res.alpha
				forward = np.abs(coeff[0, :, 0])**2

				return np.sum(forward)/np.sum(src_flux)
			

	if not run_meep and do_plots:
		sim.plot2D(fields=mp.Ez,
			eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
			field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
			boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
			output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)))

		plt.show()


"""
offsets = np.linspace(-5,5,41)
[np.float64(0.030349811445287814), np.float64(0.03328358245158583), np.float64(0.036304559702988085), np.float64(0.03939857078938648), np.float64(0.04255066705083255), np.float64(0.045745262088528464), np.float64(0.04896627469403276), np.float64(0.05219727416030071), np.float64(0.05542162599696516), np.float64(0.05862263617034765), np.float64(0.06178369212108903), np.float64(0.06488839897278853), np.float64(0.06792070952709553), np.float64(0.07086504666010171), np.float64(0.07370641813335992), np.float64(0.07643052066393524), np.float64(0.07902383542593608), np.float64(0.08147371298155252), np.float64(0.08376844767000557), np.float64(0.08589734153733151), np.float64(0.0878507578255884), np.float64(0.08962016421345323), np.float64(0.0911981660972766), np.float64(0.0925785302764318), np.float64(0.09375619946522658), np.float64(0.09472729808467614), np.float64(0.09548912981579522), np.float64(0.09604016739744978), np.float64(0.09638003515102275), np.float64(0.09650948470187334), np.float64(0.09643036435078371), np.float64(0.09614558255445406), np.float64(0.0956590657520411), np.float64(0.09497571146546374), np.float64(0.09410133617040803), np.float64(0.09304261919641589), np.float64(0.09180704249577411), np.float64(0.09040282683123879), np.float64(0.08883886474147326), np.float64(0.08712465068960408), np.float64(0.08527020877601442)]

pads = np.linspace(0.3725,0.3775,5)
[np.float64(0.10110430134250178), np.float64(0.10138044307438494), np.float64(0.10160857904883383), np.float64(0.10124620533198485), np.float64(0.10093752621688946)]

prism_factor = np.linspace(0.25,0.35,6)
[np.float64(0.059625347816201595), np.float64(0.0874101160235997), np.float64(0.09957728003021563), np.float64(0.09450362750670986), np.float64(0.07466947879140647), np.float64(0.04844309278009033)]



"""
if __name__ == "__main__":
	main()

"""
	is_ok = True

	try:
		eff = []
		widths = np.linspace(0.8,1.2,9)
		for o in widths:
			print("EXECUTING MEEP WITH WG WIDTH = ", o)
			eff.append(main(o))
	
	except Exception:
		is_ok = False	

	print(eff)

	if is_ok and mp.am_master():
		plt.figure()
		plt.plot(widths, eff)
		plt.xlabel(r'offs _deg$')
		plt.title("Waveguide Efficiency")
		plt.grid(True, which="both", alpha=0.3)
		plt.show()
"""