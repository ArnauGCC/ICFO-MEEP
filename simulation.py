import meep as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from packages.figures import *
from matplotlib.ticker import MultipleLocator
from scipy.optimize import brentq


"""
Adds Flux Regions to check where is the division between copled_fr and refl_fr
"""
def do_check_fr(sim, sx, fr_pad, sy, refl_fr_v_size, wg_fr_y, wg_fr_size):
	check = [
		mp.FluxRegion(
				center=mp.Vector3(sx/2-dpml - fr_pad, sy/2 - fr_pad - dpml - refl_fr_v_size),
				size=mp.Vector3(x=100)
				),
	
		mp.FluxRegion(center=mp.Vector3(sx/2 - dpml - fr_pad, wg_fr_y + wg_fr_size/2), size=mp.Vector3(x=100))
		]
	sim.add_flux(1, 0, 1, *check)


"""
Function to compute n_eff of the waveguide fundamental mode
"""
def f(x, n, d):
    return (
        np.pi * d * np.sqrt(n**2 - x**2)
        - np.arctan(
            np.sqrt((x**2 - 1)/(n**2 - x**2))
        )
    )


"""
Plots the graphics of the fluxes
"""
def do_plots_(df, freqs, coupled_flux, refl_flux, src_flux, forward):
	if compute_power:
		plt.figure()
		plt.plot(freqs,coupled_flux)
		plt.xlabel(r'frequency $f (kHz)$')
		plt.ylabel("Flux")
		plt.title("Coupled SPD")
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

	if compute_power:
		plt.figure()
		plt.plot(freqs,np.divide(refl_flux, src_flux))
		plt.xlabel(r'frequency $f (kHz)$')
		plt.ylabel("Flux")
		plt.title("Reflected/Source SPD")
		plt.grid(True, which="both", alpha=0.3)
		plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

		plt.figure()
		plt.plot(freqs,np.divide(coupled_flux, src_flux))
		plt.xlabel(r'frequency $f (kHz)$')
		plt.ylabel("Flux")
		plt.title("Coupled/Source SPD")
		plt.grid(True, which="both", alpha=0.3)
		plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

	if compute_modes_coeff:
		plt.figure()
		plt.plot(freqs,np.divide(forward, src_flux))
		plt.xlabel(r'frequency $f (kHz)$')
		plt.ylabel("Power")
		plt.title("WG_fundamental_mode/Source SPD")
		plt.grid(True, which="both", alpha=0.3)
		plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))
	
	plt.show()


"""
Computes if the region given has reached a stable (converged) state, 
this is, there are no significant fluctuations of the energy inside the region 
in some time (window)

There must be some fields to reach a converged state
"""
def make_stop_when_converged(center, size,	# Energy object
							 window=10,
							 tolerance=1e-6,
							 min_field=1e-6):

	energy_history = []

	def stop_when_field_const(sim):
		Ez = sim.get_array(center=center, size=size, component=mp.Ez)
		energy = np.sum(np.abs(Ez)**2)

		energy_history.append(energy)

		if len(energy_history) < 2*window:
			return False

		avg1 = np.mean(energy_history[-2*window:-window], axis=0)
		avg2 = np.mean(energy_history[-window:], axis=0)

		# Do not check convergence before energy has arrived
		if np.max(np.abs(avg2)) < min_field:
			return False

		error = np.mean(np.abs(avg2 - avg1) / (np.abs(avg1) + 1e-20))
#		if i%10==0:	print("ERROR: ", error)
			
		return error < tolerance

	return stop_when_field_const


"""EXECUTION"""
run_meep = True
gaussian_src = False
compute_power = True
do_plots = False
compute_modes_coeff = False
set_manual_inc_angle = True
set_same_pad = True
compute_src_time = False
compute_src_power = False


""" PARAMETERS """
alpha_deg = 45                          # triangle angle (base - side) degrees
res_factor = 16							# Number of pixels for wavelength in the highest refraction index 
n = 1.50                                # index waveguide
pad = 3.00                              # pad between prism and waveguidem is set if not set_same_pad
prism_length = 15						# length of prism
dpml = prism_length/30					# thickness of PML, prism_length/30 should work for any frequency
offsx = -1.2*prism_length/2             # offset of prism from the center of the cell (== 0 --> left vertex of the prism in the center of the cell)
offsy = -prism_length/3                 # offset y-axis
offs_deg = +0.00                        # offset in degrees from critical angle
df_factor = 0.3							# fwidth of source
theta_inc_deg = 45						# if set_manual_inc_angle theta_inc is set to theta_inc_deg
src_time = 80							# if not gaussian_src and not compute_src_time, continuous source stops after src_time time units
src_flux= [-1]							# if not compute_src_power and compute power, then must specify src_flux
beam_w0 = prism_length/2				# Beam waist, prism_length/2 should work for any prism_length
fr_division_y_factor = 1/24				# !!!IMPORTANT PARAMETER IF COMPUTING POWER!!! sets the division between the flux region of 
										# the reflected wave and the flux region of the coupled_wave, to visualize the division:
										# run_meep = False and do_plots = True the division must match with the simulation


"""	
Computes the time necessary to have the source turned on until the simulation 
has converged (there are no fluctuations --> the energy is constant)
"""
def compute_src_time_(sim, stop_condition):
	sim.run(
		# mp.at_beginning(mp.output_epsilon),
		# mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)), 
		until=stop_condition)

	src_time = sim.meep_time()
	sim.reset_meep()
	return src_time


"""
Crea les Flux Regions on calcular Poyinting o energia:
	coupled_fr:	Flux que es correspon amb l'energia acoblada (canviar el paràmetre fr_division_y per ajustar
				els camps que s'agafen com a radiats o com a reflectats)

	wg_fr:		Flux que es transmet per la guia d'ona, s'agafa tota l'ona evanescent (fins a 3 vegades 
				l'amplada de la guia)

	refl_fr:	Flux que es correspon amb l'ona reflectada
"""
def create_flux_regions(sx, sy, wg_y, wg_width, sim):
	# Aprox l'ona evanescent arriba a 3 vegades l'amplada de la guia d'ona
	fr_pad = prism_length/100				# frame region pad

	# fr_division_y must be always lower than sy/2 
	fr_division_y = sy*fr_division_y_factor

	# Tros que entra de la Flux Region de la guia d'ona dins el prisma
	wg_fr_in_prism = fr_division_y - (wg_y + wg_width/2 + pad)

	# Tamany de la Flux Region reflexada vertical
	refl_fr_v_size = sy/2-offsy - fr_pad - dpml - wg_fr_in_prism

	wg_fr_size = fr_division_y - (wg_y - 3*wg_width)
	wg_fr_y = wg_y - 3*wg_width + wg_fr_size/2

	coupled_fr = mp.FluxRegion(center=mp.Vector3(sx/2 - dpml - fr_pad, wg_fr_y), size=mp.Vector3(y=wg_fr_size))
	wg_fr = mp.FluxRegion(center=mp.Vector3(sx/2 - dpml - fr_pad, wg_y), size=mp.Vector3(y=3*wg_width))
	refl_fr = [
				mp.FluxRegion(
					center=mp.Vector3(sx/2-dpml - fr_pad, sy/2 - fr_pad - dpml - refl_fr_v_size/2),
					size=mp.Vector3(y=refl_fr_v_size)
					),
				mp.FluxRegion(
					center=mp.Vector3(sx/4 - dpml - fr_pad, sy/2 - fr_pad - dpml),
					size=mp.Vector3(x=sx/2)
				)]

	# Mostra una Flux Region que marca la divisió entre la Flux Region de l'ona acoblada i de l'ona reflexada
	if not run_meep and do_plots and compute_power: 
		do_check_fr(sim, sx, fr_pad, sy, refl_fr_v_size, wg_fr_y, wg_fr_size)


	return coupled_fr, wg_fr, refl_fr


def compute_initial_parameters(alpha, wg_width, fcen, df, n_p):
	if set_manual_inc_angle:
		theta_inc = math.radians(theta_inc_deg)	
	else:
		n_eff = brentq(f, 1.01, n-0.01, args=(n, wg_width))
		theta_inc = math.asin( math.sin(math.asin(n_eff/n_p) - alpha) * n_p) + alpha
		theta_inc += math.radians(offs_deg)
	
	k = mp.Vector3(1).rotate(mp.Vector3(0, 0, -1), theta_inc)

	sx = prism_length + 2*dpml                              # cell size x-axis
	sy = int(1.5 * prism_length * math.atan(alpha))       	# cell size y-axis (with 1.5 scale margin)
	prism = create_ideal_prism(alpha_deg, n_p, mp.Vector3(offsx, offsy), sx, sy)

	wg_y = -pad - wg_width/2 + offsy
	wg = create_h_waveguide(wg_y, wg_width, n)

	src_size = mp.Vector3(y=sy/6 - offsy  if offsx > -sx/2 
				else (sy/2 - offsy - dpml - (-offsx - sx/2 + dpml - offsy)*math.tan(alpha))*0.65)
	src_center = mp.Vector3(-sx/2 + dpml, sy/3 - src_size.y/2 - dpml)

	src = [mp.GaussianBeam2DSource(
		src=mp.GaussianSource(fcen, fwidth=df) if gaussian_src else mp.ContinuousSource(fcen),
		center=src_center,
		size=src_size,
		beam_x0=src_center + sx*k/4,					# relatiu al centre de la font
		beam_kdir=k,
		beam_w0=beam_w0,										# beam waist
		beam_E0=mp.Vector3(0, 0, 1),
		)]

	return k, sx, sy, prism, wg_y, wg, src_size, src_center, src


def main(fcen, wg_factor, src_time=src_time, src_flux=src_flux):
	#fcen = 0.4
	if set_same_pad:	
		global pad
		pad = 1/(fcen*6)

	wg_width = wg_factor/(2*fcen*n)
	n_p = n + 0.30							# index prism
	df = fcen*df_factor
	resolution = int(n_p*res_factor*fcen)

	alpha = math.radians(alpha_deg)         # triangle angle

	k, sx, sy, prism, wg_y, wg, src_size, src_center, src = compute_initial_parameters(alpha, wg_width, fcen, df, n_p)

	
	"""If not gaussian Source, get time until converged state"""
	if not gaussian_src and run_meep:
		if compute_src_time:
			sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
							geometry=[prism, wg],
							sources=src,
							resolution=resolution,
							boundary_layers=[mp.PML(dpml)])

			# Region of convergence is the reflected wave (the wg might not have any fiels inside --> no convergence)
			src_time = compute_src_time_(sim,  make_stop_when_converged(	mp.Vector3(sx/4, (sy/2 - offsy)/2), 
																		mp.Vector3(sx/2, (sy/2 - offsy)), 
																		window=int(2*2*resolution/fcen)))

		src = [mp.GaussianBeam2DSource(
				src=mp.ContinuousSource(fcen, end_time=src_time),
				center=src_center,
				size=src_size,
				beam_x0=src_center + sx*k/4,                     # relatiu al centre de la font
				beam_kdir=k,
				beam_w0=beam_w0,                                      # beam waist
				beam_E0=mp.Vector3(0, 0, 1),
				)]


	src_fr = mp.FluxRegion(center=src_center + mp.Vector3(x = 0.5), size=src_size*1.8)
	nfreq = 200 if gaussian_src else 1
	if not gaussian_src:	df = 0

	"""RESET to get Source Flux"""
	if run_meep and compute_src_power:
		sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
							sources=src,
							resolution=resolution,
							boundary_layers=[mp.PML(dpml)]
							)
			
		src_region = sim.add_flux(fcen, df, nfreq, src_fr)

		sim.run(until_after_sources=2)

		src_flux = mp.get_fluxes(src_region)
		sim.reset_meep()


	sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
				geometry=[prism, wg],
				sources=src,
				resolution=resolution,
				boundary_layers=[mp.PML(dpml)]
				)


	coupled_fr, wg_fr, refl_fr = create_flux_regions(sx, sy, wg_y, wg_width, sim)
	

	if compute_power:
		coupled_region = sim.add_flux(fcen, df, nfreq, coupled_fr)
		refl_region = sim.add_flux(fcen, df, nfreq, *refl_fr)

	if compute_modes_coeff:
		wg_region = sim.add_flux(fcen, df, nfreq, wg_fr)


	"""RUN SIMULATION"""
	if run_meep:
		sim.run(mp.at_beginning(mp.output_epsilon),
				mp.at_every(10, mp.to_appended(f"f{fcen}-ez", mp.output_efield_z)), 
				until=mp.stop_when_energy_decayed(dt=int(5/fcen), decay_by=1e-6 if gaussian_src else 1e-2))


		if gaussian_src:

			if compute_modes_coeff:
				freqs = mp.get_flux_freqs(wg_region)
				res = sim.get_eigenmode_coefficients(wg_region, bands=[1])
				coeff = res.alpha
				forward = np.abs(coeff[0, :, 0])**2

			if compute_power:
				coupled_flux = mp.get_fluxes(coupled_region)
				refl_flux = mp.get_fluxes(refl_region)
				freqs = mp.get_flux_freqs(refl_region)

			if compute_power and compute_modes_coeff:
				np.savez(
					f"DIFF_TMP-off{offs_deg}_fcen{fcen:.2f}_w{wg_width:.2f}_al{alpha_deg}_n{n}_pad{pad}_df{df_factor}.npz",
					freqs=freqs,
					src_flux=src_flux,
					wg_flux_1st_mode=forward,
					coupled_flux=coupled_flux,
					refl_flux=refl_flux,
					resolution=resolution,
					fcen=fcen,
					df=df,
					prism_length=prism_length
				)


			if mp.am_master() and do_plots:
				do_plots_(df, freqs, coupled_flux, refl_flux, src_flux, forward)
				
		else: 
			if compute_modes_coeff:
				res = sim.get_eigenmode_coefficients(wg_region, bands=[1])
				coeff = res.alpha
				forward = np.abs(coeff[0, :, 0])**2

				return forward[0]

			if compute_power:
				coupled_flux = mp.get_fluxes(coupled_region)
				refl_flux = mp.get_fluxes(refl_region)

				return coupled_flux[0]/src_flux[0]
			

	if not run_meep or do_plots:
		sim.plot2D(fields=mp.Ez,
			eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
			field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
			boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
			output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)))

		plt.show()


"""
Writes output to the file output.txt
If file not exists is created in working directory
"""
def write_output(text):
	if mp.am_master():
		with open("output.txt", "a") as file:
			file.write(str(text) + "\n")


"""
This function doesn't have to be called
Computes the maximum efficiency of a given frequency between two wg_factors with known source time and flux 
"""
def find_max_efficiency_(min_s, max_s, n_steps, freq, stage, src_time, src_flux):
	write_output("")
	write_output("STAGE:  " + str(stage))

	eff = []
	steps = np.linspace(min_s, max_s, n_steps)
	write_output("STEPS:")
	write_output(steps)

	for s in steps:
		print("EXECUTING MEEP WITH WG_FACTOR = ", s)
		eff.append(main(freq, s, src_time, src_flux))

	write_output("EFF:")
	write_output(eff)

	max_v= max(eff)
	max_i = eff.index(max_v)

	second_v = max(n for n in eff if n != max_v)
	second_i = eff.index(second_v)

	# The second max is not a neighbour
	if abs(max_i - second_i) > 1:
		return find_max_efficiency(min_s, max_s, 2*n_steps, freq, stage+1, src_time, src_flux)

	# Differnce is lower than 1%
	if abs(max_v - second_v) < 0.01:
		return max_v, steps[max_i]

	else:
		step_size = (max_s - min_s) / n_steps

		if max_i == n_steps - 1:
			return find_max_efficiency(steps[max_i-1], steps[max_i] + step_size, n_steps, freq, stage+1, src_time, src_flux)

		elif max_i == 0:
			return find_max_efficiency(steps[0] - step_size, steps[1], n_steps, freq, stage+1, src_time, src_flux)
			
		elif max_i - second_i > 0:
			return find_max_efficiency(steps[second_i], steps[max_i] + step_size,
								n_steps if abs(max_v - second_v) > 0.015 or int(n_steps/2) <= 2 
								else int(n_steps/2), 
								freq, stage+1, src_time, src_flux)
		else:
			return find_max_efficiency(steps[max_i] - step_size, steps[second_i], 
										n_steps if abs(max_v - second_v) > 0.015 or int(n_steps/2) <= 2 
										else int(n_steps/2), 
										freq, stage+1, src_time, src_flux)
	
"""
Computes the maximum efficiency of a given frequency between two wg_factors
Partial results are printed in the created file output.txt
Returns the maximum efficieny with the corresponding wg_factor  
"""
def find_max_efficiency(wg_min, wg_max, n_steps, freq):
	global compute_src_time
	global compute_src_power
	global run_meep
	global gaussian_src
	global compute_power
	global do_plots
	global compute_modes_coeff
	global set_manual_inc_angle
	global set_same_pad

	run_meep = True
	gaussian_src = False
	compute_power = True
	do_plots = False
	compute_modes_coeff = False
	set_manual_inc_angle = True
	set_same_pad = True
	compute_src_time = False
	compute_src_power = False
	alpha = math.radians(alpha_deg)
	
	try:
		"""Get Source time (worst case --> max wg width)"""
		resolution = int((n+0.3)*res_factor*freq)
		k, sx, sy, prism, wg_y, wg, src_size, src_center, src = compute_initial_parameters(alpha, wg_max/(2*freq*n), freq, 0, n+0.3)
		sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
									geometry=[prism, wg],
									sources=src,
									resolution=resolution,
									boundary_layers=[mp.PML(dpml)])
		
		src_time = compute_src_time_(sim,  make_stop_when_converged(	mp.Vector3(sx/4, (sy/2 - offsy)/2), 
																				mp.Vector3(sx/2, (sy/2 - offsy)), 
																				window=int(2*2*resolution/freq)))
		
		src = [mp.GaussianBeam2DSource(
						src=mp.ContinuousSource(freq, end_time=src_time),
						center=src_center,
						size=src_size,
						beam_x0=src_center + sx*k/4,                     # relatiu al centre de la font
						beam_kdir=k,
						beam_w0=40,                                      # beam waist
						beam_E0=mp.Vector3(0, 0, 1),
						)]

		"""Compute Source flux (always the same)"""
		sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
									sources=src,
									resolution=resolution,
									boundary_layers=[mp.PML(dpml)]
									)
		
		src_fr = mp.FluxRegion(center=src_center + mp.Vector3(x = 0.5), size=src_size*1.8)			
		src_region = sim.add_flux(freq, 0, 1, src_fr)

		sim.run(until_after_sources=2)

		src_flux = mp.get_fluxes(src_region)

		"""Compute max efficiency"""
		return find_max_efficiency_(wg_min, wg_max, n_steps, freq, 0, src_time, src_flux)



	except Exception as exc:
		print("---------------------------EXCEPTION------------------------------")
		write_output(str(exc))


	write_output("--------------------------------------------------------------------")	
	


"""
In a reange of frequencies, for each frequency (resolution or f_steps must be specified):
	Computes the maximum efficiency of a given frequency between two wg_factors
	Results are printed in the created file: output.txt (partial results are also printed)
"""
def compute_max_eff_freq_range(f_min, f_max, f_res, wg_min=0.5, wg_max=1.5, wg_steps=8, f_steps = None):	
	eff = []
	wg_factors = []

	if f_steps is None:	f_steps = int((f_max - f_min)/f_res) +1 if f_min < f_max else 1
	freqs = np.linspace(f_min, f_max, f_steps)
	

	for freq in freqs:
		write_output("FREQUENCY = " + str(freq))

		try:
			
			e, wg_f = find_max_efficiency(wg_min, wg_max, wg_steps, freq)

			eff.append(e)
			wg_factors.append(wg_f)

		except Exception as exc:
			print("---------------------------EXCEPTION------------------------------")
			write_output(str(exc))


		write_output("--------------------------------------------------------------------")	

	
	write_output("FINAL RESULT:")
	write_output("FREQS:")
	write_output(freqs)
	write_output("EFF:")
	write_output(eff)
	write_output("WG_FACTORS:")
	write_output(wg_factors)
	

if __name__ == "__main__":
	compute_max_eff_freq_range(0.7, 0.7, 1, wg_min=0.69, wg_max=7, wg_steps=2)
