import meep as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from packages.figures import *
from packages.utils import *
from matplotlib.ticker import MultipleLocator
from scipy.optimize import brentq


class PrismCoupler:

	def do_check_fr(self, sim, sx, fr_pad, sy, refl_fr_v_size, wg_fr_y, wg_fr_size):
		"""
		Adds Flux Regions to check where is the division between copled_fr and refl_fr
		"""
		check = [
			mp.FluxRegion(
					center=mp.Vector3(sx/2-self.dpml - fr_pad, sy/2 - fr_pad - self.dpml - refl_fr_v_size),
					size=mp.Vector3(x=100)
					),
		
			mp.FluxRegion(center=mp.Vector3(sx/2 - self.dpml - fr_pad, wg_fr_y + wg_fr_size/2), size=mp.Vector3(x=100))
			]
		sim.add_flux(1, 0, 1, *check)


	def f(self, x, n, d):
		"""
		Function to compute n_eff (x) of the waveguide fundamental mode, assuming the 
		waveguide is symmetric (only a core with n > 1 and surrounded by air). d: wg_width
		"""
		return (
			np.pi * d * np.sqrt(n**2 - x**2)
			- np.arctan(
				np.sqrt((x**2 - 1)/(n**2 - x**2))
			)
		)


	def do_plots_(self, df, freqs, coupled_flux, refl_flux, src_flux, forward):
		"""
		Plots the graphics of the fluxes, used with gaussian source
		"""
		if self.compute_power:
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

		if self.compute_modes_coeff:
			plt.figure()
			plt.plot(freqs,forward)
			plt.xlabel(r'frequency $f (kHz)$')
			plt.ylabel("Power band 1")
			plt.title("Energy for fundamental mode")
			plt.grid(True, which="both", alpha=0.3)
			plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

		plt.show()

		if self.compute_power:
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

		if self.compute_modes_coeff:
			plt.figure()
			plt.plot(freqs,np.divide(forward, src_flux))
			plt.xlabel(r'frequency $f (kHz)$')
			plt.ylabel("Power")
			plt.title("WG_fundamental_mode/Source SPD")
			plt.grid(True, which="both", alpha=0.3)
			plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))
		
		plt.show()


	""" EXECUTION (Parameters for main() function) """
	run_meep = True							# Runs the simulation, False to check geometry
	gaussian_src = False					# True: Gaussian Source, False: Continuous Source
	compute_power = True					# Computes the power of the coupled wave important to change -> fr_division_y_factor
	get_refl_power = False					# If compute_power and get_refl_power: retrurn refl_power eff, else: return coupled_power eff
	do_plots = True							# If gaussian_src and run_meep: plots efficiency graphics, else: plots geometry
	show_fr_division = True					# If do_plots shows the division between coupled waveguide and radiated energy of wg
	compute_modes_coeff = True				# If True computes the efficiency only in the waveguide
	set_manual_inc_angle = True				# If True the incident angle is set to theta_inc_deg, else is computed an aproximation of best efficiency angle
	set_same_pad = False					# If True eff_params['pad'] = 1/(freq*6) 
	compute_src_time = True					# If True src_time is computed in main(), else: src_time is set manual
	compute_src_power = True				# If True main() computes the source power, else: src_flux is set manual


	""" PARAMETERS """
	eff_params = {							# Variables that can be used find max efficiency 
		"wg_factor":	0.9,				# Sets the wg_width = wg_factor/(2*fcen*n)
		"pad": 			0.15,				# Pad between prism and waveguidem is set if not self.set_same_pad
		"offs_deg": 	+2.50,				# Offset in degrees from critical angle
	}

	alpha_deg = 45                          # Triangle angle (base - side) degrees
	res_factor = 16							# Number of pixels for wavelength in the highest refraction index 
	n = 1.50                                # Waveguide refraction index
	prism_length = 15						# Length of prism
	dpml = prism_length/30					# Thickness of PML, self.prism_length/30 should work for any frequency
	offsx = -1.2*prism_length/2             # Offset of prism from the center of the cell (== 0 --> left vertex of the prism in the center of the cell)
	offsy = -prism_length/3                 # Offset y-axis
	df_factor = 0.3							# Fwidth [of source] = df_factor * fcen
	theta_inc_deg = 45						# If self.set_manual_inc_angle theta_inc is set to self.theta_inc_deg
	src_time = 80							# If not gaussian_src (continuous source) and not compute_src_time, is assumed a converged state after src_time time units
	src_flux= [-1]							# If not compute_src_power and compute power, then must specify src_flux
	n_wlengths_power_measure = 30			# Number of wavelengths spent to measure the efficiency: self.compute_power or self.compute_modes_coeff
	beam_w0 = prism_length/2				# Beam waist, prism_length/2 should work for any prism_length
	fr_division_y_factor = 1/20				# !!!IMPORTANT PARAMETER IF COMPUTING POWER!!! sets the division between the flux region of 
											# the reflected wave and the flux region of the coupled_wave, to visualize the division:
											# self.run_meep = False and self.do_plots = True the division must match with the simulation


	def compute_src_time_(self, sim, stop_condition):
		"""	
		Computes the time necessary to have the source turned on until the simulation 
		has converged (there are no fluctuations --> the energy is constant)
		"""
		sim.run(
			# mp.at_beginning(mp.output_epsilon),
			# mp.at_every(1, mp.to_appended("ez", mp.output_efield_z)), 
			until=stop_condition)

		src_time = sim.meep_time()
		return src_time


	def create_flux_regions(self, sx, sy, wg_y, wg_width, sim):
		"""
		Crea les Flux Regions on calcular Poyinting o energia:
			coupled_fr:	Flux que es correspon amb l'energia acoblada (canviar el paràmetre fr_division_y per ajustar
						els camps que s'agafen com a radiats o com a reflectats)
	
			wg_fr:		Flux que es transmet per la guia d'ona, s'agafa tota l'ona evanescent (fins a 3 vegades 
						l'amplada de la guia)
	
			refl_fr:	Flux que es correspon amb l'ona reflectada
		"""

		# Aprox l'ona evanescent arriba a 3 vegades l'amplada de la guia d'ona
		fr_pad = self.prism_length/100				# frame region pad

		# fr_division_y must be always lower than sy/2 
		fr_division_y = sy*self.fr_division_y_factor

		# Tros que entra de la Flux Region de la guia d'ona dins el prisma
		wg_fr_in_prism = fr_division_y - (wg_y + wg_width/2 + self.eff_params["pad"])

		# Tamany de la Flux Region reflexada vertical
		refl_fr_v_size = sy/2-self.offsy - fr_pad - self.dpml - wg_fr_in_prism

		wg_fr_size = fr_division_y - (wg_y - 3*wg_width)
		wg_fr_y = wg_y - 3*wg_width + wg_fr_size/2

		coupled_fr = mp.FluxRegion(center=mp.Vector3(sx/2 - self.dpml - fr_pad, wg_fr_y), size=mp.Vector3(y=wg_fr_size))
		wg_fr = mp.FluxRegion(center=mp.Vector3(sx/2 - self.dpml - fr_pad, wg_y), size=mp.Vector3(y=3*wg_width))
		refl_fr = [
					mp.FluxRegion(
						center=mp.Vector3(sx/2-self.dpml - fr_pad, sy/2 - fr_pad - self.dpml - refl_fr_v_size/2),
						size=mp.Vector3(y=refl_fr_v_size)
						),
					mp.FluxRegion(
						center=mp.Vector3(sx/4 - self.dpml - fr_pad, sy/2 - fr_pad - self.dpml),
						size=mp.Vector3(x=sx/2)
					)]

		# Mostra una Flux Region que marca la divisió entre la Flux Region de l'ona acoblada i de l'ona reflexada
		if (not self.run_meep and self.do_plots and self.compute_power) or self.show_fr_division: 
			self.do_check_fr(sim, sx, fr_pad, sy, refl_fr_v_size, wg_fr_y, wg_fr_size)


		return coupled_fr, wg_fr, refl_fr


	def compute_initial_parameters(self, alpha, wg_width, fcen, df, n_p):
		if self.set_manual_inc_angle:
			theta_inc = math.radians(self.theta_inc_deg)	
		else:
			n_eff = brentq(self.f, 1.01, self.n-0.01, args=(self.n, wg_width))
			theta_inc = math.asin( math.sin(math.asin(n_eff/n_p) - alpha) * n_p) + alpha
			theta_inc += math.radians(self.eff_params["offs_deg"])
		
		k = mp.Vector3(1).rotate(mp.Vector3(0, 0, -1), theta_inc)

		sx = self.prism_length + 2*self.dpml                              # cell size x-axis
		sy = int(1.75 * self.prism_length * math.atan(alpha))       	# cell size y-axis (with 1.5 scale margin)
		prism = create_ideal_prism(self.alpha_deg, n_p, mp.Vector3(self.offsx, self.offsy), sx, sy)

		wg_y = -self.eff_params["pad"] - wg_width/2 + self.offsy
		wg = create_h_waveguide(wg_y, wg_width, self.n)

		src_size = mp.Vector3(y=sy/6 - self.offsy  if self.offsx > -sx/2 
					else (sy/2 - self.offsy - self.dpml - (-self.offsx - sx/2 + self.dpml - self.offsy)*math.tan(alpha))*0.65)
		src_center = mp.Vector3(-sx/2 + self.dpml, sy/3 - src_size.y/2 - self.dpml)

		src = [mp.GaussianBeam2DSource(
			src=mp.GaussianSource(fcen, fwidth=df) if self.gaussian_src else mp.ContinuousSource(fcen),
			center=src_center,
			size=src_size,
			beam_x0=sx*k/4,					# relatiu al centre de la font
			beam_kdir=k,
			beam_w0=self.beam_w0,										# beam waist
			beam_E0=mp.Vector3(0, 0, 1),
			)]

		return sx, sy, prism, wg_y, wg, src_size, src_center, src


	def main(self, fcen, src_time=src_time, src_flux=src_flux):
		if not self.gaussian_src and self.compute_power and self.compute_modes_coeff: 
			Warning("If not gaussian_src (continuous source) and compute_power = compute_modes_coeff = True "\
		   				"returns always the power only in the waveguide. The global variable: compute_power "\
		   				"or compute_modes_coeff must be False.")

		if self.set_same_pad:		self.eff_params["pad"] = 1/(fcen*6)
		if not self.gaussian_src:	df = 0
		
		wg_width = self.eff_params["wg_factor"]/(2*fcen*self.n)
		n_p = self.n + 0.30							# index prism
		df = fcen*self.df_factor
		resolution = int(n_p*self.res_factor*fcen)

		alpha = math.radians(self.alpha_deg)         # triangle angle

		sx, sy, prism, wg_y, wg, src_size, src_center, src = self.compute_initial_parameters(alpha, wg_width, fcen, df, n_p)

		src_fr = mp.FluxRegion(center=src_center + mp.Vector3(x = 0.5), size=src_size*1.8)
		nfreq = 200 if self.gaussian_src else 1

		
		"""RESET to get Source Flux"""
		if self.run_meep and self.compute_src_power:
			sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
								sources=src,
								resolution=resolution,
								boundary_layers=[mp.PML(self.dpml)]
								)
			if self.gaussian_src:
				src_region = sim.add_flux(fcen, df, nfreq, src_fr)
				sim.run(until_after_sources=2)


			else:
				# finestra de 5 longituds d'ona
				stop_cond = make_stop_when_converged(src_center + mp.Vector3(x=sx/30),
													mp.Vector3(x=sy/4) + 1.8*src_size,
													window=int(2*5*resolution/fcen))
				sim.run(until=stop_cond)

				src_region = sim.add_flux(fcen, df, nfreq, src_fr)	
				sim.run(until=self.n_wlengths_power_measure/fcen)

			src_flux = mp.get_fluxes(src_region)
			sim.reset_meep()


		sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
					geometry=[prism, wg],
					sources=src,
					resolution=resolution,
					boundary_layers=[mp.PML(self.dpml)]
					)

		"""If not gaussian Source, get time until converged state"""
		if not self.gaussian_src and self.run_meep:
			if self.compute_src_time:
				# Region of convergence is the reflected wave (the wg might not have any fiels inside --> no convergence)
				stop_cond = make_stop_when_converged(mp.Vector3(sx/4, (sy/2 - self.offsy)/2), mp.Vector3(sx/2, (sy/2 - self.offsy)), window=int(2*2*resolution/fcen))
				self.compute_src_time_(sim, stop_cond)

			sim.run(until=src_time)


		coupled_fr, wg_fr, refl_fr = self.create_flux_regions(sx, sy, wg_y, wg_width, sim)
		

		if self.compute_power:
			coupled_region = sim.add_flux(fcen, df, nfreq, coupled_fr)
			refl_region = sim.add_flux(fcen, df, nfreq, *refl_fr)

		if self.compute_modes_coeff:
			wg_region = sim.add_flux(fcen, df, nfreq, wg_fr)


		"""RUN SIMULATION"""
		if self.run_meep:
			sim.run(#mp.at_beginning(mp.output_epsilon),
					mp.at_every(1, mp.to_appended(f"f{fcen}-ez", mp.output_efield_z)), 
					until=mp.stop_when_energy_decayed(dt=int(5/fcen), decay_by=1e-6) if self.gaussian_src 
					else self.n_wlengths_power_measure/fcen)

			if self.do_plots:
				sim.plot2D(fields=mp.Ez,
						eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
						field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
						boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
						output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)))
			
				plt.show()

			if self.gaussian_src:
				if self.compute_modes_coeff:
					freqs = mp.get_flux_freqs(wg_region)
					res = sim.get_eigenmode_coefficients(wg_region, bands=[1])
					coeff = res.alpha
					forward = np.abs(coeff[0, :, 0])**2

				if self.compute_power:
					coupled_flux = mp.get_fluxes(coupled_region)
					refl_flux = mp.get_fluxes(refl_region)
					freqs = mp.get_flux_freqs(refl_region)

				if self.compute_power and self.compute_modes_coeff:
					np.savez(
						f"DIFF_TMP-off{self.eff_params["offs_deg"]}_fcen{fcen:.2f}_w{wg_width:.2f}_\
							al{self.alpha_deg}_n{self.n}_pad{self.eff_params["pad"]}_df{self.df_factor}.npz",
						freqs=freqs,
						src_flux=src_flux,
						wg_flux_1st_mode=forward,
						coupled_flux=coupled_flux,
						refl_flux=refl_flux,
						resolution=self.resolution,
						fcen=fcen,
						df=df,
						prism_length=self.prism_length
						)


				if mp.am_master() and self.do_plots:
					self.do_plots_(df, freqs, coupled_flux, refl_flux, src_flux, forward)
					
			else: 
				if self.compute_modes_coeff:
					return mp.get_fluxes(wg_region)[0]/src_flux[0]

					res = sim.get_eigenmode_coefficients(wg_region, bands=[1])
					coeff = res.alpha
					forward = np.abs(coeff[0, :, 0])**2

					return forward[0]/src_flux[0]

				if self.compute_power:
					if self.get_refl_power:
						refl_flux = mp.get_fluxes(refl_region)
						return refl_flux[0]/src_flux[0]
					
					else:
						coupled_flux = mp.get_fluxes(coupled_region)
						return coupled_flux[0]/src_flux[0]
				

		else:
			sim.plot2D(fields=mp.Ez,
				eps_parameters={'alpha':0.8, 'cmap':'binary', 'interpolation':'none'},
				field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'spline36'},
				boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3}, 
				output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)))

			plt.show()


	
	def find_max_efficiency(self, param, param_min, param_max, n_steps, freq):
		"""
		Computes the maximum efficiency of a given frequency between two values of a parameter in eff_parameter
		Partial results are printed in the created file output.txt
		Returns the maximum efficieny with the corresponding parameter  
		"""
		if param not in self.eff_params:	raise Exception(f"The parameter ({param}) must be one of the self.eff_params keys: {list(self.eff_params.keys())}.")


		self.run_meep = True
		self.gaussian_src = False
		self.do_plots = False
		self.compute_src_time = param == "offs_deg"
		self.compute_src_power = False
		self.alpha = math.radians(self.alpha_deg)
		
		try:
			
			resolution = int((self.n+0.3)*self.res_factor*freq)
			sx, sy, prism, wg_y, wg, src_size, src_center, src = self.compute_initial_parameters(self.alpha, 
										self.eff_params["wg_factor"]/(2*freq*self.n), freq, 0, self.n+0.3)			
			src_time = -1

			"""Get Source time (if self.compute_src_time is False, means that main() function won't compute it)"""
			if not self.compute_src_time:
				sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
											geometry=[prism, wg],
											sources=src,
											resolution=resolution,
											boundary_layers=[mp.PML(self.dpml)])
				
				src_time = self.compute_src_time_(sim, make_stop_when_converged(mp.Vector3(sx/4, (sy/2 - self.offsy)/2), 
																					mp.Vector3(sx/2, (sy/2 - self.offsy)), 
																					window=int(2*2*resolution/freq)))
				sim.reset_meep()


			"""Compute Source flux (always the same)"""
			sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
								sources=src,
								resolution=resolution,
								boundary_layers=[mp.PML(self.dpml)]
								)
			
			src_fr = mp.FluxRegion(center=src_center + mp.Vector3(x = 0.5), size=src_size*1.8)			

			stop_cond = make_stop_when_converged(src_center + mp.Vector3(x=sx/30),
															mp.Vector3(x=sy/4) + 1.8*src_size,
															window=int(2*5*resolution/freq))
			sim.run(until=stop_cond)

			src_region = sim.add_flux(freq, 0, 1, src_fr)	
			sim.run(until=self.n_wlengths_power_measure/freq)


			src_flux = mp.get_fluxes(src_region)
			sim.reset_meep()

			"""Compute max efficiency"""
			return find_max_efficiency_(self, param, param_min, param_max, n_steps, freq, 0, src_time, src_flux)



		except Exception as exc:
			print("---------------------------EXCEPTION------------------------------")
			write_output(str(exc))


		write_output("--------------------------------------------------------------------")	
		

	def raise_exc_list_param_not_empty(self, param):
		raise Exception(f"If searching maximum efficiency in function of {param}, then {param}s must be empty.")

	def raise_exc_list_len_not_match(self, param, f_steps, len):
		raise Exception(f"The length of frequencies ({f_steps}) is different from the length of {param}s ({len}).")

	def check_exceptions(self, param, f_steps, wg_factors_len, pads_len, offs_degs_len):
		match param:
			case "wg_factor":
				if wg_factors_len != 0: self.raise_exc_list_param_not_empty(param)

			case "pad":
				if pads_len != 0: self.raise_exc_list_param_not_empty(param)

			case "offs_deg":
				if offs_degs_len != 0: self.raise_exc_list_param_not_empty(param)

			case _:
				raise Exception(f"The parameter ({param}) must be one of the self.eff_params keys: {list(self.eff_params.keys())}.")

		
		if wg_factors_len != 0 and  f_steps != wg_factors_len:
			self.raise_exc_list_len_not_match(param, f_steps, wg_factors_len)

		elif pads_len != 0 and  f_steps != pads_len:
			self.raise_exc_list_len_not_match(param, f_steps, pads_len)
		
		elif offs_degs_len != 0 and  f_steps != offs_degs_len:
			self.raise_exc_list_len_not_match(param, f_steps, offs_degs_len)


	
	def compute_max_eff_freq_range(self, f_min, f_max, f_res, param, param_min, param_max, param_steps,
								f_steps = None, wg_factors=[], pads=[], offs_degs=[],
								coupling=True, same_inc_angle=True):
		"""
		In a reange of frequencies, for each frequency (resolution or f_steps must be specified):
			Computes the maximum efficiency of a given frequency between two values of a parameter in eff_parameter
			Results are printed in the created file: output.txt (partial results are also printed)
		"""

		if (not same_inc_angle or param == "offs_deg") and coupling:
			raise Exception("The flux of the coupled wave can't be computed is same_inc_angle = False or " \
							"param = 'offs_deg', as the self.fr_division_y_factor changes for every angle. " \
							"This is:\n same_inc_angle = False --> coupling = False.") 
		

		self.compute_power = coupling
		self.compute_modes_coeff = not coupling
		self.set_manual_inc_angle = same_inc_angle
		self.set_same_pad = len(pads) == 0 and param != "pad"


		eff = []
		parameters = []

		if f_steps is None:	f_steps = round((f_max - f_min)/f_res) +1 if f_min < f_max else 1
		freqs = np.linspace(f_min, f_max, f_steps)

		self.check_exceptions(param, f_steps, len(wg_factors), len(pads), len(offs_degs))


		i = 0
		for freq in freqs:
			write_output("FREQUENCY = " + str(freq))

			# To compute correctly the src_time --> wg_factor must be maximum
			if param == "wg_factor":
				self.eff_params[param] = param_max

			elif len(wg_factors) != 0:
				self.eff_params["wg_factor"] = wg_factors[i]

			if len(pads) != 0:
				self.eff_params["pad"] = pads[i]

			if len(offs_degs) != 0:
				self.eff_params["offs_deg"] = offs_degs[i]
			

			try:
				
				e, p = self.find_max_efficiency(param, param_min, param_max, param_steps, freq)

				eff.append(e)
				parameters.append(p)

			except Exception as exc:
				print("---------------------------EXCEPTION------------------------------")
				write_output(str(exc))


			write_output("--------------------------------------------------------------------")	
			i += 1

		
		write_output("FINAL RESULT:")
		write_output("FREQS:")
		write_output(freqs)
		write_output("EFF:")
		write_output(eff)
		write_output(f"PARAMETERS ({param}):")
		write_output(parameters)
		write_output("")
		write_output("")


	def get_variables(self):
		return list(self.eff_params.keys())

	def set_global_param(self, param, value):
		self.eff_params[param] = value



if __name__ == "__main__":
	f_min = 0.7
	f_max = 3
	f_res = 0.1
	f_steps = round((f_max - f_min)/f_res) +1

	freqs = np.linspace(f_min, f_max, f_steps)

	wg_factors = [np.float64(0.7045454545454546), np.float64(0.7522727272727273), np.float64(0.7763636363636364), np.float64(0.8045454545454546), np.float64(0.8045454545454546), np.float64(0.8477272727272728), np.float64(0.8763636363636365), np.float64(0.8954545454545455), np.float64(0.9045454545454545), np.float64(0.9045454545454545), np.float64(0.9236363636363636), np.float64(0.9310227272727273), np.float64(0.9524999999999999), np.float64(0.9572727272727273), np.float64(0.9763636363636363), np.float64(0.9763636363636363), np.float64(0.9954545454545454), np.float64(0.9954545454545454), np.float64(0.9854545454545455), np.float64(1.0045454545454546), np.float64(1.0045454545454546), np.float64(1.0045454545454546), np.float64(1.0045454545454546), np.float64(1.0045454545454546)]
	pads = [np.float64(0.42363636363636364), np.float64(0.36363636363636365), np.float64(0.34272727272727277), np.float64(0.2854545454545455), np.float64(0.2572727272727273), np.float64(0.24938579545454548), np.float64(0.22276859504132235), np.float64(0.20454545454545459), np.float64(0.20454545454545459), np.float64(0.18545454545454548), np.float64(0.18545454545454548), np.float64(0.16636363636363638), np.float64(0.1572727272727273), np.float64(0.14902892561983472), np.float64(0.13818181818181818), np.float64(0.1525), np.float64(0.14272727272727273), np.float64(0.12363636363636366), np.float64(0.12363636363636366), np.float64(0.12450413223140497), np.float64(0.11912396694214877), np.float64(0.11365702479338843), np.float64(0.12363636363636366), np.float64(0.10454545454545455)]
	offs_degs = []

	pc = PrismCoupler()
	params = list(pc.eff_params.keys())		# params = ["wg_factor", "pad", "offs_deg"]


	"""
	f1 = 0.7
	f2 = 3

	indx1 = np.where(np.isclose(freqs, f1))[0][0]
	indx2 = np.where(np.isclose(freqs, f2))[0][0]

#	print(indx2, indx1)

	compute_max_eff_freq_range(	f1, f2, 0.1, params[2], -10, 10, 11,
								f_steps=indx2 - indx1 + 1,
								wg_factors=wg_factors[indx1:indx2+1],
								pads=pads[indx1:indx2+1],
#								offs_degs=pad_wg[indx1:indx2+1],
								coupling=False, same_inc_angle=False
								)
"""
	pc.compute_max_eff_freq_range()

	freq = 0.8
	indx = np.where(np.isclose(freqs, freq))[0][0]
	pc.eff_params["wg_factor"] = wg_factors[indx]
	pc.eff_params["pad"] = pads[indx]

	print(pc.main(freq))
