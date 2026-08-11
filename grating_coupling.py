import meep as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from packages.figures import *
from packages.utils import *

class GratingCoupler:

	@dataclass
	class Interval:
		min_s: 		float
		max_s: 		float
		n_steps: 	float

	""" EXECUTION (Parameters for main() function) """
	run_meep = True							# Runs the simulation, False to check geometry
	do_plots = False						# True:	Plots the geometry (and fields if run_meep)
	end_src = False							# Computes efficencies turning on and off the source
	compute_eff = True						# Computes the efficiency of the coupled waveguide
	compute_only_1st_mode = False			# Computes the efficiency of the first mode of the waveguide
	compute_src_time = True					# True: src_time is computed in main(), else: src_time is set manual
	compute_src_power = True				# True: main() computes the source power, else: src_power is set manual
	show_region_converged_state = False		# Shows the region used to determinate a converged state (if do_plots)
	h5_file_transistent = False				# Creates an h5 file with the transistent state of the simulation


	""" PARAMETERS """
	freq = 1
	eff_params = {							# Variables that can be used find max efficiency 
		"theta_deg":            16, 		# Degrees of inclination from normal incidence to the waveguide
		"n_cells":				25,			# Number of cells of the grating
		"grat_period_factor":   1, 			# The grating period is gp = eff_params['grat_period_factor'] / freq
		"wg_width_factor":      1,			# The waveguide widyh is wg_width = eff_params["wg_width_factor"] / (2*freq)
		"grat_height_factor":   1/2,		# The grating height is gh = eff_params['grat_height_factor'] * wg_width
		"grat_duty_cycle":		0.5,		# Sets the grating duty cycle
	}

	n_wg = 1.5								# waveguide refraction index
	n_wlengths_power_measure = 30			# Number of wavelengths spent to measure the efficiency: self.compute_power or self.compute_modes_coeff
	bottom_layers: list[Layer]=[]			# Layers that can be added at the bottom of the waveguide
	n_default = 1							# Default refraction index for the simulation
	src_freq = freq							# Frequency of the source (in units of 1/um)

	res_factor = 20							# Number of pixels for wavelength in the highest refraction index 
	src_time = -1							# if not compute_src_time, is assumed a converged state after src_time time units to compute efficiency
	src_power = -1							# if not compute_src_power, then must specify src_flux to compute efficiency
	h5_frames_per_wlength = 4				# Number of frames per wavelength in the h5 files
	h5_name = ""


	def n_max(self):
		if len(self.bottom_layers) == 0: return self.n_wg
		return max(max(l.index for l in self.bottom_layers), self.n_wg)


	def compute_initial_parameters(self, freq, gp, wg_width, gh, gdc):
		pad_src_wg = 5/freq
		pad_inf = 3/freq
		dpml = 1/freq


		sx = int(self.eff_params["n_cells"]*gp*2.25)
		sy = wg_width+2*dpml+pad_src_wg
		for l in self.bottom_layers:
			if l.width != mp.inf:	break
			sy += l.width
		sy = int(sy+pad_inf)

		theta = math.radians(90 - self.eff_params["theta_deg"])
		wg_y = sy/2 - dpml - pad_src_wg - wg_width/2
		geometry = create_h_grating(gp, gh, gdc, self.eff_params["n_cells"], wg_width, self.n_wg, mp.Vector3(-sx/4 + pad_src_wg*0.9*math.tan(math.radians(self.eff_params["theta_deg"])), wg_y), n_ext=self.n_default)

		bottom = wg_y - wg_width/2
		for l in self.bottom_layers:
			if l.width == mp.inf:
				geometry.append(create_h_waveguide(bottom - pad_inf/2, pad_inf, l.index))
				break	
			geometry.append(create_h_waveguide(bottom - l.width/2, l.width, l.index))
			bottom -= l.width
		
		src_size = mp.Vector3(x = sx/2)
		src_center = mp.Vector3(-sx/2 + src_size.x/2, sy/2 - dpml)
		k = mp.Vector3(1).rotate(mp.Vector3(0, 0, -1), theta)
		beam_w0 = sx/2

		#top = wg_y + wg_width/2
		#geometry.append(create_h_waveguide(top + dpml/60, 2*dpml/60, 2.44))

		src = [mp.GaussianBeam2DSource(
				src=mp.ContinuousSource(self.src_freq),
				center=src_center,
				size=src_size,
				beam_x0=sy*k/4,                 # relatiu al centre de la font
				beam_kdir=k,
				beam_w0=beam_w0,                # beam waist
				beam_E0=mp.Vector3(0, 0, 1),
				)]
		
		sx += 2*dpml
		return dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src


	def main(self, freq=None, src_time=None, src_power=None):
		if freq is None:
			freq = self.freq

		if not self.compute_src_power and src_power is None:
			src_power = self.src_power

		if not self.compute_src_time and src_time is None:
			src_time = self.src_time


		wg_width = 	self.eff_params["wg_width_factor"] / (2*self.n_wg*freq)
		gp = 		self.eff_params['grat_period_factor'] / freq
		gdc = 		self.eff_params['grat_duty_cycle']
		gh = 		self.eff_params['grat_height_factor'] * wg_width

		df = 0
		nfreq =  1
		if self.h5_name != "":	self.h5_name = self.h5_name + '-'

		dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src = self.compute_initial_parameters(freq, gp, wg_width, gh, gdc)
		resolution = int(self.res_factor * self.n_max() * self.src_freq)

		src_fr = mp.FluxRegion(center=src_center - mp.Vector3(y=1), size=src_size*1.2)

		"""Compute source power to get efficiency (only in stable state)"""
		if self.compute_src_power and not self.end_src and self.run_meep:
			sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
									sources=src,
									resolution=resolution,
									boundary_layers=[mp.PML(dpml)],
									default_material=mp.Medium(index=self.n_default)
									)

			sim_src.run(until=10)

			src_region = sim_src.add_flux(self.src_freq, df, nfreq, src_fr)
			sim_src.run(until= self.n_wlengths_power_measure/self.src_freq)
			src_power = mp.get_fluxes(src_region)[0]


		sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
							geometry=geometry,
							sources=src,
							resolution=resolution,
							boundary_layers=[mp.PML(dpml)], 
							default_material=mp.Medium(index=self.n_default)
							)

		conv_region_size = mp.Vector3(sx/4, wg_width)
		conv_region_center = mp.Vector3((sx - conv_region_size.x)/2, wg_y)

		wg_fr = mp.FluxRegion(center=mp.Vector3(sx/2 - 2*dpml, wg_y), size=mp.Vector3(y=wg_width*2))
		
		if self.run_meep:
			stop_cond = make_stop_when_converged(center=conv_region_center, size=conv_region_size,
												window=int(2*10*resolution/(self.src_freq*self.n_wg)), print_err=False, 
												tolerance=1e-5, err_rate=5000)

			if self.end_src and self.compute_eff:
				wg_region = sim.add_flux(self.src_freq, df, nfreq, wg_fr)

			run_args = [mp.at_beginning(mp.with_prefix(f"{self.h5_name}", mp.output_epsilon))]
			if self.h5_file_transistent:
				run_args.append(mp.at_every(1/(self.h5_frames_per_wlength*self.src_freq), mp.to_appended(f"{self.h5_name}stp_cond-ez", mp.output_efield_z)))

			sim.run(*run_args,
					until=stop_cond if self.compute_src_time else src_time)


			if self.end_src:
				time = sim.meep_time()
				sim.change_sources([])

			if self.compute_eff:
				
				if not self.end_src:
					wg_region = sim.add_flux(self.src_freq, df, nfreq, wg_fr)

				sim.run(
#	                    mp.at_beginning(mp.with_prefix(f"{self.h5_name}", mp.output_epsilon)),
	                    mp.at_every(1/(self.h5_frames_per_wlength*self.src_freq), mp.to_appended(f"{self.h5_name}ez", mp.output_efield_z)),
						until= mp.stop_when_energy_decayed(dt=int(5/self.src_freq), decay_by=1e-4) if self.end_src 
						else self.n_wlengths_power_measure/self.src_freq)

				"""Compute source power to get efficiency (turning on and off)"""
				if self.compute_src_power and self.end_src:
					src = [mp.GaussianBeam2DSource(
								src=mp.ContinuousSource(self.src_freq, end_time=time),
								center=src_center,
								size=src_size,
								beam_x0=sy*k/4,                 # relatiu al centre de la font
								beam_kdir=k,
								beam_w0=beam_w0,                # beam waist
								beam_E0=mp.Vector3(0, 0, 1),
								)]

					sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
													sources=src,
													resolution=resolution,
													boundary_layers=[mp.PML(dpml)],
													default_material=mp.Medium(index=self.n_default)
													)
					
					src_region = sim_src.add_flux(self.src_freq, df, nfreq, src_fr)
					sim_src.run(
#						mp.at_every(1, mp.to_appended(f"src-ez", mp.output_efield_z)),
						until=mp.stop_when_energy_decayed(dt=int(5/self.src_freq), decay_by=1e-2))
					src_power = mp.get_fluxes(src_region)[0]

				if self.compute_only_1st_mode:
					"""Compute the efficiency of the first mode of the waveguide"""
					mode = sim.get_eigenmode_coefficients(wg_region, bands=[1])
					coeff = mode.alpha
					forward = np.abs(coeff[0, :, 0])**2
					eff = -forward[0]/src_power

				else:
					eff = -mp.get_fluxes(wg_region)[0]/src_power
				return eff

				

		if self.do_plots:
			if self.show_region_converged_state:
				fr = mp.FluxRegion(center=conv_region_center, size=conv_region_size)
				sim.add_flux(self.src_freq, 0, 1, fr)

			plot_params = {
				'eps_parameters': {'alpha': 0.8, 'cmap': 'binary', 'interpolation': 'none'},
				'field_parameters': {'alpha': 0.8, 'cmap': 'RdBu', 'interpolation': 'spline36'},
				'boundary_parameters': {
					'hatch': 'o', 'linewidth': 1.5, 'facecolor': 'y',
					'edgecolor': 'b', 'alpha': 0.3
				},
				'output_plane': mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy))
			}

			if self.run_meep:
				plot_params['fields'] = mp.Ez

			sim.plot2D(**plot_params)
			
			plt.show()


	def set_global_param(self, param, value):
		self.eff_params[param] = value


	def get_variables(self):
		return list(self.eff_params.keys())


	def find_max_efficiency(self, param, min_s, max_s, n_steps, freq):
		"""
		Computes the maximum efficiency of a given frequency (global parameter) between two values of a 
		variable in self.eff_parameter (to get the list of variables use: get_variables())
		Partial results are printed in the created file output.txt
		Returns the maximum efficieny with the corresponding value of the parameter  
		"""

		if param not in self.eff_params:	raise Exception(f"The variable ({param}) must be one of the "\
			"self.eff_params keys: {self.get_variables()}. To get the list use: get_variables()")

		self.run_meep = True
		self.do_plots = False
		self.end_src = False
		self.compute_eff = True
		self.compute_src_power = param in ["theta_deg",  "grat_period_factor", "n_cells"] 
		self.compute_src_time = True
		self.show_region_converged_state = False
		# The param "grat_period_factor" changes the size --> changes src flux --> src_power changes 

		src_power = None
		"""Compute Source Power"""
		if not self.compute_src_power:
			wg_width = 	self.eff_params["wg_width_factor"] / (2*self.n_wg*freq)
			gp = 		self.eff_params['grat_period_factor'] / freq
			gdc = 		self.eff_params['grat_duty_cycle']
			gh = 		self.eff_params['grat_height_factor'] * wg_width
			resolution = int(self.res_factor * self.n_max() * self.src_freq)
	
			dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src = self.compute_initial_parameters(freq, gp, wg_width, gh, gdc)

			sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
									sources=src,
									resolution=resolution,
									boundary_layers=[mp.PML(dpml)],
									default_material=mp.Medium(index=self.n_default)
									)

			sim_src.run(until=10)
			src_fr = mp.FluxRegion(center=src_center - mp.Vector3(y=1), size=src_size*1.2)

			src_region = sim_src.add_flux(self.src_freq, 0, 1, src_fr)
			sim_src.run(until= self.n_wlengths_power_measure/self.src_freq)
			src_power = mp.get_fluxes(src_region)[0]

		time = None
		if param  in ["wg_width_factor"]:
			self.compute_src_time = False

			"""Time for converged state is maximum for larger wg_widths"""
			if param == "wg_width_factor":
				self.set_global_param(param, max_s)

			sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
										geometry=geometry,
										sources=src,
										resolution=resolution,
										boundary_layers=[mp.PML(dpml)],
										default_material=mp.Medium(index=self.n_default)
										)

			conv_region_size = mp.Vector3(sx/4, wg_width)
			conv_region_center = mp.Vector3((sx - conv_region_size.x)/2, wg_y)

			stop_cond = make_stop_when_converged(center=conv_region_center, size=conv_region_size,
													window=int(2*10*resolution/(self.src_freq*self.n_wg)), tolerance=1e-5)

			sim.run(until=stop_cond)
			time = sim.meep_time()

		return find_max_efficiency_(self, param, min_s, max_s, n_steps, freq, 0, time, src_power, do_ints=param=="n_cells")


	def optimize_all_parameters(self, freq, theta_deg=None, n_cells=None, grat_period_factor=None, wg_width_factor=None, grat_height_factor=None, grat_duty_cycle=None):
		"""
		Computes the maximum efficiency of a given frequency (global parameter) between two values of a 
		variable in self.eff_parameter (to get the list of variables use: get_variables())
		Partial results are printed in the created file output.txt
		Returns the maximum efficieny with the corresponding value of the parameter  
		"""

		vars = self.get_variables()
		for i, var in enumerate(vars):
			if i == 0:
				if theta_deg is not None:
					if isinstance(theta_deg, int) or isinstance(theta_deg, float):
						self.set_global_param(var, theta_deg)
						continue
					elif isinstance(theta_deg, self.Interval):
						min_s, max_s, n_steps = theta_deg.min_s, theta_deg.max_s, theta_deg.n_steps
					else:
						raise Exception(f"theta_deg must be an int or an instance of GratingCoupler.Interval, not {type(theta_deg)}")
				else:
					min_s, max_s, n_steps = 5, 50, 16
			elif i == 1:
				if n_cells is not None:
					if isinstance(n_cells, int) or isinstance(n_cells, float):
						self.set_global_param(var, round(n_cells))
						continue
					elif isinstance(n_cells, self.Interval):
						min_s, max_s, n_steps = n_cells.min_s, n_cells.max_s, n_cells.n_steps
					else:
						raise Exception(f"n_cells must be an int or an instance of GratingCoupler.Interval, not {type(n_cells)}")
				else:
					min_s, max_s, n_steps = 10, 40, 7
			elif i == 2:
				if grat_period_factor is not None:
					if isinstance(grat_period_factor, int) or isinstance(grat_period_factor, float):
						self.set_global_param(var, grat_period_factor)
						continue
					elif isinstance(grat_period_factor, self.Interval):
						min_s, max_s, n_steps = grat_period_factor.min_s, grat_period_factor.max_s, grat_period_factor.n_steps
					else:
						raise Exception(f"grat_period_factor must be an int or an instance of GratingCoupler.Interval, not {type(grat_period_factor)}")	
				else:
					min_s, max_s, n_steps = 0.5, 2, 16
			elif i == 3:
				if wg_width_factor is not None:
					if isinstance(wg_width_factor, int) or isinstance(wg_width_factor, float):
						self.set_global_param(var, wg_width_factor)
						continue
					elif isinstance(wg_width_factor, self.Interval):
						min_s, max_s, n_steps = wg_width_factor.min_s, wg_width_factor.max_s, wg_width_factor.n_steps
					else:
						raise Exception(f"wg_width_factor must be an int or an instance of GratingCoupler.Interval, not {type(wg_width_factor)}")
				else:
					min_s, max_s, n_steps = 0.5, 1.5, 11
			elif i == 4:
				if grat_height_factor is not None:
					if isinstance(grat_height_factor, int) or isinstance(grat_height_factor, float):
						self.set_global_param(var, grat_height_factor)
						continue
					elif isinstance(grat_height_factor, self.Interval):
						min_s, max_s, n_steps = grat_height_factor.min_s, grat_height_factor.max_s, grat_height_factor.n_steps
					else:
						raise Exception(f"grat_height_factor must be an int or an instance of GratingCoupler.Interval, not {type(grat_height_factor)}")
				else:
					min_s, max_s, n_steps = 0.2, 0.8, 7
			elif i == 5:
				if grat_duty_cycle is not None:
					if isinstance(grat_duty_cycle, int) or isinstance(grat_duty_cycle, float):
						self.set_global_param(var, grat_duty_cycle)
						continue
					elif isinstance(grat_duty_cycle, self.Interval):
						min_s, max_s, n_steps = grat_duty_cycle.min_s, grat_duty_cycle.max_s, grat_duty_cycle.n_steps
					else:
						raise Exception(f"grat_duty_cycle must be an int or an instance of GratingCoupler.Interval, not {type(grat_duty_cycle)}")
				else:
					min_s, max_s, n_steps = 0.2, 0.8, 7

			self.set_global_param(var, self.find_max_efficiency(var, min_s, max_s, n_steps, freq)[1])

		write_output(vars)
		write_output([self.eff_params[v] for v in vars])


	def get_eff_in_freq_range(self, freq_min, freq_max, n_freqs=None, f_res=None):
		if n_freqs is None and f_res is None:
			raise Exception("Must specify n_freqs or f_res")
		if n_freqs is None:
			n_freqs = round((freq_max - freq_min)/f_res) + 1


		self.run_meep = True
		self.do_plots = False
		self.end_src = False
		self.compute_eff = True
		self.compute_src_time = True
		self.compute_src_power = True
		
		freqs = np.linspace(freq_min, freq_max, n_freqs)
		effs = []

		for f in freqs:
			self.src_freq = f
			effs.append(self.main())
			write_output(f"freq: {f}, eff: {effs[-1]}")
		return freqs, effs


if __name__ == "__main__":

	"""
	PARAMS_OPT:
	theta_deg:				16.12142857142857		--> 16
	grat_period_factor:		1.0009713799102498		--> 1
	wg_width_factor:		1.0045454545454546		--> 1
	grat_height_factor:		0.4877379115226338		--> 0.5
	grat_duty_cycle:		0.4928571428571429		--> 0.5
	
	PARAMS_INTERV:
	theta_deg			5, 50
	grat_period_factor: 0.5, 2
	wg_width_factor:	0.5, 1.5
	grat_height_factor:	0.2, 0.8
	grat_duty_cycle:	0.2, 0.8
	"""