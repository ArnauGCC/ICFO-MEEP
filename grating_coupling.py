from __future__ import annotations
import meep as mp
import math
import matplotlib.pyplot as plt
import numpy as np
from packages.figures import *
from packages.utils import *


class GratingCoupler:
	class NamedMatrix:
		def __init__(self, data, dims):
			self.data = data
			self.dims = dims

		def get_dims(self):
			return self.dims

		def at(self, **kwargs):
			indices = [kwargs[dim] for dim in self.dims]
			return self.data[tuple(indices)]

	@dataclasses.dataclass
	class Interval:
		min_s: float
		max_s: float
		n_steps: int | None = None
		res_s: float | None = None

		def __init__(
			self,
			min_s: float,
			max_s: float,
			n_steps: int | None = None,
			res_s: float | None = None,
		):
			if n_steps is None and res_s is None:
				raise ValueError(
					"Either n_steps or res_s must be specified."
				)

			if n_steps is not None and res_s is not None:
				if n_steps < 2:
					raise ValueError("n_steps must be at least 2.")
				if res_s <= 0:
					raise ValueError("res_s must be greater than 0.")

				self.n_steps = n_steps
				self.res_s = res_s

			elif n_steps is not None:
				if n_steps < 2:
					raise ValueError("n_steps must be at least 2.")

				self.n_steps = n_steps
				self.res_s = (max_s - min_s) / (n_steps - 1)

			else:
				if res_s <= 0:
					raise ValueError("res_s must be greater than 0.")

				self.res_s = res_s
				self.n_steps = round((max_s - min_s) / res_s) + 1

				if self.n_steps < 2:
					self.n_steps = 2

			self.min_s = min_s
			self.max_s = max_s

		def get_interval(self) -> np.ndarray:
			return np.linspace(self.min_s, self.max_s, self.n_steps, )


	@dataclasses.dataclass
	class Mode:
		alpha:		np.array
		vgrp:		np.array
		kpoints:	list
		kdom:		list
		cscale:		np.array


	@dataclasses.dataclass
	class EffParams:						# Variables that can be used find max efficiency
		"""
		This parameters must be set to a simple type (int, float) for executing the simulation.
		They also can be set to an instance of GratingCoupler.Interval to find the maximum efficiency in a range of values.
		"""
		theta_deg:			float	| GratingCoupler.Interval | None = None		# Degrees of inclination from normal incidence to the waveguide
		n_cells:			int		| GratingCoupler.Interval | None = None		# Number of cells of the grating
		grat_period:		float	| GratingCoupler.Interval | None = None		# The grating period is gp = eff_params['grat_period'] / freq
		wg_width_factor:	float	| GratingCoupler.Interval | None = None		# The waveguide widyh is wg_width = eff_params["wg_width_factor"] / (2*freq)
		grat_height_factor:	float	| GratingCoupler.Interval | None = None		# The grating height is gh = eff_params['grat_height_factor'] * wg_width
		grat_duty_cycle:	float	| GratingCoupler.Interval | None = None		# Sets the grating duty cycle


	""" EXECUTION (Parameters for main() function) """
	run_meep = True							# Runs the simulation, False to check geometry
	do_plots = False						# True:	Plots the geometry (and fields if run_meep)
	end_src = False							# Computes efficencies turning on and off the source
	compute_eff_by_modes = True				# True: Computes the efficiency of the first n_modes_to_compute of the waveguide
	compute_eff_by_power = True				# True: Computes the efficiency with the total power
	compute_src_time = True					# True: src_time is computed in main(), else: src_time is set manual
	compute_src_power = True				# True: main() computes the source power, else: src_power is set manual
	show_region_converged_state = False		# Shows the region used to determinate a converged state (if do_plots)
	h5_file_transistent = False				# Creates an h5 file with the transistent state of the simulation
	set_manual_resolution = False			# True: Sets resolution = manual_res, useful to get smoother results in scripts with multiple simulations (changing resolution between sims makes results noisy --> only useful for low resolution of the results)
	do_simple_type_return = False			# True:  main() only returns a simple value (not a dictionary nor a list), if compute_eff_by_modes returns the efficiency of the last mode
	compute_eff = compute_eff_by_modes or compute_eff_by_power		# Computes the efficiency of the coupled waveguide

	""" PARAMETERS """
	eff_params = EffParams(16, 25, 1, 1, 1/2, 0.5)		
	#eff_params = {							 
	#	"theta_deg":            16, 		
	#	"n_cells":				25,			
	#	"grat_period":		  	1, 			
	#	"wg_width_factor":      1,			
	#	"grat_height_factor":   1/2,		
	#	"grat_duty_cycle":		0.5,		
	#}

	n_wg = 1.5								# waveguide refraction index
	n_wlengths_power_measure = 30			# Number of wavelengths spent to measure the efficiency: self.compute_power or self.compute_modes_coeff
	bottom_layers: list[Layer]=[]			# Layers that can be added at the bottom of the waveguide
	n_default = 1							# Default refraction index for the simulation
	src_freq = 1							# Frequency of the source (in units of 1/um)
	n_modes_to_compute = 1					# If compute_by_modes_not_power, then efficiency is computed by the first n_modes_to_compute

	res_factor = 20							# Number of pixels for wavelength in the highest refraction index 
	manual_res = 70							# if set_manual_resolution, resolution = manual_res
	src_time = -1							# if not compute_src_time, is assumed a converged state after src_time time units to compute efficiency
	src_power = -1							# if not compute_src_power, then must specify src_flux to compute efficiency
	width_sim_scale = 1						# Factor to make simulation wider (useful to avoid the apparition of reflection and radtiated fields in measuring regions)
	h5_frames_per_wlength = 4				# Number of frames per wavelength in the h5 files
	h5_name = ""							# output h5 files have h5_name


	def n_max(self):
		"""
		Returns the maximum refractive index in the simulation
		"""
		if len(self.bottom_layers) == 0: return self.n_wg
		return max(max(l.index for l in self.bottom_layers), self.n_wg)


	def comp_resolution(self, res_factor=res_factor, freq=src_freq):
		"""
		Computes the resolution needed to guarantee res_factor pixels in the highest refractive index of the simulation
		"""
		return int(res_factor * self.n_max() * freq)


	def get_valid_modes(self, modes_result, cscale_max=100.0, cscale_min=1e-4, indexs=None):
		"""
		Given the result of get_eigenmode_coefficients() returns the valid modes based on its 
		cscale, kpoint and vgrp (as its a 2D sim with an horitzontal wg). 
		IMPORTANT: if valid_indices of the result will be used, then indexs must be specified
		as the bands used to compute get_eigenmode_coefficients().
		"""
		if indexs is None:
			indexs = range(1, len(modes_result.vgrp)+1)
		elif len(indexs) != len(modes_result.vgrp):
			raise Exception("Indexs must have as many elements as modes have modes_result.")


		cscale = np.asarray(modes_result.cscale)
		vgrp = np.asarray(modes_result.vgrp)
		kx = np.array([kp.x for kp in modes_result.kpoints])

		# Boolean mask for valid physical modes
		valid_mask = (
			(cscale > cscale_min) & (cscale < cscale_max) & (vgrp > 0) & (kx > 0)
		)

		valid_indices = np.where(valid_mask)[0].tolist()

		# Slice all native attributes
		new_alpha = modes_result.alpha[valid_indices, :, :]
		new_vgrp = modes_result.vgrp[valid_indices]
		new_kpoints = [modes_result.kpoints[i] for i in valid_indices]
		new_kdom = [modes_result.kdom[i] for i in valid_indices]
		new_cscale = modes_result.cscale[valid_indices]

		# Re-construct a matching object with all fields preserved
		filtered_result = type(modes_result)(
			alpha=new_alpha,
			vgrp=new_vgrp,
			kpoints=new_kpoints,
			kdom=new_kdom,
			cscale=new_cscale,
		)

		valid_indices =  [indexs[i] for i in valid_indices]

		return filtered_result, valid_indices


	def compute_modes_power(self, sim, wg_region, forward_fields=True):
		"""
		Computes the modes power of a FluxRegion in a waveguide (wg_region) [on a 2D sim with an horitzontal wg]
		"""
		tolerances = [5e-1, 1e-1, 5e-2, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8]
		indxs = range(1, self.n_modes_to_compute+1)

		# IMPORTANT STEP!!! Avoids divergence of MPB solver (specially with large values of n_modes_to_compute and regions with weak, few fields)
		for t in tolerances:
			print("MPB SOLVER WITH TOLERANCE:  ", t)
			m1 = sim.get_eigenmode_coefficients(wg_region, bands=indxs, direction=mp.X, eig_tolerance=t)
			m11, indxs = self.get_valid_modes(m1, indexs=indxs)

		modes = sim.get_eigenmode_coefficients(wg_region, bands=indxs, direction=mp.X)
		modes, i1 = self.get_valid_modes(modes, indxs)
		
		coeffs = modes.alpha
		forward = np.abs(coeffs[:, :, 0 if forward_fields else 1])**2

		"""					
		print("POWER: ", mp.get_fluxes(wg_region)[0])
		print("")
		print("")
		print(modes)
		print("")
		print("")

		alpha = modes.alpha[:, 0, :]

		P_forward = np.abs(alpha[:, 0])**2
		P_backward = np.abs(alpha[:, 1])**2

		print("forward:", P_forward)
		print("backward:", P_backward)

		P_net = np.sum(P_forward - P_backward)

		print("net modal power:", P_net)
		"""

		return forward


	def compute_initial_parameters(self, gp, wg_width, gh, gdc):
		"""
		Computes the parameters used to start the simulation
		"""
		pad_src_wg = 5
		pad_inf = 3
		dpml = 1


		sx = int(self.eff_params.n_cells*gp*2.25)
		sy = wg_width+2*dpml+pad_src_wg
		for l in self.bottom_layers:
			if l.width != mp.inf:	break
			sy += l.width
		sy = int(sy+pad_inf)

		theta = math.radians(90 - self.eff_params.theta_deg)
		wg_y = sy/2 - dpml - pad_src_wg - wg_width/2
		geometry = create_h_grating(gp, gh, gdc, self.eff_params.n_cells, wg_width, self.n_wg, mp.Vector3(-sx/4 + pad_src_wg*0.9*math.tan(math.radians(self.eff_params.theta_deg)), wg_y), n_ext=self.n_default)

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

		src = [mp.GaussianBeam2DSource(
				src=mp.ContinuousSource(self.src_freq),
				center=src_center,
				size=src_size,
				beam_x0=sy*k/4,                 # relatiu al centre de la font
				beam_kdir=k,
				beam_w0=beam_w0,                # beam waist
				beam_E0=mp.Vector3(0, 0, 1),
				)]
		
		# Used to avoid the apparition of reflection and radtiated fields in measuring regions
		if self.eff_params.n_cells < 10:
			sx *= 3
		else:
			sx += 2*dpml

		sx *= abs(self.width_sim_scale)
		sx = int(sx)

		return dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src


	def create_flux_regions_and_stop_conds(self, sx, dpml, src_size, src_center, resolution, wg_y, wg_width):
		src_fr = mp.FluxRegion(center=src_center - mp.Vector3(y=0.1), size=src_size*1.2)

		src_conv_region_size = mp.Vector3(src_size.x*1.2, 0.4)
		src_conv_region_center = mp.Vector3(src_center.x, src_center.y - src_conv_region_size.y/2 - 0.1)
		src_stop_cond = make_stop_when_converged(center=src_conv_region_center, size=src_conv_region_size,
												 window=int(2*10*resolution/(self.src_freq*self.n_default)), 
												 tolerance=1e-4, min_field=50, max_time=25)


		conv_region_size = mp.Vector3(sx/4, wg_width)
		conv_region_center = mp.Vector3((sx - conv_region_size.x)/2, wg_y)
		stop_cond = make_stop_when_converged(center=conv_region_center, size=conv_region_size,
											 window=int(2*10*resolution/(self.src_freq*self.n_wg)), print_err=False, 
											 tolerance=1e-5, min_field=50, err_rate=5000, max_time=100)


		wg_fr_left = mp.FluxRegion(center=mp.Vector3(-sx/2 + 2*dpml, wg_y), size=mp.Vector3(y=wg_width*3), direction=mp.X)
		wg_fr_right = mp.FluxRegion(center=mp.Vector3(sx/2 - 2*dpml, wg_y), size=mp.Vector3(y=wg_width*3), direction=mp.X)

		
		return src_fr, src_stop_cond, stop_cond, wg_fr_left, wg_fr_right
		

	def main(self, freq=None, src_time=None, src_power=None):
		"""
		Creates and executes a simulation based on the global parameters values
		"""

		if self.do_simple_type_return and self.compute_eff_by_power == self.compute_eff_by_modes:
			raise Exception("Ambiguous behaviour, if do_simple_type_return, then self.compute_eff_by_power != self.compute_eff_by_modes")

		if freq is not None:
			self.src_freq = freq
		else:			# If no frequency is specified then it's used the global param src_freq (default=1)
			freq = self.src_freq

		if not self.compute_src_power and src_power is None:
			src_power = self.src_power

		if not self.compute_src_time and src_time is None:
			src_time = self.src_time


		wg_width = 	self.eff_params.wg_width_factor / (2*self.n_wg)		# It's assumed central frquency = 1
		gp = 		self.eff_params.grat_period
		gdc = 		self.eff_params.grat_duty_cycle
		gh = 		self.eff_params.grat_height_factor * wg_width


		dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src = self.compute_initial_parameters(gp, wg_width, gh, gdc)
		resolution = self.manual_res if self.set_manual_resolution else self.comp_resolution()
		src_fr, src_stop_cond, stop_cond, wg_fr_left, wg_fr_right = self.create_flux_regions_and_stop_conds(sx, dpml, src_size, src_center, resolution, wg_y, wg_width)
		df = 0
		nfreq =  1
		if self.h5_name != "":	self.h5_name = self.h5_name + '-'


		"""Compute source power to get efficiency (only in stable state)"""
		if self.compute_src_power and not self.end_src and self.run_meep:
			sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
									sources=src,
									resolution=resolution,
									boundary_layers=[mp.PML(dpml)],
									default_material=mp.Medium(index=self.n_default)
									)

			sim_src.run(until=src_stop_cond)

			src_region = sim_src.add_flux(freq, df, nfreq, src_fr)
			sim_src.run(until= self.n_wlengths_power_measure/freq)
			src_power = -mp.get_fluxes(src_region)[0]


		sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
							geometry=geometry,
							sources=src,
							resolution=resolution,
							boundary_layers=[mp.PML(dpml)], 
							default_material=mp.Medium(index=self.n_default)
							)

		
		if self.run_meep:
			if self.end_src and self.compute_eff:
				wg_region_left = sim.add_mode_monitor(freq, df, nfreq, wg_fr_left)
				wg_region_right = sim.add_mode_monitor(freq, df, nfreq, wg_fr_right)

			run_args = [mp.at_beginning(mp.with_prefix(f"{self.h5_name}", mp.output_epsilon))]
			if self.h5_file_transistent:
				run_args.append(mp.at_every(1/(self.h5_frames_per_wlength*freq), mp.to_appended(f"{self.h5_name}stp_cond-ez", mp.output_efield_z)))

			sim.run(*run_args,
					until=stop_cond if self.compute_src_time else src_time)


			if self.end_src:
				time = sim.meep_time()
				sim.change_sources([])

			if self.compute_eff:
				
				if not self.end_src:
					wg_region_left = sim.add_mode_monitor(freq, df, nfreq, wg_fr_left)
					wg_region_right = sim.add_mode_monitor(freq, df, nfreq, wg_fr_right)


				sim.run(
	                    mp.at_every(1/(self.h5_frames_per_wlength*freq), mp.to_appended(f"{self.h5_name}ez", mp.output_efield_z)),
						until= mp.stop_when_energy_decayed(dt=int(5/freq), decay_by=1e-4) if self.end_src 
						else self.n_wlengths_power_measure/freq)

				"""Compute source power to get efficiency (turning on and off)"""
				if self.compute_src_power and self.end_src:
					src = [mp.GaussianBeam2DSource(
								src=mp.ContinuousSource(self.src_freq, end_time=time),
								center=src_center,
								size=src_size,
								beam_x0=sy*k/4,
								beam_kdir=k,
								beam_w0=beam_w0,
								beam_E0=mp.Vector3(0, 0, 1),
								)]

					sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
													sources=src,
													resolution=resolution,
													boundary_layers=[mp.PML(dpml)],
													default_material=mp.Medium(index=self.n_default)
													)
					
					src_region = sim_src.add_flux(freq, df, nfreq, src_fr)
					sim_src.run(until=mp.stop_when_energy_decayed(dt=int(5/freq), decay_by=1e-4))
					src_power = -mp.get_fluxes(src_region)[0]

				result = {'S': 	src_power}

				if self.compute_eff_by_modes:
					"""Compute the efficiency of the first mode of the waveguide"""
					forward_L = self.compute_modes_power(sim, wg_region_left, forward_fields=False)
					forward_R = self.compute_modes_power(sim, wg_region_right)

					if self.do_simple_type_return:
						return forward_R[-1, 0]/src_power

					for mode in range(forward_L.shape[0]):
						eff = forward_L[mode, 0]/src_power
						result[f"{mode}L"] = eff

					for mode in range(forward_R.shape[0]):
						eff = forward_R[mode, 0]/src_power
						result[f"{mode}R"] = eff
						

				if self.compute_eff_by_power:
					eff = mp.get_fluxes(wg_region_right)[0]/src_power

					if self.do_simple_type_return: return eff
					result["PR"] = eff
					result["PL"] = -mp.get_fluxes(wg_region_left)[0]/src_power
				

		if self.do_plots:
			if self.show_region_converged_state:
				fr = mp.FluxRegion(center=stop_cond.center, size=stop_cond.size)
				fr_src = mp.FluxRegion(center=src_stop_cond.center, size=src_stop_cond.size)
				sim.add_flux(freq, 0, 1, fr)
				sim.add_flux(freq, 0, 1, fr_src)

			# src_region has to be added because it's in another simulation
			src_region = sim.add_flux(freq, 0, 1, src_fr)	
			if not self.run_meep:
				sim.add_flux(freq, 0, 1, wg_fr_left)
				sim.add_flux(freq, 0, 1, wg_fr_right)

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


		return result


	def set_global_param(self, param: str, value):
		"""
		Sets the value of a global parameter
		"""
		setattr(self.eff_params, param, value)


	def get_variables(self):
		"""
		Returns a list of the parameters that changes the efficiency of the simulation
		"""
		return [field.name for field in dataclasses.fields(self.EffParams)]


	def get_src_time_power(self, param, freq, max_step):
		self.src_freq = freq

		"""Compute Source Power"""
		# The param "grat_period" changes the size --> changes src flux --> src_power changes 
		if param not in ["theta_deg",  "grat_period", "n_cells"]:
			self.compute_src_power = False

			wg_width = 	self.eff_params.wg_width_factor / (2*self.n_wg)
			gp = 		self.eff_params.grat_period
			gdc = 		self.eff_params.grat_duty_cycle
			gh = 		self.eff_params.grat_height_factor * wg_width
			resolution = self.manual_res if self.set_manual_resolution else int(self.res_factor * self.n_max() * freq)
	
			dpml, sx, sy, wg_y, geometry, src_size, src_center, k, beam_w0, src = self.compute_initial_parameters(gp, wg_width, gh, gdc)
			src_fr, src_stop_cond, stop_cond, wg_fr_left, wg_fr_right = self.create_flux_regions_and_stop_conds(sx, dpml, src_size, src_center, resolution, wg_y, wg_width)

			sim_src = mp.Simulation(cell_size=mp.Vector3(sx, sy),
									sources=src,
									resolution=resolution,
									boundary_layers=[mp.PML(dpml)],
									default_material=mp.Medium(index=self.n_default)
									)
			sim_src.run(until=src_stop_cond)

			src_region = sim_src.add_flux(freq, 0, 1, src_fr)
			sim_src.run(until= self.n_wlengths_power_measure/freq)
			src_power = -mp.get_fluxes(src_region)[0]

		else:
			self.compute_src_power = True
			src_power = None


		if param  in ["wg_width_factor"]:
			self.compute_src_time = False

			"""Time for converged state is maximum for larger wg_widths"""
			if param == "wg_width_factor":
				self.set_global_param(param, max_step)

			sim = mp.Simulation(cell_size=mp.Vector3(sx, sy),
										geometry=geometry,
										sources=src,
										resolution=resolution,
										boundary_layers=[mp.PML(dpml)],
										default_material=mp.Medium(index=self.n_default)
										)

			sim.run(until=stop_cond)
			time = sim.meep_time()
		else:
			self.compute_src_time = True
			time = None

		return time, src_power


	def find_max_efficiency(self, param, min_s, max_s, n_steps, freq):
		"""
		Computes the maximum efficiency of a given frequency (global parameter) between two values of a 
		variable in self.eff_parameter (to get the list of variables use: get_variables())
		Partial results are printed in the created file output.txt
		Returns the maximum efficieny with the corresponding value of the parameter  
		"""

		if param not in self.get_variables():	raise Exception(f"The variable ({param}) must be one of the "\
			"self.eff_params keys: {self.get_variables()}. To get the list use: get_variables()")

		self.run_meep = True
		self.do_plots = False
		self.end_src = False
		self.compute_eff = True
		self.compute_eff_by_modes = False
		self.show_region_converged_state = False
		self.do_simple_type_return = True

		time, src_power = self.get_src_time_power(param, freq, max_s)

		return find_max_efficiency_(self, param, min_s, max_s, n_steps, freq, 0, time, src_power, do_ints=param=="n_cells")


	default_intervals_search = EffParams(
		Interval(5, 50, 16),		# Interval for max eff search of theta_deg
		Interval(5, 40, 8),			# Interval for max eff search of n_cells
		Interval(0.5, 2, 16),		# Interval for max eff search of grat_period
		Interval( 0.5, 1.5, 11),	# Interval for max eff search of wg_width_factor
		Interval(0.2, 0.8, 7),		# Interval for max eff search of grat_height_factor
		Interval(0.2, 0.8, 7)		# Interval for max eff search of grat_duty_cycle
	)

	def optimize_all_parameters(self, freq, opt_params: EffParams=EffParams()):
		"""
		Computes the maximum efficiency of a given frequency (global parameter) between two values of a 
		variable in self.eff_parameter (to get the list of variables use: get_variables())
		
		opt_params can be used to fix a value when optimizing the other ones or to sepcifiy manually an 
		interval of search for a parameter otherwise, self.default_intervals_search will be used.
		
		Partial results are printed in the created file output.txt
		Returns the maximum efficieny with the corresponding value of the parameter  
		"""

		for param in dataclasses.fields(opt_params):
			value = getattr(opt_params, param.name)
			if value is not None and (isinstance(value, int) or isinstance(value, float)):
				self.set_global_param(param.name, value)

		for param in dataclasses.fields(opt_params):
			value = getattr(opt_params, param.name)
			if value is None:
				interv = getattr(self.default_intervals_search, param.name)
				min_s, max_s, n_steps = interv.min_s, interv.max_s, interv.n_steps
			elif isinstance(value, self.Interval):
				min_s, max_s, n_steps = value.min_s, value.max_s, value.n_steps
			else:
				continue

			self.set_global_param(param.name, self.find_max_efficiency(param.name, min_s, max_s, n_steps, freq)[1])


		opt_values = [getattr(opt_params, param.name) for param in dataclasses.fields(opt_params)]
		write_output(self.get_variables())
		write_output(opt_values)

		return opt_values


	def scan_all_parameters(self, freq=src_freq, opt_params: EffParams=EffParams()):

		self.src_freq = freq
		fields = dataclasses.fields(opt_params)
		i = 0

		parameters = []
		for param in fields:
			value = getattr(opt_params, param.name)
			if value is not None and (isinstance(value, int) or isinstance(value, float)):
				self.set_global_param(param.name, value)
				i += 1
			else:
				parameters.append(param.name)


		if i >= len(fields):	raise Exception("There is no parameter to scan, at least one parameter of opt_params must be an Interval or None (default interval)")
		end_of_recursivity = i == len(fields) - 1


		for param in fields:
			value = getattr(opt_params, param.name)
			if value is None:
				interv = getattr(self.default_intervals_search, param.name)
				min_s, max_s, n_steps = interv.min_s, interv.max_s, interv.n_steps
			elif isinstance(value, self.Interval):
				min_s, max_s, n_steps = value.min_s, value.max_s, value.n_steps
			else:
				continue

			
			steps = np.linspace(min_s, max_s, n_steps)
			if end_of_recursivity:
				first = True
				time, power = self.get_src_time_power(param, freq, max_s)

				for step in steps:
					write_output("  "*i, f"Executing with {param.name}: ", step)
					self.set_global_param(param.name, step)
					res = self.main(src_time=time, src_power=power)

					if first:
						first = False
						result = {key: [] for key in res}
					else:
						for key in res:
							if key not in result:
								result[key] = [0]*len(result[next(iter(result))])
		
					for key in result:
						result[key].append(res[key] if key in res else np.float64(0))
						write_output("  "*i, f" {key}: {result[key][-1]}")

				write_output(result)
				return self.NamedMatrix(result, parameters), self.NamedMatrix(steps, parameters)

			else:
				result = []
				result_steps = []
				for step in steps:
					write_output("  "*i, f"Executing with {param.name}: ", step)
					setattr(opt_params, param.name, step)
					data, steps = self.scan_all_parameters(freq, opt_params)
					result.append(data.data)
					result_steps.append(steps.data)

				return self.NamedMatrix(result, parameters), self.NamedMatrix(result_steps, parameters) 




	def get_eff_in_freq_range(self, freq_min, freq_max, n_freqs=None, f_res=None, set_same_resolution=False):
		"""
		Given a frequency range (f_min, f_max) and its resolution computes the efficiency of the simulation
		in this frequencies.

		set_same_resolution is an important parameter to obtain good results with high resolution (f_res < 0.01), 
		otherwise the results won't be smooth.
		"""

		self.run_meep = True
		self.do_plots = False
		self.end_src = False
		self.compute_eff = True
		self.compute_src_time = True
		self.compute_src_power = True

		if set_same_resolution:
			self.set_manual_resolution = True
			self.manual_res = self.comp_resolution(freq=freq_max)
		
		freqs = compute_freqs(freq_min, freq_max, f_res, n_freqs)
		first = True

		for f in freqs:
			write_output(f"Freq: {f}:")
			res = self.main(freq=f)

			if first:
				first = False
				result = {key: [] for key in res}
			else:
				for key in res:
					if key not in result:
						result[key] = [0]*len(result[next(iter(result))])

			for key in result:
				result[key].append(res[key] if key in res else np.float64(0))
				write_output(f"    {key}: {result[key][-1]}")


		write_output("")
		write_output("")
		write_output(freqs)
		for key in result:
			write_output(f"{key}:")
			write_output(result[key])

		write_output("")


		return freqs, result


if __name__ == "__main__":
	pass