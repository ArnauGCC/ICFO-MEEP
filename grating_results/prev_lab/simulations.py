from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import *
import meep as mp


gc = GratingCoupler()
gc.n_grat = 3.47


vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']
gc.set_global_param(vars[0], 25)
gc.set_global_param(vars[1], 12)
gc.set_global_param(vars[2], 0.5)
gc.set_global_param(vars[3], 0.2521613832853026)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)


gc.bottom_layers.append(Layer(height=2.27, index=1.44))
gc.bottom_layers.append(Layer(height=mp.inf, index=3.47))

gc.n_modes_to_compute = 6
gc.width_sim_scale = 1.5

gc.src_freq = 0.875

gc.h5_file_transistent = True
gc.h5_frames_per_wlength=20
gc.do_plots = True
gc.show_region_converged_state = True

#gc.src_time = 243/4
#gc.compute_src_power = False
#gc.src_power = 951.6590633073919
#gc.compute_src_power = False

print(gc.main())


"""
write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output("SIM N_WG=3.47 n_ext = 1 theta = 25º n_cells = 12 wgw = 1.75 freq = 0.875 [Layer Stack, 2 MODES]")


parameters = gc.EffParams(25, 
                          12, 
                          0.5, 
                          1.75, 
                          gc.Interval(0, 1, res_s=0.025), 
                          gc.Interval(0, 1, res_s=0.025))


result, steps = gc.scan_all_parameters(0.875, parameters)
    

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)
"""