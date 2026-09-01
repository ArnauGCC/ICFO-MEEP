from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import *
import meep as mp


gc = GratingCoupler()
vars = gc.get_variables()       # ['theta_deg', 'n_cells', 'grat_period', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']


theta = 45
n_cells = 25
grat_period = 0.41
gc.src_freq = 1/0.6321
grat_depth = 50e-3             # 50, 100, 120
thin_ito_layer = 20e-3          # 18 - 20
pe_layer = 1.126                # 0.8 - 1.2  (1.126)


# ito layer --> grating
gc.n_grat = 1.7
gc.grat_width = 100

gc.bottom_layers.append(Layer(height=pe_layer, index=2.30245))
gc.bottom_layers.append(Layer(height=mp.inf, index=2.2024))



gc.set_global_param(vars[0], theta)
gc.set_global_param(vars[1], n_cells)
gc.set_global_param(vars[2], grat_period)
gc.set_global_param(vars[3], grat_depth+thin_ito_layer)
gc.set_global_param(vars[4], grat_depth/(grat_depth+thin_ito_layer))
gc.set_global_param(vars[5],0.5)

gc.width_sim_scale = 1.3
gc.n_modes_to_compute = 6

parameters = gc.EffParams(theta_deg=theta, 
                          n_cells=n_cells, 
                          grat_period=grat_period, 
                          grat_height=100e-3, 
                          grat_depth_factor=gc.Interval(0, 1, res_s=0.025), 
                          grat_duty_cycle=gc.Interval(0, 1, res_s=0.025))


#grat_depths = [50e-3, 100e-3, 120e-3, 300e-3]
#for grat_depth in grat_depths:
write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output(f"SIM LAB GRATS grat_height=100e-3 n_cells={n_cells}")
#gc.h5_name = str(grat_depth)

#parameters.grat_height = grat_depth+thin_ito_layer
#parameters.grat_depth_factor = grat_depth/(grat_depth+thin_ito_layer)
result, steps = gc.scan_all_parameters(gc.src_freq, parameters)
#write_output(result)
write_output(result.dims)
write_output(result.data)
write_output("\n")
write_output(steps.data)






#gc.run_meep = False
#gc.h5_file_transistent = True
#gc.do_plots = True
#gc.show_region_converged_state = True
#print(gc.main())
