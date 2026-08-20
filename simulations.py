from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import *
import meep as mp


gc = GratingCoupler()
gc.n_wg = 3.47


vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period_factor', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']
gc.set_global_param(vars[0], 25)
gc.set_global_param(vars[1], 12)
gc.set_global_param(vars[2], 0.5)
gc.set_global_param(vars[3], 1.5)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)


gc.bottom_layers.append(Layer(2.27, 1.44))
gc.bottom_layers.append(Layer(mp.inf, 3.47))

gc.n_modes_to_compute = 6
gc.width_sim_scale = 1.5


write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output("SIM N_WG=3.47 FORCE 2 MODES n_ext = 1 theta = 25º n_cells = 12 [Layer Stack]")

result = []
for wg_w in compute_freqs(1, 2, 0.05):
    write_output("WG WIDTH: ", wg_w)
    gc.set_global_param(vars[3], wg_w)
    freqs, res = gc.get_eff_in_freq_range(0.8, 1.8, f_res=0.025)
    result.append(res)
    write_output("")

write_output(result)
