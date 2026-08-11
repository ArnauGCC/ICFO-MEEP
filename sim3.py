from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import write_output
import numpy as np
import meep as mp


gc = GratingCoupler()
vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period_factor', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']

gc.set_global_param(vars[0], 10)
gc.set_global_param(vars[1], 15)
gc.set_global_param(vars[2], 0.5)
gc.set_global_param(vars[3], 1.5)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)

# Optimal values for src_freq = freq (= 1)
#vals = [np.float64(10), np.int64(15), np.float64(0.5215107960457857), np.float64(1.3727272727272728), np.float64(0.430952380952381), np.float64(0.447215536902577)]
#for i in range(len(vals)):
#    gc.set_global_param(vars[i], vals[i])

gc.n_wg = 3.47
gc.res_factor = 20
gc.bottom_layers.append(Layer(2.27, 1.44))
gc.bottom_layers.append(Layer(mp.inf, 3.47))
gc.n_default = 1.46


#freqs, effs = gc.get_eff_in_freq_range(0.8, 2, f_res=0.1)
#print(freqs)
#print(effs)

gc.h5_file_transistent = True
gc.compute_only_1st_mode = True
gc.src_freq = 1.48
#gc.end_src = True
print(gc.main())


#write_output("SIMULATION 3")
#write_output("")
#
#
#gc.optimize_all_parameters(1)





"""
gc.set_global_param(vars[0], 16.12142857142857)

ef1, fact1 = gc.find_max_efficiency(vars[1], 0.5, 2, 16, 1)
gc.set_global_param(vars[1], fact1)

ef2, fact2 = gc.find_max_efficiency(vars[2], 0.5, 1.5, 11, 1)
gc.set_global_param(vars[2], fact2)

ef3, fact3 = gc.find_max_efficiency(vars[3], 0.2, 0.8, 7, 1)
gc.set_global_param(vars[3], fact3)

ef4, fact4 = gc.find_max_efficiency(vars[4], 0.2, 0.8, 7, 1)

write_output(vars[1:])
write_output(ef1, ef2, ef3, ef4)
write_output(fact1, fact2, fact3, fact4)


"""

