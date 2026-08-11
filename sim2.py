from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import write_output
import meep as mp
import numpy as np

gc = GratingCoupler()
vars = gc.get_variables()           # ['theta_deg', 'grat_period_factor', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']

"""
gc.set_global_param(vars[0], 37.325)
gc.set_global_param(vars[1], 25)
gc.set_global_param(vars[2], 1)
gc.set_global_param(vars[3], 1.5)
gc.set_global_param(vars[4], 0.5)
gc.set_global_param(vars[5],0.5)
"""

vals = [np.float64(42.9375), np.int64(21), np.float64(1.0033333333333334), np.float64(1.4572727272727273), np.float64(0.5284328546233308), np.float64(0.4928571428571429)]
for i in range(len(vals)):
    gc.set_global_param(vars[i], vals[i])


gc.n_wg = 2
gc.res_factor = 16


gc.bottom_layers.append(Layer(2.27, 1.5))
gc.bottom_layers.append(Layer(mp.inf, 2))

gc.h5_file_transistent = True

gc.run_meep = False
gc.do_plots = True
print(gc.main())


#write_output("SIMULATION 2")
#write_output("")


#gc.optimize_all_parameters(1)

#·write_output("###########################################################################")
#write_output("###########################################################################")
#write_output("")
#write_output("")




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

