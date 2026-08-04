from grating_coupling import GratingCoupler
from packages.utils import write_output

gc = GratingCoupler()
vars = gc.get_variables()           # ['theta_deg', 'grat_period_factor', 'wg_width_factor', 'grat_height_factor', 'grat_duty_cycle']

gc.set_global_param(vars[0], 16.12142857142857)

gc.h5_name = "0.9"
gc.set_global_param(vars[1], 0.9)
ef1 = gc.main()

gc.h5_name = "0.99"
gc.set_global_param(vars[1], 0.99)
ef2 = gc.main()

gc.h5_name = "1"
gc.set_global_param(vars[1], 1.01)
ef1 = gc.main()


gc.h5_name = "1.01"
gc.set_global_param(vars[1], 1.01)
ef1 = gc.main()












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

