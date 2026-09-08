from grating_coupling import GratingCoupler
from packages.figures import Layer
from packages.utils import *
import meep as mp
#import numpy as np
import os
import time
import math

"""
pid = 29341  # PID of the process you are waiting for

while True:
    try:
        os.kill(pid, 0)  # Check whether the process still exists
        time.sleep(10)
    except ProcessLookupError:
        # Process has ended
        break
time.sleep(60)
"""

gc = GratingCoupler()
vars = gc.get_variables()           # ['theta_deg', 'n_cells', 'grat_period', 'grat_height', 'grat_depth_factor', 'grat_duty_cycle']

gc.grating_excaved = False
gc.set_global_param(vars[0], 45)
gc.set_global_param(vars[1], 20)
gc.set_global_param(vars[2], 0.41)
gc.set_global_param(vars[3], 0.2)
gc.set_global_param(vars[4], 1-20e-3/gc.eff_params.grat_height)
gc.set_global_param(vars[5],0.5)

gc.src_freq = 1/0.6321

gc.n_grat = 1.7
gc.bottom_layers.append(Layer(height=1.1261, index=2.30245))
gc.bottom_layers.append(Layer(height=mp.inf, index=2.2024))

"""
theta = 34
gc.set_global_param(vars[0], theta)

period = (0.6321/(2.30245 - math.sin(math.radians(theta))))
write_output(period)
gc.set_global_param(vars[2], period)
"""

#gc.manual_res = int(25 * 3.47 * 0.975)
#gc.set_manual_resolution = True
#gc.h5_frames_per_wlength = 20
#gc.n_wlengths_power_measure = 5

gc.width_sim_scale = 1.2
gc.width_expansion = gc.Region.RIGHT
gc.compute_eff_by_modes = True
gc.n_modes_to_compute = 8
gc.res_factor = 50

gc.run_meep = False
#gc.h5_file_transistent = True
gc.do_plots = True
gc.show_region_converged_state = True

write_output(gc.main())
"""
write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output("RESOLUTION CHECKING")

res_factors, resolutions, hist_res = gc.check_simulation_convergence_changing_res()

gc.res_factor = res_factors[-1]

write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output("CELLS CHECKING")

cells, hist_res = gc.check_simulation_convergence_changing_n_cells()
"""



"""
hist_res = [0]
res = gc.main()['PR']
write_output(f"RES_FACTOR: {gc.res_factor} --> RESOLUTION = {gc.comp_resolution()}")
write_output(f"EFFICIENCY: {res}")

res_factors = [gc.res_factor]
resolutions = [gc.comp_resolution()]

while abs(hist_res[-1] - res)/res > 0.001:
    hist_res.append(res)
    gc.res_factor += 5
    res = gc.main()['PR']
    write_output("----------------------------------")
    write_output(f"RES_FACTOR: {gc.res_factor} --> RESOLUTION = {gc.comp_resolution()}")
    write_output(f"EFFICIENCY: {res}\n\n")
    res_factors.append(gc.res_factor)
    resolutions.append(gc.comp_resolution())


hist_res.append(res)
hist_res.pop(0)

write_output("----------------------------------")
write_output("\n\nFINAL RESULT")
write_output(f"EFFS: {hist_res}")
write_output(f"RES_FACTORs: {res_factors}")
write_output(f"RESOLUTIONs: {resolutions}")
"""

    

"""
thetas = compute_freqs(1, 80, f_res=1)

h = 120e-3
gc.set_global_param(vars[3], h+20e-3)
gc.set_global_param(vars[4], h/gc.eff_params.grat_height)
gc.do_offset = True
first = True

write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output(r"SIM LAB Grating n_cells = 20 crest_height = 120nm changing $\theta$ and adjusting period")

for theta in thetas:
    write_output(f"Theta: {theta}:")    
    gc.set_global_param(vars[0], theta)
    gc.set_global_param(vars[2], 0.6321/(2.30245 - math.sin(math.radians(theta))))
    res = gc.main()

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
write_output(thetas)
for key in result:
    write_output(f"{key}:")
    write_output(result[key])



thetas = compute_freqs(1, 80, f_res=1)
thetas = thetas*-1

h = 120e-3
gc.set_global_param(vars[3], h+20e-3)
gc.set_global_param(vars[4], h/gc.eff_params.grat_height)
gc.do_offset = False
first = True

write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output(r"SIM LAB Grating n_cells = 20 crest_height = 120nm changing $\theta$ and adjusting period")

for theta in thetas:
    write_output(f"Theta: {theta}:")    
    gc.set_global_param(vars[0], theta)
    gc.set_global_param(vars[2], 0.6321/(2.30245 - math.sin(math.radians(theta))))
    res = gc.main()

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
write_output(thetas)
for key in result:
    write_output(f"{key}:")
    write_output(result[key])

"""
"""
PR = result['PR']
theta_max = max(PR)
gc.set_global_param(vars[0], theta_max)


write_output("\n\n\n*******************************************************************")
write_output("******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing height - depth (general)")

parameters = gc.EffParams(theta_deg=theta_max, 
                          n_cells=20, 
                          grat_period=0.6321/(2.30245 - math.sin(theta_max)), 
                          grat_height=gc.Interval(0.05, 1.1, res_s=0.05), 
                          grat_depth_factor=gc.Interval(0.5, 1, res_s=0.05), 
                          grat_duty_cycle=0.5)

result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)
"""
"""

def compute_effs_in_indxs_range():
    indxs = compute_freqs(1, 2.2, f_res=0.05)
    first = True

    for i in indxs:
        write_output(f"INDEX: {i}:")
        gc.bottom_layers[1] = Layer(height=mp.inf, index=i)
        res = gc.main()

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
    write_output(indxs)
    for key in result:
        write_output(f"{key}:")
        write_output(result[key])

    write_output("")
    return indxs, result



write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing substrate index - crest height")


crest_heights = compute_freqs(125e-3, 350e-3, f_res=5e-3)
src_time, src_power = gc.get_src_time_power(vars[4], gc.src_freq, 140e-3)
gc.compute_src_time = True
gc.src_power = src_power
print(f"SRC_POWER: {src_power}")
effs = []

for h in crest_heights:
    write_output(f"Height: {h*1000}nm")
    gc.set_global_param(vars[3], h+20e-3)
    gc.set_global_param(vars[4], h/gc.eff_params.grat_height)
    effs.append(compute_effs_in_indxs_range())

print("FINAL RESULT")
print(effs)
"""




"""

parameters = gc.EffParams(theta_deg=45, 
                          n_cells=20, 
                          grat_period=gc.Interval(0.3, 1.6, res_s=0.1), 
                          grat_height=gc.Interval(40e-3, 170e-3, res_s=10e-3), 
                          grat_depth_factor=gc.Interval(0.65, 0.95, res_s=0.025), 
                          grat_duty_cycle=0.5)
write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing period - height - depth")
result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)
"""

"""
write_output("\n\n\n*******************************************************************")
write_output("*******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing depth - duty cycle")

parameters = gc.EffParams(theta_deg=45, 
                          n_cells=20, 
                          grat_period=0.41, 
                          grat_height=0.7,          #gc.Interval(50e-3, 1.5, res_s=50e-3) 
                          grat_depth_factor=gc.Interval(0, 1, res_s=0.025), 
                          grat_duty_cycle=gc.Interval(0, 1, res_s=0.025))

result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)

"""


             #------------------------------------------------------------------------------------------------------------------
write_output("\n\n\n*******************************************************************")
write_output("******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing height - depth (precise)")

parameters = gc.EffParams(theta_deg=45, 
                          n_cells=20, 
                          grat_period=0.41, 
                          grat_height=gc.Interval(70e-3, 170e-3, res_s=5e-3), 
                          grat_depth_factor=gc.Interval(0.65, 0.95, res_s=0.0125), 
                          grat_duty_cycle=0.5)

result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)


write_output("\n\n\n*******************************************************************")
write_output("******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing height - depth (general)")

parameters = gc.EffParams(theta_deg=45, 
                          n_cells=20, 
                          grat_period=0.41, 
                          grat_height=gc.Interval(0.1, 1.1, res_s=0.05), 
                          grat_depth_factor=gc.Interval(0, 1, res_s=0.05), 
                          grat_duty_cycle=0.5)

result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)


"""

write_output("\n\n\n*******************************************************************")
write_output("******************************************************************")
write_output("SIM LAB Grating n_cells = 20 changing period - height")

parameters = gc.EffParams(theta_deg=45, 
                          n_cells=20, 
                          grat_period=gc.Interval(0.3, 1.3, res_s=0.025), 
                          grat_height=gc.Interval(0.2, 1.2, res_s=0.025), 
                          grat_depth_factor=0.5, 
                          grat_duty_cycle=0.5)

result, steps = gc.scan_all_parameters(gc.src_freq, parameters)

write_output("\n\n")
write_output(result.dims)
write_output(result.data)
write_output("\n\n")
write_output(steps.data)"""
