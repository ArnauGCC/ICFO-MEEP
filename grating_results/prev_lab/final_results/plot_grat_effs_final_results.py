import numpy as np
from packages.utils import *
import matplotlib.pyplot as plt
import h5py


plot_freq_range = True
plot_max_eff = True
plot_R =        False
plot_L =        True
plot_power =    True
plot_modes =    True
plot_valids =   True


sim = "/force_2_modes_freq-grat_width"
with h5py.File("grating_coupler_results.h5", "r") as f: 
    x_vals = f[sim]["x_values"][:]
    y_vals = f[sim]["y_values"][:]
    X, Y = np.meshgrid(x_vals, y_vals)


def matrix_with_mask(matrix, matrix_mask, tolerance=0):
    mask = (matrix_mask <= tolerance)
    return np.where(mask, matrix, 0)


def find_max(matrix, n_th_max = 1):
    flat = matrix.ravel()

    idx = np.argpartition(flat, -n_th_max)[-n_th_max]
    max_pos = np.unravel_index(idx, matrix.shape)

    wgw =  y_vals[max_pos[0]]
    freq = x_vals[max_pos[1]]
    eff = matrix[max_pos]

    return max_pos, freq, wgw, eff


def do_plot(matrix, title, xlabel="Frequency", ylabel="Waveguide Width"):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, matrix, cmap='viridis')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel('Efficiency')
    t = "Grating Coupler N_WG=3.47 "
    if sim != "/force_2_modes_freq-grat_width":
        t += "at freq = 0.875 and grat_height = 0.252 "

    plt.title(f"{t}{title}")


def plot_region(reg, maxs_to_find=0, tolerance=0, get_second_mode=True, xlabel="Frequency", ylabel="Waveguide Width"):
    with h5py.File("grating_coupler_results.h5", "r") as f: 
        power = f[sim][f"total_power_{reg}"][:]
        mode0 = f[sim][f"mode0{reg}"][:]
        mode1 = f[sim][f"mode1{reg}"][:]
        if get_second_mode: mode2 = f[sim][f"mode2{reg}"][:]

    if plot_power:
        do_plot(power, "Total Power RIGHT" if reg == 'R' else "Total Power LEFT", xlabel, ylabel)

    if plot_modes:
        do_plot(mode0, "Fundamental Mode Power RIGHT" if reg == 'R' else "Fundamental Mode Power LEFT", xlabel, ylabel)
        do_plot(mode1, "First Mode Power RIGHT" if reg == 'R' else "First Mode Power LEFT", xlabel, ylabel)
        if get_second_mode: do_plot(mode2, "Second Mode Power RIGHT" if reg == 'R' else "Second Mode Power LEFT", xlabel, ylabel)

    if plot_valids:
        if get_second_mode:
            valid_power = matrix_with_mask(power, mode2, tolerance=tolerance)
            valid_modes = matrix_with_mask(np.add(mode0, mode1), mode2, tolerance=tolerance)

            do_plot(valid_power, "Valid Power RIGHT" if reg == 'R' else "Valid Power LEFT", xlabel, ylabel)
            do_plot(valid_modes, "Valid Modes 0 + 1 RIGHT" if reg == 'R' else "Valid Modes 0 + 1 LEFT", xlabel, ylabel)
        else:
            valid_modes = np.add(mode0, mode1)
            do_plot(valid_modes, "Valid Modes 0 + 1 RIGHT" if reg == 'R' else "Valid Modes 0 + 1 LEFT", xlabel, ylabel)


        for max in range(1, maxs_to_find+1):
            print(f"{max}-th max value found in valid modes:")
            max_pos, freq, wgw, eff = find_max(valid_modes, n_th_max=max)
            if get_second_mode:
                print(f"Position: {max_pos},   freq: {freq},   wg_height: {wgw},   eff: {eff}")
            else:
                print(f"Position: {max_pos},   Duty Cycle: {freq},   Depth Factor: {wgw},   eff: {eff}")

    plt.show()



def main():
    if plot_freq_range:
        if plot_R:
            plot_region('R', maxs_to_find=2)

        if plot_L:
            plot_region('L', maxs_to_find=2)

    if plot_max_eff:
        global sim, X, Y, x_vals, y_vals
        sim = "/depth-duty_cycle_MAX_freq0-875_height0-25"
        with h5py.File("grating_coupler_results.h5", "r") as f: 
            x_vals = f[sim]["x_values"][:]
            y_vals = f[sim]["y_values"][:]
            X, Y = np.meshgrid(x_vals, y_vals)

        if plot_R:
            plot_region('R', maxs_to_find=2, xlabel="Duty Cycle", ylabel="Depth Factor", get_second_mode=False)

        if plot_L:
            plot_region('L', maxs_to_find=2, xlabel="Duty Cycle", ylabel="Depth Factor", get_second_mode=False)


if __name__ == "__main__":
    main()