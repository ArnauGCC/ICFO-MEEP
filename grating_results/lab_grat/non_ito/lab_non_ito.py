import numpy as np
from packages.utils import *
import matplotlib.pyplot as plt
import h5py


def do_plot(X, Y, matrix, title, xlabel='Duty Cycle', ylabel='Depth Factor'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, matrix, cmap='viridis')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel('Efficiency')
    plt.title(title)


def plot_depth_duty():
    sim = "/lab_non_ito_chng_depth-duty"
    with h5py.File("grating_coupler_results.h5", "r") as f: 
        X, Y = np.meshgrid(f[sim]['x_values'][:], f[sim]['y_values'][:])

        t = "lab grating without ITO with n_cells=20 grat_height (PE layer) = 0.7 "
        do_plot(X, Y, f[sim]["total_power_R"], t+"TOTAL POWER Right")
        do_plot(X, Y, f[sim]["mode0R"], t+"MODE 0R")
        do_plot(X, Y, f[sim]["mode1R"], t+"MODE 1R")
        do_plot(X, Y, f[sim]["mode2R"], t+"MODE 2R")
        do_plot(X, Y, f[sim]["total_power_L"], t+"TOTAL POWER Left")

    plt.show()


def find_max(matrix, x_vals, y_vals, n_th_max = 1):
    flat = matrix.ravel()

    idx = np.argpartition(flat, -n_th_max)[-n_th_max]
    max_pos = np.unravel_index(idx, matrix.shape)

    x_val = x_vals[max_pos[1]]
    y_val =  y_vals[max_pos[0]]
    eff = matrix[max_pos]

    return max_pos, x_val, y_val, eff


def plot_height_depth():
    sim = "/lab_non_ito_chng_height-depth"
    with h5py.File("grating_coupler_results.h5", "r") as f: 
        x_vals = f[sim]['x_values'][:]
        y_vals = f[sim]['y_values'][:]
        X, Y = np.meshgrid(x_vals, y_vals)

        PR = f[sim]["total_power_R"]
        _0R = f[sim]["mode0R"]
        _1R = f[sim]["mode1R"]
        _2R = f[sim]["mode2R"]
        _3R = f[sim]["mode3R"]
        SUMA = np.add(np.add(np.add(_0R, _1R), _2R), _3R)

        t = "lab grating without ITO with n_cells=20 "

        X, Y = np.meshgrid(x_vals, y_vals)
        do_plot(X, Y, PR, t+"Power", ylabel='Height', xlabel='Depth Factor')
        do_plot(X, Y, SUMA, t+"SUMA", ylabel='Height', xlabel='Depth Factor')
        do_plot(X, Y, _0R, t+"MODE 0R", ylabel='Height', xlabel='Depth Factor')
        do_plot(X, Y, _1R, t+"MODE 1R", ylabel='Height', xlabel='Depth Factor')
        do_plot(X, Y, _2R, t+"MODE 2R", ylabel='Height', xlabel='Depth Factor')
        do_plot(X, Y, _3R, t+"MODE 3R", ylabel='Height', xlabel='Depth Factor')


        max_pos, depth, height, eff = find_max(SUMA, x_vals, y_vals)
        print(f"Position: {max_pos},   depth: {depth},   height: {height},   eff: {eff}")


    plt.show()


def plot_period_height():
    sim = "/lab_non_ito_chng_period-height"
    with h5py.File("grating_coupler_results.h5", "r") as f: 
        x_vals = f[sim]['x_values'][:]
        y_vals = f[sim]['y_values'][:]
        X, Y = np.meshgrid(x_vals, y_vals)

        PR = f[sim]["total_power_R"]
        _0R = f[sim]["mode0R"]
        _1R = f[sim]["mode1R"]
        _2R = f[sim]["mode2R"]
        _3R = f[sim]["mode3R"]
        SUMA = np.add(np.add(np.add(_0R, _1R), _2R), _3R)


        t = "lab grating without ITO with n_cells=20 "

        do_plot(X, Y, PR, t+"Power R", ylabel='Period', xlabel='Height')
        do_plot(X, Y,  SUMA, t+"SUMA", ylabel='Period', xlabel='Height')
        do_plot(X, Y, _0R, t+"MODE 0R", ylabel='Period', xlabel='Height')
        do_plot(X, Y, _1R, t+"MODE 1R", ylabel='Period', xlabel='Height')
        do_plot(X, Y, _2R, t+"MODE 2R", ylabel='Period', xlabel='Height')
        do_plot(X, Y, _3R, "MODE 3R", ylabel='Period', xlabel='Height')

        max_pos, x_val, y_val, eff = find_max(SUMA, x_vals, y_vals)
        print(f"Position: {max_pos},   height: {x_val},   period: {y_val},   eff: {eff}")


        plt.show()


if __name__ == "__main__":
    plot_depth_duty()
    plot_height_depth()
    plot_period_height()
