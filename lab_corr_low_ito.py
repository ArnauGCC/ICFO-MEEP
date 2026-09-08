import numpy as np
from packages.utils import *
import matplotlib.pyplot as plt
import h5py


PR = [0.020558836802575912, 0.020927094244921028, 0.020718853914459906, 0.0003076256199174047, 0.01938297435886498, 0.00032376905160389977, 0.01965828740481319, 0.01985206540886831, 0.01935342267485922, 0.0003205825277677899, 0.01922808855406256, 0.01888732284732248, 0.01929706728495259, 0.01905581940628704]
S = [387.5381856483799, 387.97414038588215, 388.29936212075853, 388.5118383985217, 388.5756571994773, 388.56044653144, 388.6346310519957, 388.7579957037281, 388.73711160332164, 388.75872173620206, 388.77876114774085, 388.7478210619961, 388.74690126504015, 388.8398316830705]
PL = [0.001856850619431004, 0.00231374219755479, 0.0023033698360108087, 0.0005176421176283302, 0.002037007323057423, 0.00050233848332246, 0.0019764103526252524, 0.002100650312515662, 0.0020713288055138335, 0.000521605024316329, 0.002061233394163144, 0.0019173587476186923, 0.002038600116321716, 0.0019868757783339714]
res_factors = compute_freqs(20, 85, f_res=5)

def relative_errors(values):
    return [
        abs(values[i] - values[i - 1]) / abs(values[i - 1])
        for i in range(len(values))
    ]


plt.figure()
plt.plot(res_factors, (PR))
plt.xlabel('res_factor')
plt.ylabel(r'value')
plt.title("POWER Right")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(res_factors, PL)
plt.xlabel('res_factor')
plt.ylabel(r'value')
plt.title("POWER Left")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(res_factors, S)
plt.xlabel('res_factor')
plt.ylabel(r'value')
plt.title("Source")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(res_factors, relative_errors(PR))
plt.xlabel('res_factor')
plt.ylabel(r'value')
plt.title("RE POWER Right")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(res_factors, relative_errors(S))
plt.xlabel('res_factor')
plt.ylabel(r'value')
plt.title("RE Source")
plt.grid(True, which="both", alpha=0.3)


plt.show()
