from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np

file_name = "/home/arnaugc/ICFO/MEEP/prism_results/DIFF_TMP-off2.25_fcen1.00_w0.28_al45_n1.5_pad0.375_df0.05.npz"
data = np.load(file_name)

wg_flux = data["wg_flux"]
src_flux = data["src_flux"]
refl_flux = data["refl_flux"]
freqs = data["freqs"]
df = data["df"]

plt.figure()
plt.plot(freqs,wg_flux)
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Flux")
plt.title("Waveguide SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

plt.figure()
plt.plot(freqs,refl_flux)
plt.xlabel(r'frequency $f (kHz)$')
#s = len(freqs)
#plt.xlim(freqs[int(4*s/10)], freqs[int(6*s/10)])
plt.ylabel("Flux")
plt.title("Reflected SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

plt.figure()
plt.plot(freqs,src_flux)
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Flux")
plt.title("Source SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))


plt.figure()
plt.plot(freqs,np.divide(refl_flux, src_flux))
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Flux")
plt.title("Reflected/Source SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

plt.figure()
plt.plot(freqs,np.divide(wg_flux, src_flux))
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Flux")
plt.title("Waveguide/Source SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))

plt.figure()
plt.plot(freqs,np.divide(wg_flux, src_flux) + np.divide(refl_flux, src_flux))
plt.xlabel(r'frequency $f (kHz)$')
plt.ylabel("Flux")
plt.title("Reflected/Source + Waveguide/Source SPD")
plt.grid(True, which="both", alpha=0.3)
plt.gca().xaxis.set_major_locator(MultipleLocator(df*0.1))


plt.show()