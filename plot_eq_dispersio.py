import matplotlib.pyplot as plt
import numpy as np

f_min = 0.7
f_max = 2.5
f_res = 0.1
f_steps = round((f_max - f_min)/f_res) +1

n = 1.50                                # index waveguide

freqs = np.linspace(f_min, f_max, f_steps)
wg_factors = [np.float64(0.7857142857142857), np.float64(0.7576530612244897), np.float64(0.7755102040816326), np.float64(0.8137755102040817), np.float64(0.8137755102040817), np.float64(0.8714285714285714), np.float64(0.8622448979591837), np.float64(0.9005102040816326), np.float64(0.90625), np.float64(0.9183673469387755), np.float64(0.9183673469387755), np.float64(0.9566326530612245), np.float64(0.9566326530612245), np.float64(0.9566326530612245), np.float64(0.9948979591836735), np.float64(0.9948979591836735), np.float64(0.9964285714285714), np.float64(0.9964285714285714), np.float64(0.9955357142857143)]


wg_widths = wg_factors/(2*freqs*n)
eff_wg_width = [0.6882152669262575, 0.6974634137839019, 0.7014771913905148, 0.7049594411871288, 0.7075875802012354, 0.709520078659464, 0.7107515556515785, 0.7116505204116244, 0.7064925173008624, 0.7057978774299665, 0.7017332109972281, 0.694126457986922, 0.6866706801799823, 0.679922625277651, 0.6762803267971592, 0.6723638441510565, 0.6667557245319612, 0.6611362236561644, 0.6544778701219338]


plt.figure()
plt.plot(wg_widths, freqs)
plt.xlabel(r'Waveguide width')
plt.ylabel(r'Frequency')
plt.title("Equació de dispersió")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(wg_widths, eff_wg_width)
plt.xlabel(r'Waveguide width')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(freqs, eff_wg_width)
plt.xlabel(r'Frequency')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.show()


eff_pad_coupled = [0.7065039342466871, 0.7020189836842021, 0.7044647122850498, 0.7067162459103635, 0.7064381688120173, 0.7104333103306286, 0.7134411398202064, 0.712716942654223, 0.7058705669459218, 0.7055911792050616, 0.7029797957621692, 0.6937224114558426, 0.6869139882276475, 0.6821187779342313, 0.6761792325697554, 0.672712535803712, 0.667015122002146, 0.6613151594163749, 0.6604138605501378]
pad_coupled = [np.float64(0.18444444444444447), np.float64(0.19444444444444448), np.float64(0.17272222222222222), np.float64(0.16333333333333336), np.float64(0.14222222222222225), np.float64(0.1366666666666667), np.float64(0.1202716049382716), np.float64(0.11672839506172841), np.float64(0.10487654320987655), np.float64(0.10487654320987655), np.float64(0.10165432098765431), np.float64(0.0956172839506173), np.float64(0.09116049382716052), np.float64(0.08983055555555555), np.float64(0.07911666666666667), np.float64(0.07661728395061729), np.float64(0.07265432098765431), np.float64(0.07265432098765431), np.float64(0.07265432098765431)]

freqs = np.delete(freqs, 0)
eff_pad_coupled.pop(0)
pad_coupled.pop(0)

plt.figure()
plt.plot(pad_coupled, freqs)
plt.xlabel(r'Pad prism - waveguide')
plt.ylabel(r'Frequency')
plt.title("Max eff of coupled power")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(pad_coupled, eff_pad_coupled)
plt.xlabel(r'Pad prism - waveguide')
plt.ylabel(r'Efficiency')
plt.title("Coupled Power Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(freqs, eff_pad_coupled)
plt.xlabel(r'Frequency')
plt.ylabel(r'Efficiency')
plt.title("Coupled Power Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.show()


f_max = 1.9
f_steps = round((f_max - f_min)/f_res) +1
freqs = np.linspace(f_min, f_max, f_steps)


eff_pad_wg = [np.float64(0.11574369600569094), np.float64(0.1156749951459681), np.float64(0.11715390996197543), np.float64(0.1251989841841299), np.float64(0.1331709759029413), np.float64(0.1423440685585467), np.float64(0.146982330523643), np.float64(0.14687609031828383), np.float64(0.14841975956498726), np.float64(0.14404794620583763), np.float64(0.14430616219624184), np.float64(0.15703901441599122), np.float64(0.1564529104872858)]
pad_wg = [np.float64(0.3625277777777778), np.float64(0.36524999999999996), np.float64(0.33394444444444443), np.float64(0.28494444444444444), np.float64(0.25363888888888886), np.float64(0.2563611111111111), np.float64(0.22777777777777775), np.float64(0.20781481481481484), np.float64(0.19919444444444445), np.float64(0.20191667), np.float64(0.1932963), np.float64(0.17030864), np.float64(0.15881481)]

#eff_pad_wg.pop(0)
#pad_wg.pop(0)

plt.figure()
plt.plot(pad_wg, freqs)
plt.xlabel(r'Pad prism - waveguide')
plt.ylabel(r'Frequency')
plt.title("Max eff of waveguide power")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(pad_wg, eff_pad_wg)
plt.xlabel(r'Pad prism - waveguide')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Power Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.figure()
plt.plot(freqs, eff_pad_wg)
plt.xlabel(r'Frequency')
plt.ylabel(r'Efficiency')
plt.title("Waveguide Power Efficiency")
plt.grid(True, which="both", alpha=0.3)

plt.show()


"""
max_effs = []
max_pads = []
for i in range(0, len(freqs)):
    if eff_pad_wg[i] > eff_pad_coupled[i]:
        max_effs.append(eff_pad_wg[i])
        max_pads.append(pad_wg[i])
    else:
        max_effs.append(eff_pad_coupled[i])
        max_pads.append(pad_coupled[i])

print(str(max_effs))
print("\n \n")
print(str(max_pads))"""