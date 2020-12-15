#!/usr/bin/ipython3
import material_table as _mtbl
import numpy as _np
import matplotlib.pyplot as plt
from TMM import TransferMatrix
from matplotlib.widgets import Slider
from scipy.optimize import curve_fit


air = _mtbl.MaterialTable.fromMaterial('Vacuum')
BK7 = _mtbl.MaterialTable.fromMaterial('BK7')

# Task 8
plt.ion()

fitfunc = lambda l, n, d : TransferMatrix([air,_mtbl.MaterialTable.fromValue(n),BK7],_np.abs(d)).solveTransmission(wls,0,False)[:,0]

wls = _np.arange(0.4,1,0.005)

fit, cov = curve_fit(fitfunc,wls,0,p0=[1.24,0.14])

fig,ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
plt.set_cmap('gnuplot')
plt.plot(wls,fitfunc(wls,*fit))
plt.ylim([0,0.2])
plt.ylabel('Reflectance')
plt.xlabel("$\lambda$ [$\mu$m]")
plt.grid()
plt.show()
print(fit)
print(_np.sqrt(_np.diag(cov)))
