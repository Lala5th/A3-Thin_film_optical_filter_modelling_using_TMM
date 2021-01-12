#!/usr/bin/ipython3
import numpy as _np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from TMM import TransferMatrix
from material_table import MaterialTable
import scipy.constants as _const
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
plt.ion()


AlAs = MaterialTable.fromCsv('AlAs.txt',delimiter='\t',skiprows=1,wavelength_multiplier=1e-3)
SiO2 = MaterialTable.fromCsv('SiO2.txt',delimiter='\t',skiprows=1,wavelength_multiplier=1e-3)
air = MaterialTable.fromMaterial('air')

TM = TransferMatrix([air,AlAs,SiO2,air],[0.01,0.25])

ref = _np.loadtxt('air-AlAs_10-SiO2_250-air.csv',skiprows=1,delimiter='\t')
ref[:,0]*=1e-3
wl = ref[:,0]

r,t = TM.solveTransmission(wl,0,True).transpose()
fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
ax2.plot(wl,t-ref[:,1])
ax2.set_ylabel('Difference [a.u.]')
ax2.set_ylim([-5e-4,5e-4])
ax1.plot(wl,ref[:,1])
ax1.plot(wl,t)
ax1.set_ylabel('Transmittance [a.u.]')
ax2.set_xlabel("$\lambda$ [$\mu$m]")
ax1.set_ylim([0.5,1])
ax2.grid()
ax1.grid()
plt.show()
