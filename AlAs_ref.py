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

wl = _np.arange(0.4,1,0.001)
ref = _np.loadtxt('air-AlAs_10-SiO2_250-air.csv',skiprows=1,delimiter='\t')
ref[:,0]*=1e-3
reference = interp1d(ref[:,0],ref[:,1])

r,t = TM.solveTransmission(wl,0,True).transpose()
plt.plot(wl,t-reference(wl))
plt.ylabel('Difference')
plt.xlabel("$\lambda$ [$\mu$m]")
plt.ylim([-1e-2,1e-2])
plt.show()
