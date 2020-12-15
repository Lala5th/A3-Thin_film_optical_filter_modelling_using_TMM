#!/usr/bin/ipython3
import numpy as _np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from TMM import TransferMatrix
from material_table import MaterialTable
import scipy.constants as _const
from scipy.optimize import curve_fit
plt.ion()

air = MaterialTable.fromMaterial('air',0.4,1)
BK7 = MaterialTable.fromMaterial('BK7',0.4,1)
MgF2 = MaterialTable.fromMaterial('MgF2',0.4,1)
Ta2O5 = MaterialTable.fromMaterial('Ta2O5',0.4,1)

initStack = [Ta2O5,MgF2]


wl = _np.arange(0.5,0.7,0.005)

def fit(layers = 2):
	widths = []
	for i in range(layers):
		stack = initStack*(i+1)
		stack.append(BK7)
		stack.insert(0,air)
		fitfunc = lambda l,d1,d2 : TransferMatrix(stack,[_np.abs(d1)]+[_np.abs(d2)]+widths).solveTransmission(l,0,True)[:,0]
		fit,cov = curve_fit(fitfunc,wl,[1]*len(wl),p0=[0.08,0.08],sigma=1-fitfunc(wl,0.08,0.08))
		widths.append(_np.abs(fit[0]))
		widths.append(_np.abs(fit[1]))

	return (stack,widths)



fig2 = plt.figure()
stack,d = fit()
plt.plot(wl,TransferMatrix(stack,d).solveTransmission(wl,0,True)[:,0])
plt.ylim([0,1])
plt.show()
