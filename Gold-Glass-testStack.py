#!/usr/bin/ipython3
import numpy as _np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from TMM import TransferMatrix
from material_table import MaterialTable
import scipy.constants as _const
from scipy.optimize import curve_fit
import matplotlib.colors as colors

plt.ion()

air = MaterialTable.fromMaterial('air',0.4,1)
BK7 = MaterialTable.fromMaterial('BK7',0.4,1)
Au = MaterialTable.fromMaterial('Au',0.4,1)

wl = _np.arange(0.4,1,0.001)


depths = _np.linspace(0,6,1000)
Is = TransferMatrix([air] + [BK7,Au]*5 + [BK7],[1,0.01]*5).getIntensityProfile(wl,0,False,xs=depths,function=False)
cmap = plt.get_cmap('gnuplot')
norm = colors.LogNorm()
plt.figure()
plt.pcolormesh(wl,depths,Is[:,:,0].transpose(),cmap=cmap,norm=norm)
plt.colorbar()
plt.show()
plt.figure()
plt.pcolormesh(wl,depths,Is[:,:,1].transpose(),cmap=cmap,norm=norm)
plt.colorbar()
plt.show()
