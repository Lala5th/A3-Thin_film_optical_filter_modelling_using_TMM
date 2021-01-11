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

a = MaterialTable.fromValue(1)
b = MaterialTable.fromValue(2)

wl = _np.arange(0.4,1,0.001)


depths = _np.linspace(0,6,1000)
Is = TransferMatrix([a] + [a]*2 + [b]*30 + [b],[1]*2 + [0.1]*30).getIntensityProfile(wl,0,False,zs=depths,function=False)
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
