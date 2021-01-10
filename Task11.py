#!/usr/bin/ipython3
import numpy as _np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from TMM import TransferMatrix
from material_table import MaterialTable
import scipy.constants as _const
plt.ion()

air = MaterialTable.fromMaterial('air',0.4,1)
BK7 = MaterialTable.fromMaterial('BK7',0.4,1)
MgF2 = MaterialTable.fromMaterial('MgF2',0.4,1)
Ta2O5 = MaterialTable.fromMaterial('Ta2O5',0.4,1)

stack = [Ta2O5,MgF2]*2
stack.append(BK7)
stack.insert(0,air)

ds = _np.arange(0.01,0.2,0.001)

wl = _np.arange(0.4,1,0.001)

phases = _np.array(_np.mod([[_np.angle(TransferMatrix(stack,[d]*4).trackSingleBeam([-1,0,1,2,3],w,0,False)) for w in wl] for d in ds],2*_const.pi))

fig,ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
l = plt.plot(wl,phases[5])
plt.legend([-1,0,1,2,3])
axd = plt.axes([0.25, 0.15, 0.65, 0.03])
sd = Slider(axd, 'd',ds[0],ds[-1],valinit=ds[5],valstep=0.001)

def update(val):
	did = _np.where(_np.abs(ds - sd.val)<=1e-6)[0][0]
	[i.set_ydata(phases[did,:,l.index(i)]) for i in l]
	fig.canvas.draw_idle()

sd.on_changed(update)

plt.show()
