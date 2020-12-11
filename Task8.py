#!/usr/bin/ipython3
import material_table as _mtbl
import numpy as _np
import matplotlib.pyplot as plt
from TMM import TransferMatrix
from matplotlib.widgets import Slider


air = _mtbl.MaterialTable.fromMaterial('air')
BK7 = _mtbl.MaterialTable.fromMaterial('BK7')

# Task 8

ns = _np.arange(1,2,0.01)
thicknesses = _np.arange(0.01,0.5,0.01)
wls = _np.arange(0.4,1,0.005)

customMaterials = [_mtbl.MaterialTable.fromValue(n) for n in ns]
transmission = []
for i in range(ns.size):
	temp_d_data = []
	for d in thicknesses:
		temp_d_data.append([d,TransferMatrix([air,customMaterials[i],BK7],d).solveTransmission(wls,0,False)[:,0]])
	transmission.append(temp_d_data)
transmission = _np.array(transmission)

fig,ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
plt.set_cmap('gnuplot')
l = plt.plot(wls,transmission[5,2,1])
plt.ylim([0,0.2])
plt.grid()
axn = plt.axes([0.25, 0.1, 0.65, 0.03])
axd = plt.axes([0.25, 0.15, 0.65, 0.03])
sn = Slider(axn, 'n',ns[0],ns[-1],valinit=ns[5],valstep=0.01)
sd = Slider(axd, 'd',thicknesses[0],thicknesses[-1],valinit=thicknesses[2],valstep=0.01)

def update(val):
	nid = _np.where(_np.abs(ns - sn.val)<=1e-6)[0][0]
	did = _np.where(_np.abs(thicknesses - sd.val)<=1e-6)[0][0]
	l[0].set_ydata(transmission[nid,did,1])
	fig.canvas.draw_idle()

sn.on_changed(update)
sd.on_changed(update)
plt.show()
