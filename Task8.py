#!/usr/bin/ipython3
import material_table as _mtbl
import numpy as _np
import matplotlib.pyplot as plt
from TMM import TransferMatrix
from matplotlib.widgets import Slider


air = _mtbl.MaterialTable.fromMaterial('Vacuum')
BK7 = _mtbl.MaterialTable.fromMaterial('BK7')

# Task 8
plt.ion()


wls = _np.arange(0.4,1,0.005)

fig,ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
plt.set_cmap('gnuplot')
l = plt.plot(wls,TransferMatrix([air,_mtbl.MaterialTable.fromValue(1.1),BK7],0.15).solveTransmission(wls,0,False)[:,0])
plt.ylim([0,0.2])
plt.grid()
axn = plt.axes([0.25, 0.1, 0.65, 0.03])
axd = plt.axes([0.25, 0.15, 0.65, 0.03])
sn = Slider(axn, 'n',1,2,valinit=1.1,valstep=0.01)
sd = Slider(axd, 'd',0.01,0.5,valinit=0.15,valstep=0.01)

def update(val):
	n = sn.val
	d = sd.val
	l[0].set_ydata(TransferMatrix([air,_mtbl.MaterialTable.fromValue(n),BK7],d).solveTransmission(wls,0,False)[:,0])
	fig.canvas.draw_idle()

sn.on_changed(update)
sd.on_changed(update)
plt.show()
