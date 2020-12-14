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

stack = [Ta2O5,MgF2]

initStack= stack*2
initStack.append(BK7)
initStack.insert(0,air)

wl = _np.arange(0.4,1,0.005)

fig,ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
TM = TransferMatrix(initStack,[0.01,0.01]*2)
l = plt.plot(wl,_np.mod(_np.angle(TM.trackSingleBeam([-1,0,1,2,3],wl,0,False)),2*_const.pi))
plt.legend([-1,0,1,2,3])
plt.ylim([0,2*_const.pi])
axn = plt.axes([0.25, 0.12, 0.65, 0.03])
axd1 = plt.axes([0.25, 0.07, 0.65, 0.03])
axd2 = plt.axes([0.25, 0.02, 0.65, 0.03])
sd1 = Slider(axd1, 'd1',0.001,0.5,valinit=0.01,valstep=0.001)
sd2 = Slider(axd2, 'd2',0.001,0.5,valinit=0.01,valstep=0.001)
sn = Slider(axn, 'n',2,100,valinit=2,valstep=1)

plt.show()

fig2 = plt.figure()

reflectance, = plt.plot(wl,TM.solveTransmission(wl,0,True)[:,0])
plt.ylim([0,1])
plt.show()

def update(val):
	d1 = sd1.val
	d2 = sd2.val
	n = int(sn.val)
	s = stack*n
	s.append(BK7)
	s.insert(0,air)
	TM = TransferMatrix(s,[d1,d2]*n)
	yval = _np.mod(_np.angle(TM.trackSingleBeam([-1,0,1,2,3],wl,0,False)),2*_const.pi)
	[curve.set_ydata(yval[:,i]) for i,curve in enumerate(l)]
	reflectance.set_ydata(TM.solveTransmission(wl,0,True)[:,0])
	fig.canvas.draw_idle()
	fig2.canvas.draw_idle()

sd1.on_changed(update)
sd2.on_changed(update)
sn.on_changed(update)
