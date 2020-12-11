#!/usr/bin/ipython3
import material_table as _mtbl
import numpy as _np
import matplotlib.pyplot as plt
import scipy.constants as _const
from TMM import TransferMatrix
from matplotlib.widgets import Slider


air = _mtbl.MaterialTable.fromMaterial('air')
BK7 = _mtbl.MaterialTable.fromMaterial('BK7')

# Task 8

ns = _np.arange(1,2,0.01)
thicknesses = _np.arange(0.01,0.5,0.01)
wls = _np.arange(0.4,1,0.005)
thetas = _np.arange(0,_const.pi,0.1)

customMaterials = [_mtbl.MaterialTable.fromValue(n) for n in ns]
transmission = []
for i in range(ns.size):
	temp_d_data = []
	for d in thicknesses:
		temp_theta_data = []
		for t in thetas:
			temp_theta_data.append([d,TransferMatrix([air,customMaterials[i],BK7],d).solveTransmission(wls,t,False)[:,0]])
		temp_d_data.append(temp_theta_data)
	transmission.append(temp_d_data)
transmission = _np.array(transmission)

fig,ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
plt.set_cmap('gnuplot')
l = plt.plot(wls,transmission[5,2,7,1])
plt.ylim([0,0.2])
plt.ylabel('Reflectance')
plt.xlabel('$\lambda$ [$\mu$m]')
plt.grid()
axn = plt.axes([0.25, 0, 0.65, 0.03])
axd = plt.axes([0.25, 0.05, 0.65, 0.03])
axt = plt.axes([0.25, 0.10, 0.65, 0.03])
sn = Slider(axn, 'n',ns[0],ns[-1],valinit=ns[5],valstep=0.01)
sd = Slider(axd, 'd',thicknesses[0],thicknesses[-1],valinit=thicknesses[2],valstep=0.01)
st = Slider(axt, 'theta',thetas[0],thetas[-1],valinit=thetas[7],valstep=0.1)

def update(val):
	nid = _np.where(_np.abs(ns - sn.val)<=1e-6)[0][0]
	did = _np.where(_np.abs(thicknesses - sd.val)<=1e-6)[0][0]
	tid = _np.where(_np.abs(thetas - st.val)<=1e-6)[0][0]
	l[0].set_ydata(transmission[nid,did,tid,1])
	fig.canvas.draw_idle()

sn.on_changed(update)
sd.on_changed(update)
st.on_changed(update)
plt.show()
