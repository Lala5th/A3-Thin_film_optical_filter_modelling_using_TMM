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

wl = _np.arange(0.4,1,0.001)

def getThicknesses(l=0.5,periods = 2):
	widths = []
	stack = [air] + [Ta2O5,MgF2]*periods + [BK7]
	TM = TransferMatrix(stack,[1,1]*periods)
	IM = [TM.getInterfacialMatrix(i,l,0,True) for i in range(-1,periods*2)]
	p0 = _np.mod(_np.angle(-IM[0][1,0]/IM[0][1,1]),2*_const.pi) # phase from immediate reflection
	phase = _np.mod(_np.angle(1-(IM[0][1,0]/IM[0][1,1])**2),2*_const.pi) # Phase change of thransmission through the first interface and back
	for i in range(periods*2):
		phase_diff = _np.mod(p0 - _np.mod(phase + _np.angle(-IM[i+1][1,0]/IM[i+1][1,1]),2*_const.pi),2*_const.pi)
		if (phase_diff < 1e-6):
			phase_diff = 2*_const.pi
		width = l*phase_diff/(4*_const.pi*TM.getLayer(i)[0].getRealN(l)) # 4pi/l *n*d = phase_diff
		widths.append(width)
		phase = _np.mod(phase+phase_diff+_np.angle(1-(IM[i+1][1,0]/IM[i+1][1,1])**2),2*_const.pi) # Phase change of thransmission through the ith interface and back

	return (stack,widths)



fig2 = plt.figure()
i = 1
while True:
	stack,d = getThicknesses(0.633,i)
	RT = TransferMatrix(stack,d).solveTransmission(0.633,0,True)
	Rprime = RT[0]/(RT[0]+RT[1])
	if Rprime > 0.9999:
		print(Rprime,'at',i)
		break
	i += 1
transmissions = TransferMatrix(stack,d).solveTransmission(wl,0,True)
reflectance = transmissions[:,0]/(transmissions[:,0] + transmissions[:,1])
plt.plot(wl,reflectance)
plt.ylabel('Reflectance')
plt.xlabel("$\lambda$ [$\mu$m]")
plt.ylim([0,1])
plt.show()

# Task 15

Au = MaterialTable.fromMaterial('Au')

phase = _np.angle(TransferMatrix(stack,d).solveTransmission(0.633,0,True,True)[0])
phase_Au = _np.angle(TransferMatrix([air,Au]).solveTransmission(0.633,0,True,True)[0])

print(phase,':',phase_Au)
