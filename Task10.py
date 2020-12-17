#!/usr/bin/ipython3
import material_table as _mtbl
import numpy as _np
import matplotlib.pyplot as plt
from TMM import TransferMatrix

air = _mtbl.MaterialTable.fromMaterial('air')
gold = _mtbl.MaterialTable.fromMaterial('Au')
BK7 = _mtbl.MaterialTable.fromMaterial('BK7')


# Task 10

transferredI = []
wavelengths = _np.arange(0.4,1,0.001)
for i in range(1,21):
	TM = TransferMatrix([air,gold,BK7],i*10*1e-3)
	transferredI.append([i*10,TM.solveTransmission(wavelengths,0,True)])

transferredI = _np.array(transferredI)

plt.figure()
for i in range(len(transferredI)):
	plt.plot(wavelengths,transferredI[i,1][:,1],label=str(transferredI[i,0])+'nm')
plt.legend()
plt.ylabel('Transmittance')
plt.xlabel("$\lambda$ [$\mu$m]")
plt.yscale('log')
plt.show()
plt.figure()
for i in range(len(transferredI)):
	plt.plot(wavelengths,transferredI[i,1][:,0],label=str(transferredI[i,0])+'nm')
plt.legend()
plt.ylabel('Reflectance')
plt.xlabel("$\lambda$ [$\mu$m]")
plt.show()
plt.figure()
for i in range(len(transferredI)):
	plt.plot(wavelengths,1 - transferredI[i,1][:,1] - transferredI[i,1][:,0],label=str(transferredI[i,0])+'nm')
plt.legend()
plt.ylabel('Absorption')
plt.xlabel("$\lambda$ [$\mu$m]")
plt.show()
