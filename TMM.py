#!/usr/bin/ipython3
import numpy as _np
import material_table as _mtbl
import scipy.constants as _const

class TransferMatrix:
	'''Class containing the code and parameters necessary to perform transfer
	matrix calculations on a given stack.'''

	def __init__(self,materials,thicknesses=[]):
		'''The constructor takes 2 array like properties. Materials and
		thicknesses. The first material in the materials array is assumed
		as the material of the world, while the last is taken as the material of
		the substrate. The thicknesses parameter can be omitted.'''

		# Code to sanitise input. In theory thickness can be a single value
		# For materials it is just a precaution so that no datatype mismatch
		materials = _np.array(materials,ndmin=1)
		thicknesses = _np.array(thicknesses,ndmin=1)

		# Check if materials contains exactly 2 more elements than thicnesses
		# as the world and substrate do not have thickness
		if len(materials) - 2 != len(thicknesses):
			raise ValueError('The material array must strictly be 2 larger than \
			the thickness array!')

		self._stack = _np.array([materials[1:-1],thicknesses]).transpose()
		self._world = materials[0]
		self._substrate = materials[-1]


	def getLayer(self,i):
		'''This function would return layer i in the stack. Layer -1 corresponds
		to the world while Layer len(self._stack) corresponds to the substrate'''

		if i == -1:
			return [self._world, None]
		elif i == len(self._stack):
			return [self._substrate, None]
		elif i < len(self._stack) and i >= 0:
			return self._stack[i]

		raise KeyError('The layer number should be between -1 and ' + str(len(self._stack)))

	def getKz(self,i,l,theta,k0=None):
		'''Method to get kz from k0, lambda and theta for a given layer.
		k0 is optional to supply, but due to performance it is recommended
		for large stacks'''
		if k0 == None:
			k0 = _mtbl.getWavenumber(l)
		kx = k0*_np.sin(theta)
		return _np.sqrt(k0**2*self.getLayer(i)[0].getN(l)**2-kx**2)

	def getPropagationMatrix(self,i,l,theta,k0=None):
		'''Method for getting the propagation matrix for a given stack layer,
		incidence(theta) and wavelength(l)'''
		kz = self.getKz(i,l,theta,k0)
		# Raise error if layer number i does not make sense
		if not (i >= 0 and i < len(self._stack)):
			raise KeyError('The propagation matrix can only be calculated for \
			layers of the stack (i.e. for index between 0 and ' + str(len(self._stack)-1)\
			 + ' inclusive)')
		return _np.array([[_np.exp(1j*kz*self.getLayer(i)[1]),0],
						  [0,_np.exp(-1j*kz*self.getLayer(i)[1])]])

	def getInterfacialMatrix(self,i,l,theta,p,k0=None):
		'''Method for getting the interfacial transfer matrix for a given
		interface. The interface is between the ith and the (i+1)th layer. p
		decides whether the polarisation is p (True) or s (False)'''

		kz1 = self.getKz(i,l,theta,k0)
		kz2 = self.getKz(i+1,l,theta,k0)

		def Xs(pos):

			C = _np.sqrt(kz2/kz1)
			if pos:
				return (C + 1/C)/2
			return (C - 1/C)/2

		def Xp(pos):

			n1 = self.getLayer(i)[0].getN(l)
			n2 = self.getLayer(i+1)[0].getN(l)
			C = _np.sqrt(n2**2*kz1/(kz2*n1**2))
			if pos:
				return (C + 1/C)/2
			return (C - 1/C)/2

		if p:
			return _np.array([[Xp(True),Xp(False)],
							  [Xp(False),Xp(True)]])

		return _np.array([[Xs(True),Xs(False)],
						  [Xs(False),Xs(True)]])

	def getTotalTransferMatrix(self,l,theta,p):
		'''The method to compile the total transfer matrix of the system'''

		k0 = _mtbl.getWavenumber(l)

		transferMatrix = _np.array([[1,0],[0,1]])
		transferMatrix = transferMatrix @ self.getInterfacialMatrix(-1,l,theta,p,k0)
		for i in range(len(self._stack)):
			transferMatrix = transferMatrix @ self.getPropagationMatrix(i,l,theta,k0)
			transferMatrix = transferMatrix @ self.getInterfacialMatrix(i,l,theta,p,k0)

		return transferMatrix

	def solveTransmission(self,lambdas,thetas,ps):
		'''Method for returning the transmission and reflection of the stack.
		returns a len(lambdas)*len(thetas)*len(ps)*2 shaped array. Extraneous
		nesting is discarded via numpy.squeeze (i.e. [[[5,3]]] -> [5,3]).'''

		lambdas = _np.array(lambdas,ndmin=1)
		thetas = _np.array(thetas,ndmin=1)
		ps = _np.array(ps,ndmin=1)

		ret = []

		for l in lambdas:
			thetaVals = []

			for theta in thetas:
				pVals = []

				for p in ps:
					transferMatrix = self.getTotalTransferMatrix(l,theta,p)

					r = - transferMatrix[1,0]/transferMatrix[1,1]
					t =  transferMatrix[0,0] + transferMatrix[0,1]*r

					R = _np.real(r)**2 + _np.imag(r)**2
					T = _np.real(t)**2 + _np.imag(t)**2
					pVals.append((R,T))
				thetaVals.append(pVals)
			ret.append(thetaVals)
		ret = _np.array(ret)
		return _np.squeeze(ret)


if __name__ == '__main__':

	fails = 0
	tests = 0
	def test(x,e):
		global fails,tests
		'''Function to wrap the unit test results. 1e-10 for rounding errors.
		x is the measured value while e is the expected.'''
		tests += 1
		if (_np.abs(x-e)  >= 1e-10).any():
			print('\33[31mFAILED\33[0m')
			fails += 1
		else:
			print('\33[32mPASSED\33[0m')

	vacuum = _mtbl.MaterialTable.fromMaterial('Vacuum')
	propagationMatrixTest = TransferMatrix([vacuum,vacuum,vacuum],[1])
	print('Testing propagation matrix: ',end='')
	test(propagationMatrixTest.getPropagationMatrix(0,2*_const.pi,0),_np.array([[_np.exp(1j),0],[0,_np.exp(-1j)]]))

	print()
	if fails == 0:
		print('All tests passed')
	else:
		print(fails,'tests failed out of',tests)


	import matplotlib.pyplot as plt
	plt.ion()

	# Task 10

	air = _mtbl.MaterialTable.fromMaterial('air')
	gold = _mtbl.MaterialTable.fromMaterial('Au')
	BK7 = _mtbl.MaterialTable.fromMaterial('BK7')

	transferredI = []
	wavelengths = _np.arange(0.4,1,0.001)
	for i in range(1,11):
		TM = TransferMatrix([air,gold,BK7],i*10*1e-3)
		transferredI.append([i*10,TM.solveTransmission(wavelengths,0,True)])

	transferredI = _np.array(transferredI)

	for i in range(len(transferredI)):
		plt.plot(wavelengths,transferredI[i,1][:,1],label=str(transferredI[i,0])+'nm')
	plt.legend()
	plt.show()
