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

		if (thicknesses <= 0).any():
			raise ValueError('The layer thicknesses must be strictly positive')

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
		kx = k0*self.getLayer(i)[0].getN(l)*_np.sin(theta)
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
			Xpos = Xp(True)
			Xneg = Xp(False)
		else:
			Xpos = Xs(True)
			Xneg = Xs(False)

		return _np.array([[Xpos,Xneg],
						  [Xneg,Xpos]])

	def getTotalTransferMatrix(self,l,theta,p):
		'''The method to compile the total transfer matrix of the system'''

		k0 = _mtbl.getWavenumber(l)

		transferMatrix = _np.array([[1,0],[0,1]])
		transferMatrix = transferMatrix @ self.getInterfacialMatrix(-1,l,theta,p,k0)
		for i in range(len(self._stack)):
			transferMatrix =  self.getPropagationMatrix(i,l,theta,k0) @ transferMatrix
			transferMatrix = self.getInterfacialMatrix(i,l,theta,p,k0) @ transferMatrix

		return transferMatrix

	def solveTransmission(self,lambdas,thetas,ps,complex=False):
		'''Method for returning the transmission and reflection of the stack.
		returns a len(lambdas)*len(thetas)*len(ps)*2 shaped array. Extraneous
		nesting is discarded via numpy.squeeze (i.e. [[[5,3]]] -> [5,3]).
		The complex parameter is optional and if set to True returns the
		transmitted complex wave amplitude instead of the intensity.'''

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

					if complex:
						pVals.append((r,t))
					else:
						R = _np.real(r)**2 + _np.imag(r)**2
						T = _np.real(t)**2 + _np.imag(t)**2
						pVals.append((R,T))
				thetaVals.append(pVals)
			ret.append(thetaVals)
		ret = _np.array(ret)
		return _np.squeeze(ret)

	def getIntensityProfile(self,lambdas,thetas,ps,zs=None,function=True,complex=False):
		'''A function for returning the intensity at a given depth within the stack.
		It returns the wave intensity/amplitude for a given depth x within the
		stack. It can also be asked to return a function that can be used to get
		these values for a given x later on.'''

		interface_locations = [sum(self._stack[:i,1]) for i in range(len(self._stack)+1)]

		def	getInterface(depth):
			ret = 0
			for i,loc in enumerate(interface_locations):
				if depth <= loc:
					break
				ret = i
			return ret

		# Function to generate output function
		# Must be nested otherwise namespaces are messed up
		def getFunction(layerbeams,ks):
			def func(depth):
				# Get location within stack and layer
				i = getInterface(depth)
				k = ks[i]
				distance = depth - interface_locations[i]
				propagation = _np.array([[_np.exp(1j*k*distance),0],
										 [0,_np.exp(-1j*k*distance)]])
				beam = propagation @ layerbeams[i] # This is fine as propagation is diagonal

				if complex:
					return beam
				return _np.real(beam)**2 + _np.imag(beam)**2
			return func

		lambdas = _np.array(lambdas,ndmin=1)
		thetas = _np.array(thetas,ndmin=1)
		ps = _np.array(ps,ndmin=1)
		zs = _np.array(zs,ndmin=1)

		funcs = []

		for l in lambdas:
			thetaVals = []

			for theta in thetas:
				pVals = []

				for p in ps:

					propagationMatrix = [self.getPropagationMatrix(i,l,theta) for i in range(len(self._stack))]
					interfacialMatrix = [self.getInterfacialMatrix(i,l,theta,p) for i in range(-1,len(self._stack))]
					kzs = [self.getKz(i,l,theta) for i in range(len(self._stack)+1)]

					r,t = self.solveTransmission(l,theta,p,complex=True)
					a = [1,r] # The beam at incidence

					# Get the partial values for the beams at the interfaces
					# This speeds up calculations if in function mode, but can
					# be slow if only a single value is needed from the top of
					# the stack.
					a = interfacialMatrix[0] @ a
					interface_beams = [a]
					for i in range(len(self._stack)):
						a = a @ propagationMatrix[i]
						a = a @ interfacialMatrix[i+1]
						interface_beams.append(a)
					pVals.append(getFunction(interface_beams,kzs))
				thetaVals.append(pVals)
			funcs.append(thetaVals)

		if function:
			return _np.squeeze(_np.array(funcs))

		# Calculate the actual values
		ret = []
		for i in funcs:
			jval = []
			for j in i:
				kval = []
				for k in j:
					zval = []
					for z in zs:
						zval.append(k(z))
					kval.append(zval)
				jval.append(kval)
			ret.append(jval)
		return _np.squeeze(_np.array(ret))




	def trackSingleBeam(self,interfaces,lambdas,theta,p):
		'''Tracks a single beam through the stack. The interfaces array
		contains the interfaces on which it reflects.'''

		interfaces = _np.array(interfaces,ndmin=1)

		# Check for sane input
		for interface in interfaces:
			interface = _np.array(interface,ndmin=1)
			if len(interface) != 0:
				if interface.max() >= len(self._stack) or interface.min() < -1:
					raise ValueError('Invalid interface number encountered.')
			prev = -5 # Previous interface location
			forward = True # Whether the beam is forward propagating
			for i in interface:
				if forward and prev >= i:
					raise ValueError(f'Unable to perform calculation. Invalid route: {interface}')
				elif not forward and prev <= i:
					raise ValueError(f'Unable to perform calculation. Invalid route: {interface}')
				prev = i
				forward = not forward

		ret = []
		lambdas = _np.array(lambdas,ndmin=1)
		for l in lambdas:

			k0 =  _mtbl.getWavenumber(l)
			propagationMatrix = [self.getPropagationMatrix(i,l,theta,k0) for i in range(len(self._stack))]
			interfacialMatrix = [self.getInterfacialMatrix(i,l,theta,p,k0) for i in range(len(self._stack))]
			interfacialMatrix.insert(0,self.getInterfacialMatrix(-1,l,theta,p,k0))

			beams = []
			for interface in interfaces:
				loc = -1 # To track the location
				forward = True # To track propagation direction
				beam = 1 # To track the phase and amplitude of the beam
				interface = _np.append(interface,len(self._stack)+1) # So that the escape from the stack can be evaluated
				interface = _np.array(interface,ndmin=1)
				for i in interface:
					while True:
						if loc == i:
							break
						if (loc == -1 and (not forward)) or (loc == len(self._stack) and forward):
							break
						if forward:
							r = -beam*interfacialMatrix[loc+1][1,0]/interfacialMatrix[loc+1][1,1]
							beam = beam*interfacialMatrix[loc+1][0,0] + r*interfacialMatrix[loc+1][0,1]
							loc += 1
							if loc != len(self._stack):
								beam *= propagationMatrix[loc][0,0]
						else:
							beam = beam/interfacialMatrix[loc+1][1,1]
							loc -= 1
							if loc != -1:
								beam /= propagationMatrix[loc][1,1]
						if (loc == -1 and (not forward)) or (loc == len(self._stack) and forward):
							break
					if (loc == -1 and (not forward)) or (loc == len(self._stack) and forward):
						break
					# If it was forward propagating then after the reflection
					# it will propagate backward so proper solutions used
					if forward:																		# [a b][i] [t]
						beam = -beam*interfacialMatrix[loc+1][1,0]/interfacialMatrix[loc+1][1,1]	# [c d][r] [0]
						forward = False
						if(loc != -1):
							beam /= propagationMatrix[loc][1,1]										# [a 0]  [0] [0]
							loc -= 1																# [0 a/] [t] [i]
							continue
						else:
							break
					else:																			# [a b][0] [r]
						beam = interfacialMatrix[loc+1][0,1]*beam/interfacialMatrix[loc+1][1,1]		# [c d][t] [i]
						forward = True
						if(loc != len(self._stack)):
							beam *= propagationMatrix[loc][0,0]
							loc += 1
							continue
						else:
							break

				# Get final interface transmission
				if forward:
					r = -beam*interfacialMatrix[len(self._stack)-1][1,0]/interfacialMatrix[len(self._stack)-1][1,1]
					beam = beam*interfacialMatrix[len(self._stack)-1][0,0] + r*interfacialMatrix[len(self._stack)-1][0,1]
				else:
					beam = beam/interfacialMatrix[0][1,1]

				beams.append(beam)

			ret.append(beams)

		return _np.squeeze(_np.array(ret))

__all__ = ['TransferMatrix']


if __name__ == '__main__':

	fails = 0
	tests = 0
	def test(x,e):
		global fails,tests
		'''Function to wrap the unit test results. 1e-10 for rounding errors.
		x is the measured value while e is the expected.'''
		tests += 1
		if (_np.abs(x-e)  >= 1e-6).any():
			print('\33[31mFAILED\33[0m')
			fails += 1
		else:
			print('\33[32mPASSED\33[0m')

	vacuum = _mtbl.MaterialTable.fromMaterial('Vacuum')
	propagationMatrixTest = TransferMatrix([vacuum,vacuum,vacuum],[1])
	print('Testing propagation matrix: ',end='')
	test(propagationMatrixTest.getPropagationMatrix(0,2*_const.pi,0),_np.array([[_np.exp(1j),0],[0,_np.exp(-1j)]]))
	print('\nTesting single beam tracking')
	BK7 = _mtbl.MaterialTable.fromMaterial('BK7')
	print('Single reflection BK7 substrate:',end='')
	TM0 = TransferMatrix([vacuum,BK7],[])
	test(TM0.solveTransmission(0.5,0,True)[0],_np.abs(TM0.trackSingleBeam(-1,0.5,0,True))**2)
	print('Single reflection BK7 substrate, vacuum layer:',end='')
	TM1 = TransferMatrix([vacuum,vacuum,BK7],[1])
	test(TM1.solveTransmission(0.5,0,True)[0],_np.abs(TM1.trackSingleBeam(0,0.5,0,True))**2)
	print('Single reflection BK7 layer, vacuum layer interface:',end='')
	TM2 = TransferMatrix([vacuum,vacuum,BK7,BK7],[1,1])
	test(TM2.solveTransmission(0.5,0,True)[0],_np.abs(TM2.trackSingleBeam(0,0.5,0,True))**2)
	print('Single reflection BK7 substrate, 2 vacuum layers:',end='')
	TM2A = TransferMatrix([vacuum,vacuum,vacuum,BK7],[1,1])
	test(TM2A.solveTransmission(0.5,0,True)[0],_np.abs(TM2A.trackSingleBeam(1,0.5,0,True))**2)

	print()
	if fails == 0:
		print('All tests passed')
	else:
		print(fails,'tests failed out of',tests)
