#!/usr/bin/ipython3
import yaml as _yml
import io as _io
import numpy as _np
import re as _re
from scipy.interpolate import interp1d as _interp1d
import scipy.constants as _const

class MaterialTable:
	'''This class initialises the database, while offering access to specific
	interpolations on the data. getN(wavelength) returns the n for a specific
	wavelength, while getEc(wavelength) returns the extinction coeffictent.
	In both cases wavelength is in units of um'''

	# The object that will contain the datafiles' location
	library = None

	@staticmethod
	def initDatabase():
		'''This method loads the map for the database'''

		if MaterialTable.library != None:
			print("Reloading database mapping file.")

		try:
			with open(r'./database/parsed_lib.yml') as file:

				# Loading yml mapping.
				MaterialTable.library = _yml.full_load(file)

		except FileNotFoundError:
			print('Library file at ./database/parsed_lib.yml was not found')

	@staticmethod
	def getK(l):
		'''A function for transforming wavelength in um to wavenumber
		in um^-1'''
		return 2*_const.pi/l

	def __init__(self,material,lmin=None,lmax=None,allow_interpolation=True,\
				 require_k=False):
		'''Constructor to create a material table for material. lmin and lmax
		are optional parameters specifying the limits within which the
		dataset is expected to be valid in units of um.
		allow_interpolation decides whether interpolation is acceptable for
		the refractive index data.
		require_k is a parameter specifying, whether data for k is also
		desired. Datasets with k specified as well are always preffered,
		and if found will be used regardless of require_k. k is the extinction
		coefficient in this context.'''

		if material == None or material == 'Vacuum':
			print('Loading custom material: Vacuum')
			# To remain consistent with instance variables
			self.datafile = None
			self.dataN = None
			self.dataK = None
			self.getN = lambda l : 1
			self.getEc = lambda l : 0
			return # call return so search would not be called

		self.findMaterial(material,lmin,lmax,allow_interpolation,require_k)

	def findMaterial(self,material,lmin=None,lmax=None,allow_interpolation=True,\
					 require_k=False):
		'''A simple function for parsing the library and to find the requested
		material. lmin and lmax are the bounds for which the dataset
		is expected values in units of um
		allow_interpolation decides whether interpolation is acceptable for
		the refractive index data.
		require_k is a parameter specifying, whether data for k is also
		desired. Datasets with k specified as well are always preffered,
		and if found will be used regardless of require_k. k is the extinction
		coefficient in this context.'''

		# Check so undefined behaviour is avoided
		if type(material) != str:
			raise TypeError("Input type for material must be string!")

		# The collection of data on the material
		mat = None

		for category_key in MaterialTable.library.keys():

			category = MaterialTable.library[category_key]['content']

			if material in category.keys():
				mat = category[material]
				break

		# Check if database was found
		if mat == None:
			raise KeyError('Could not find material ' + material + ' in database.')

		else:

			print('Found material database for  ' + material + ' in category: ' + \
				   category_key)
			print('Found ' + str(len(mat['content'].items())) +\
			 	  ' datasets within material database for ' + mat['name'])

			self.parseDataset(mat['content'],lmin,lmax,allow_interpolation,require_k)

	def parseDataset(self,database,lmin=None,lmax=None,allow_interpolation=True,\
					 require_k=False):
		'''This method is responsible for checking and loading the dataset
		once the material database was found in findMaterial, or otherwise.
		lmin and lmax are the optional bounds of wavelength we are interested in
		the material in units of um.
		allow_interpolation decides whether interpolation is acceptable for
		the refractive index data.
		require_k is a parameter specifying, whether data for k is also
		desired. Datasets with k specified as well are always preffered,
		and if found will be used regardless of require_k. k is the extinction
		coefficient in this context.'''

		if (lmin == None and lmax != None): # Sanitising the input
			lmin = lmax
		elif (lmin != None and lmax == None):
			lmax = lmin
		elif (lmin == None and lmax == None):
			lmin = _np.inf
			lmax = -_np.inf # This makes the check below always return True
		elif (lmin > lmax):
			raise ValueError('lmin must be strictly smaller than lmax')

		goodDatasets = []
		for dataset in database:
			with open('./database/data/' + database[dataset]['data']) as file:
				datafile = _yml.full_load(file)
			for data in datafile['DATA']:
				# Filtering out not n data
				if not (('tabulated n' in data['type'] and allow_interpolation) \
					    or 'formula' in data['type']):
					continue
				# If formula is in the data type, the range is specified
				if 'formula' in data['type']:
					range = _np.loadtxt(_io.StringIO(data['wavelength_range']))
				else:
					tabulated = _np.loadtxt(_io.StringIO(data['data']))
					range = _np.array([tabulated[0,0],tabulated[-1,0]])
				if(range[0] < lmin and range[1] > lmax):
					# Append data for the good datasets. The structure is:
					# [dataset, N-index,range-N, datafile, K-index,range-K]
					# as k is not tested yet, it is None as well as its range
					goodDatasets.append([dataset,datafile['DATA'].index(data),range,datafile,None,None])

		# Check for k data
		contains_k = []
		for d in goodDatasets:
			datafile = d[3]
			for data in datafile['DATA']:
				if not ('tabulated k' in data['type'] or 'tabulated nk' in data['type']):
					continue
				tabulated = _np.loadtxt(_io.StringIO(data['data']))
				range = _np.array([tabulated[0,0],tabulated[-1,0]])
				if(range[0] < lmin and range[1] > lmax):
					contains_k.append(d)
					contains_k[-1][4] = datafile['DATA'].index(data)
					contains_k[-1][5] = range

		if len(goodDatasets) == 0:
			raise KeyError('No datasets satisfying requirements found')
		print('Found',len(goodDatasets),'datasets, which satisfied requirements')
		if len(contains_k) != 0:
			goodDatasets = contains_k
		elif require_k:
			raise KeyError('No k data was found')
		print('Out of which',len(contains_k),'datasets contain k data')
		print('Using ' + goodDatasets[0][0])
		print('Name:',database[goodDatasets[0][0]]['name'])
		if 'division' in database[goodDatasets[0][0]].keys():
			print('Division:',database[goodDatasets[0][0]]['division'])
		print('Location: ./database/data/'+database[goodDatasets[0][0]]['data'])
		print('Dataset type:',goodDatasets[0][3]['DATA'][goodDatasets[0][1]]['type'])
		print('Dataset range:',goodDatasets[0][2][0],'um -',goodDatasets[0][2][1],'um')
		if (contains_k == goodDatasets):
			if(goodDatasets[0][1] != goodDatasets[0][4]):
				print('Dataset type:',goodDatasets[0][3]['DATA'][goodDatasets[0][4]]['type'])
				print('Dataset range:',goodDatasets[0][5][0],'um -',goodDatasets[0][5][1],'um')
		self.dataN = goodDatasets[0][1]
		self.dataK = goodDatasets[0][4]
		self.datafile = goodDatasets[0][3]

		# If k dataset is not found error is raised if getK is called
		if self.dataK == None:
			def error(l):
				raise Exception('No data for k was in dataset.')
			self.getEc = error
		# Check for the different types of data k and n might be contained in
		elif 'tabulated k' in self.datafile['DATA'][self.dataK]['type']:
			tabulated = _np.loadtxt(_io.StringIO(self.datafile['DATA'][self.dataK]['data']))
			self.getEc = _interp1d(tabulated[:,0],tabulated[:,1])
		if 'tabulated nk' in self.datafile['DATA'][self.dataN]['type']:
			tabulated = _np.loadtxt(_io.StringIO(self.datafile['DATA'][self.dataN]['data']))
			self.getN = _interp1d(tabulated[:,0],tabulated[:,1])
			self.getEc = _interp1d(tabulated[:,0],tabulated[:,2])
		elif 'tabulated n' in self.datafile['DATA'][self.dataN]['type']:
			tabulated = _np.loadtxt(_io.StringIO(self.datafile['DATA'][self.dataN]['data']))
			self.getN = _interp1d(tabulated[:,0],tabulated[:,1])
		else:
			formulaN = int(_re.search('(?<=formula\ )[0-9]',self.datafile['DATA'][self.dataN]['type']).group(0))
			coefficients =  _np.loadtxt(_io.StringIO(self.datafile['DATA'][self.dataN]['coefficients']))
			self.getN = MaterialTable.formulae[formulaN-1](coefficients,goodDatasets[0][2])

	# Additional function for error handling
	@staticmethod
	def isWithinRange(l,range):
		'''Checks whether the wavelength is within range'''
		l = _np.array(l)
		if not ((range[0] < l).all() and (range[1] > l).all()):
			raise ValueError('Some of the values are out of dataset range')

	# The different models utilised in the database for calculation of n
	@staticmethod
	def sellmeier(c,range=[-_np.inf,_np.inf]):
		'''The Sellmier dispersion formula'''
		if len(c) > 17:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(17-len(c)):
			cp = _np.append(cp,0)
		c = cp

		def func(l):
			'''The Sellmier dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = 1 + c[0]

			for i in _np.arange(8):
				n += c[2*i+1]*l**2/(l**2 - c[2*i+2]**2)
			return _np.sqrt(n)

		return func

	@staticmethod
	def sellmeier2(c,range=[-_np.inf,_np.inf]):
		'''The Sellmier2 dispersion formula'''
		if len(c) > 17:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(17-len(c)):
			cp = _np.append(cp,0)
		c = cp

		def func(l):
			'''The Sellmier2 dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = 1 + c[0]

			for i in _np.arange(8):
				n += c[2*i+1]*l**2/(l**2 - c[2*i+2])
			return _np.sqrt(n)

		return func

	@staticmethod
	def polynomial(c,range=[-_np.inf,_np.inf]):
		'''The Polynomial dispersion formula'''
		if len(c) > 17:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(17-len(c)):
			cp = _np.append(cp,0)
		c = cp
		def func(l):
			'''The Polynomial dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = c[0]

			for i in _np.arange(8):
				n += c[2*i+1]*l**c[2*i+2]
			return _np.sqrt(n)
		return func

	@staticmethod
	def refractiveIndexInfo(c,range=[-_np.inf,_np.inf]):
		'''The RefractiveIndex.INFO dispersion formula'''
		if len(c) > 17:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(17-len(c)):
			cp = _np.append(cp,0)
		c = cp

		def func(l):
			'''The RefractiveIndex.INFO dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = c[0]

			for i in _np.arange(2):
				n += c[4*i+1]*l**c[4*i+2]/(l**2 - c[4*i+3]**c[4*i+4])
			for i in _np.arange(4,8):
				n += c[2*i+1]*l**c[2*i+2]
			return _np.sqrt(n)
		return func

	@staticmethod
	def cauchy(c,range=[-_np.inf,_np.inf]):
		'''The Cauchy dispersion formula'''
		if len(c) > 11:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(11-len(c)):
			cp = _np.append(cp,0)
		c = cp
		def func(l):
			'''The Cauchy dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = c[0]

			for i in _np.arange(5):
				n += c[2*i+1]*l**c[2*i+2]
			return n
		return func

	@staticmethod
	def gases(c,range=[-_np.inf,_np.inf]):
		'''The Gases dispersion formula'''
		if len(c) > 11:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(11-len(c)):
			cp = _np.append(cp,0)
		c = cp
		def func(l):
			'''The Gases dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = c[0] + 1

			for i in _np.arange(5):
				n += c[2*i+1]/(c[2*i+2]-(1/l**2))
			return n
		return func

	@staticmethod
	def herzberger(c,range=[-_np.inf,_np.inf]):
		'''The Herzberger dispersion formula'''
		if len(c) > 6:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(6-len(c)):
			cp = _np.append(cp,0)
		c = cp

		def func(l):
			'''The Herzberger dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = c[0]
			n += c[1]*(1/(l**2 - 0.028))
			n += c[2]*(1/(l**2 - 0.028))**2
			n += c[3]*l**2 + c[4]*l**4 + c[5]*l**6
			return n
		return func

	@staticmethod
	def retro(c,range=[-_np.inf,_np.inf]):
		'''The Retro dispersion formula'''
		if len(c) > 4:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(4-len(c)):
			cp = _np.append(cp,0)
		c = cp

		def func(l):
			'''The Retro dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = c[0]
			n += c[1]*l**2/(l**2 - c[2])
			n += c[3]*l**2
			return _np.sqrt((2*n+1)/(1-n))
		return func

	@staticmethod
	def exotic(c,range=[-_np.inf,_np.inf]):
		'''The Exotic dispersion formula'''
		if len(c) > 6:
			print('Input array longer than expected parameter number, so it will be truncated')
		c = _np.array(c)
		cp = c
		for i in _np.arange(6-len(c)):
			cp = _np.append(cp,0)
		c = cp

		def func(l):
			'''The Exotic dispersion formula'''
			l = _np.array(l)
			MaterialTable.isWithinRange(l,range)
			n = c[0]
			n += c[1]/(l**2 - c[2])
			n += c[3]*(l-c[4])/((l-c[4])**2+c[5])
			return _np.sqrt(n)
		return func

# The formulae in list form to help parsing
MaterialTable.formulae = [MaterialTable.sellmeier,
						  MaterialTable.sellmeier2,
						  MaterialTable.polynomial,
						  MaterialTable.refractiveIndexInfo,
						  MaterialTable.cauchy,
						  MaterialTable.gases,
						  MaterialTable.herzberger,
						  MaterialTable.retro,
						  MaterialTable.exotic]

# Running init code on import
MaterialTable.initDatabase()

if __name__ == '__main__':

	fails = 0
	tests = 0
	def test(x,e):
		global fails,tests
		'''Function to wrap the unit test results. 1e-10 for rounding errors.
		x is the measured value while e is the expected.'''
		tests += 1
		if _np.abs(x-e)  >= 1e-10:
			print('\33[31mFAILED\33[0m')
			fails += 1
		else:
			print('\33[32mPASSED\33[0m')


	manganese = MaterialTable('Mn')
	# Unit test for the formula. The 1e-10 comparison is for floating point rounding errors.
	print('\nUnit testing the formulae:')
	print('1 Sellmeier: ',end='')
	test(MaterialTable.sellmeier([1]*17)(2),_np.sqrt(2+8*4/3))
	print('2 Sellmeier-2: ',end='')
	test(MaterialTable.sellmeier2([1]*17)(2),_np.sqrt(2+8*4/3))
	print('3 Polynomial: ',end='')
	test(MaterialTable.polynomial([1]*17)(2),_np.sqrt(1+2*8))
	print('4 refractiveIndex.INFO: ',end='')
	test(MaterialTable.refractiveIndexInfo([1]*17)(2),_np.sqrt(1+2*4+2*2/3))
	print('5 Cauchy: ',end='')
	test(MaterialTable.cauchy([1]*11)(2),11)
	print('6 Gases: ',end='')
	test(MaterialTable.gases([1]*11)(2),2+20/3)
	print('7 Herzberger: ',end='')
	test(MaterialTable.herzberger([1]*6)(2),1+(1/(4-0.028))+(1/(4-0.028))**2+4+16+64)
	print('8 Retro: ',end='')
	test(MaterialTable.retro([0.1]*4)(2),_np.sqrt(172/31))
	print('9 Exotic: ',end='')
	test(MaterialTable.exotic([1]*6)(2),_np.sqrt(11/6))

	print('Unit testing custom material: Vacuum')
	vacuum = MaterialTable('Vacuum')
	print('Testing n: ',end='')
	test(vacuum.getN(2),1)
	print('Testing Ec: ',end='')
	test(vacuum.getEc(2),0)

	if fails == 0:
		print('All tests passed')
	else:
		print(fails,'tests failed out of',tests)
