#!/usr/bin/ipython3
import yaml as _yml
import io as _io
import numpy as _np
import re as _re
from scipy.interpolate import interp1d as _interp1d
import scipy.constants as _const

class MaterialTable:
	'''This class initialises the database, while offering access to specific
	interpolations on the data. getRealN(wavelength) returns the n for a specific
	wavelength, while getEc(wavelength) returns the extinction coeffictent.
	In both cases wavelength is in units of um'''

	# The object that will contain the datafiles' location (i.e. the database manifest)
	library = None

	@staticmethod
	def initDatabase(reload=False):
		'''This method loads the map for the database'''

		if MaterialTable.library != None and reload:
			print("Reloading database mapping file.")
		elif MaterialTable.library != None:
			return # exit if reload is not required.

		try:
			with open(r'./database/parsed_lib.yml') as file:

				# Loading yml mapping.
				MaterialTable.library = _yml.full_load(file)

		except FileNotFoundError:
			print('Library file at ./database/parsed_lib.yml was not found')

	@classmethod
	def fromYml(cls,fname):
		'''Loads a yml file with the same structure as the ones contained in the
		refractiveindex.info database.'''
		table = cls() #initialise an empty MaterialTable named table

		with open(fname) as file:
			table.datafile = _yml.full_load(file)
			print('Loaded yml datafile:',fname)

		dataK_found = False
		dataN_found = False
		# Find k and n data and load them
		for datatable in table.datafile['DATA']:
			if datatable['type'] == 'tabulated nk':
				# If both found in same dataset no extra data is needed
				data = _np.loadtxt(_io.StringIO(datatable['data']))
				print('Loaded tabulated nk data.')
				print('Range:',data[0,0],'um -',data[-1,0],'um')
				if not dataN_found:
					table.getRealN = _interp1d(data[:,0],data[:,1])
				else:
					print ('\33[33mWarning:\33[0m Duplicate n data. Skipped')
				if not dataK_found:
					table.getEc = _interp1d(data[:,0],data[:,2])
				else:
					print ('\33[33mWarning:\33[0m Duplicate k data. Skipped')

				dataN_found = True
				dataK_found = True

			elif datatable['type'] == 'tabulated n':
				if dataN_found:
					print ('\33[33mWarning:\33[0m Found duplicate n data. Skipped')
					continue

				data = _np.loadtxt(_io.StringIO(datatable['data']))
				print('Loaded tabulated n data')
				print('Range:',data[0,0],'um -',data[-1,0],'um')

				dataN_found = True
				table.getRealN = _interp1d(data[:,0],data[:,1])
			elif datatable['type'] == 'tabulated k':
				if dataK_found:
					print ('\33[33mWarning:\33[0m Found duplicate k data. Skipped')
					continue

				data = _np.loadtxt(_io.StringIO(datatable['data']))
				print('Loaded tabulated k data')
				print('Range:',data[0,0],'um -',data[-1,0],'um')

				dataK_found = True
				table.getEc = _interp1d(data[:,0],data[:,1])

			elif 'formula' in datatable['type']:
				if dataN_found:
					print ('\33[33mWarning:\33[0m Found duplicate n data. Skipped')
					continue

				ID = int(_re.search('(?<=formula\ )[0-9]',datatable['type']).group(0)) - 1

				if (ID >= len(MaterialTable.formulae)):
					raise KeyError('Unknown formula:',datatable['type'])

				print('Loaded',datatable['type'],'data')
				range = _np.loadtxt(_io.StringIO(datatable['wavelength_range']))
				coeffs = _np.loadtxt(_io.StringIO(datatable['coefficients']))
				print('Range:',range[0],'um -',range[1],'um')
				table.getRealN = MaterialTable.formulae[ID](coeffs,range)

				dataN_found = True

			else:
				raise ValueError('Unknown type of datatable:',datatable['type'])


			if dataN_found and dataK_found:
				# Check if we are at the end of the dataset
				if datatable != table.datafile['DATA'][-1]:
					print('\33[33mWarning:\33[0m Some datatables were skipped \
					as sufficient data was found before end of file.')
				break # End loop if all parameters are loaded

		if not dataK_found:
			print('No k data found. Assuming lossless material')
			table.getEc = lambda l : 0

		return table


	@classmethod
	def fromValue(cls,realN,Ec=0):
		'''Constructor that allows the creation of a theoretical material
		with complex refractive index realN + 1j*Ec.'''
		material = cls()
		material.getRealN = lambda l : realN
		material.getEc = lambda l : Ec
		print('Created material with n =',n,'and k =',Ec)
		return material

	@classmethod
	def fromCsv(cls,fname,skiprows=0,delimiter=',',wavelength_multiplier=1):
		'''A constructor that allows the loading of a csv datafile.'''

		material = cls()

		data = _np.loadtxt(fname,skiprows=skiprows,delimiter=delimiter)
		data[:,0]*=wavelength_multiplier

		print('Loaded',fname)
		print('Dataset range:',data[0,0],'um -',data[-1,0],'um')

		if len(data.shape) != 2:
			raise ValueError('Loaded data was misshapen.')
		if data.shape[1] > 3 :
			print('\33[33mWarning:\33[0m Extraneous columns detected. These were \
			skipped')
		if data.shape[1] == 2:
			print('2 columns detected. Assuming lossless material.')
			material.getEc = lambda l : 0
			material.getRealN = _interp1d(data[:,0],data[:,1])
		else:
			material.getRealN = _interp1d(data[:,0],data[:,1])
			material.getEc = _interp1d(data[:,0],data[:,2])
		return material

	@classmethod
	def fromMaterial(cls,material,lmin=None,lmax=None,allow_interpolation=True,\
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

		materialTable = cls()
		materialTable.findMaterial(material,lmin,lmax,allow_interpolation,require_k)
		return materialTable

	def __init__(self):
		'''The default constructor should not be called. Use fromMaterial(material,
		lmin=None,lmax=None,allow_interpolation=True,require_k=False). '''
		MaterialTable.initDatabase()

	def getRealN(self,l):
		'''Default function for returning real refractive index. As it was not
		initialised yet it only throws an error. This function is replaced
		by proper initialisation.'''
		raise Exception('The material table was not properly initialised')

	def getEc(self,l):
		'''Default function for returning extinction coefficient. As it was not
		initialised yet it only throws an error. This function is replaced
		by proper initialisation.'''
		raise Exception('The material table was not properly initialised')

	def getN(self,l):
		'''Method to get the complex refractive index at a given wavelength l.'''
		return self.getRealN(l) + 1j*self.getEc(l)

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


		if material == None or material == 'Vacuum':
			print('Loading custom material: Vacuum')
			# To remain consistent with instance variables
			self.datafile = None
			dataN = None
			dataK = None
			self.getRealN = lambda l : 1
			self.getEc = lambda l : 0
			return # call return so search would not be called

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
		print('Datatable type:',goodDatasets[0][3]['DATA'][goodDatasets[0][1]]['type'])
		print('Datatable range:',goodDatasets[0][2][0],'um -',goodDatasets[0][2][1],'um')
		if (contains_k == goodDatasets):
			if(goodDatasets[0][1] != goodDatasets[0][4]):
				print('Datatable type:',goodDatasets[0][3]['DATA'][goodDatasets[0][4]]['type'])
				print('Datatable range:',goodDatasets[0][5][0],'um -',goodDatasets[0][5][1],'um')
		dataN = goodDatasets[0][1]
		dataK = goodDatasets[0][4]
		self.datafile = goodDatasets[0][3]

		# If k dataset is not found error is raised if getK is called
		if dataK == None:
			self.getEc = lambda l : 0
		# Check for the different types of data k and n might be contained in
		elif 'tabulated k' in self.datafile['DATA'][dataK]['type']:
			tabulated = _np.loadtxt(_io.StringIO(self.datafile['DATA'][dataK]['data']))
			self.getEc = _interp1d(tabulated[:,0],tabulated[:,1])
		if 'tabulated nk' in self.datafile['DATA'][dataN]['type']:
			tabulated = _np.loadtxt(_io.StringIO(self.datafile['DATA'][dataN]['data']))
			self.getRealN = _interp1d(tabulated[:,0],tabulated[:,1])
			self.getEc = _interp1d(tabulated[:,0],tabulated[:,2])
		elif 'tabulated n' in self.datafile['DATA'][dataN]['type']:
			tabulated = _np.loadtxt(_io.StringIO(self.datafile['DATA'][dataN]['data']))
			self.getRealN = _interp1d(tabulated[:,0],tabulated[:,1])
		else:
			formulaN = int(_re.search('(?<=formula\ )[0-9]',self.datafile['DATA'][dataN]['type']).group(0))
			coefficients =  _np.loadtxt(_io.StringIO(self.datafile['DATA'][dataN]['coefficients']))
			self.getRealN = MaterialTable.formulae[formulaN-1](coefficients,goodDatasets[0][2])

	# Additional function for error handling
	@staticmethod
	def isWithinRange(l,range):
		'''Checks whether the wavelength is within range. A supporting method
		for easier error handling in the different datamodels.'''
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

def getWavenumber(l):
	'''A function for transforming wavelength in um to wavenumber
	in um^-1'''
	return 2*_const.pi/l

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


	BK7 = MaterialTable.fromMaterial('BK7')
	BK7_Schott = MaterialTable.fromYml('./database/data/glass/schott/N-BK7.yml')
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

	print('\nUnit testing custom material: Vacuum')
	vacuum = MaterialTable.fromMaterial('Vacuum')
	print('Testing n: ',end='')
	test(vacuum.getRealN(2),1)
	print('Testing Ec: ',end='')
	test(vacuum.getEc(2),0)

	print('\nUnit testing fromYml using BK7 (SCHOTT):')
	wavelengths = [0.5,0.6,0.7,0.8,0.9,1]
	print('Testing complex n at [0.5,0.6,0.7,0.8,0.9,1]: ',end='')
	test(BK7.getN(wavelengths),BK7_Schott.getN(wavelengths))

	if fails == 0:
		print('All tests passed')
	else:
		print(fails,'tests failed out of',tests)
