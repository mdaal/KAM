'''This is the KIDs Analysis Module, KAM'''

import urllib2
import scipy.io #for loading .mat file
import os
import numpy as np


import tables
import matplotlib.pyplot as plt
import datetime
from scipy.optimize import minimize, leastsq

import numpy.ma as ma
import sys # for status percentage

import platform
mysys = platform.system()

import warnings #trying to get a warning every time rather than just the first time.
warnings.filterwarnings('always')

database_location = 'Data' + os.sep + 'My_Data_Library.h5'


try:
	execfile('KIPs_Access.txt')
except:
	print('KIPs_Access.txt not found. Create this file to download data.')
	# remote access file must have these two lines
	# username = _________  # e.g. 'johndoe'
	# password = _________  # e.g. '294hr5'
	##############

class loop:
	'''The purpose of this class is to hold data associated with resonance loops and fitting them'''
	def __init__(self):
		self.index = None
		self.z =  None
		self.freq = None

		#output of circle fit
		self.r = None
		self.a = None #fit circle center is a+i*b
		self.b = None
		self.outer_iterations = None
		self.inner_iterations = None

		#circle fit parameters
		self.s = None
		self.Gx = None
		self.Gy = None
		self.g = None
		self.sigma = None
		self.circle_fit_exit_code = None

		#loop fit estimates quantities
		self.fr_est = None
		self.FWHM_est = None
		self.depth_est = None
		self.Q_est = None

		#fit quantities
		self.Q = None 
		self.Qc = None
		self.Qi = None
		self.fr = None 
		self.FWHM = None
		self.phi = None
		self.chisquare = None
		self.pvalue = None
		self.Phase_Fit_Method = None

	def __del__(self):
		pass
		

class metadata:
	def __init__(self):
		self.Time_Created = None
		self.Atten_Added_At_NA = None # redundant if self.Atten_NA_Input and self.Atten_NA_Output are defined; should be merged somehow
		self.NA_Average_Factor = None
		self.Fridge_Base_Temp = None
		self.Box = None
		self.Ground_Plane = None
		self.LNA = None
		self.IFBW = None
		self.Test_Location = None
		self.Minimum_Q = None
		self.Notes = None
		self.Num_Points_Per_Scan = None
		self.Wait_Time = None
		self.Press = None
		self.Min_Freq_Resolution = None
		self.Run = None
		self.Sensor = None
		self.Fridge_Run_Start_Date = None
		self.Fsteps  = None
		#self.IS_Sonnet_Simulation = None
		self.Data_Source = None
		self.Atten_At_4K = None
		self.Atten_NA_Output = None # positive value in dB
		self.Atten_NA_Input = None # positive value in dB
		self.Atten_RTAmp_Input = None # positive value in dB
		self.RTAmp_In_Use = None
		self.Meaurement_Duration = None
		self.Num_Heater_Voltages = None
		self.Num_Powers = None
		self.Num_Ranges = None 
		self.Cable_Calibration = None
		self.Temperature_Calibration = None # a list of tuples [(heatervoltge1, temperature), (heatervoltage2,temperature, ...)]

class thermometry:
	def __init__(self):
		pass
	def load_MonitoringVI_file(self, filename):
		import io
		from scipy.signal import gaussian,wiener, filtfilt, butter,  freqz
		from scipy.ndimage import filters
		from scipy.interpolate import UnivariateSpline
		pos = filename.rfind(os.sep)

		try:
			with io.open(filename[:pos+1]+ 'Make_ScanData.m',mode='r') as f:
				while 1:
					line  = f.readline()
					if line == '': # End of file is reached
						break
					elif line.find('ScanData.Heater_Voltage') >= 0:
						Voltages = line[line.find('['):line.find(']')+1]
						break
		except:
			print('Unable to find or read Make_ScanData.m for list of heater voltages')
			Voltages = 'Unknown'


		temp_data = np.loadtxt(filename, dtype=np.float, comments='#', delimiter=None, converters=None, skiprows=3, usecols=None, unpack=False, ndmin=0)
		
		num_pts_in_gaussian_window = 20
		b = gaussian(num_pts_in_gaussian_window, 10)
		ga = filters.convolve1d(temp_data[:,1], b/b.sum())

		npts = temp_data[:,1].size
		end = temp_data[-1,0]
		dt = end/float(npts)
		nyf = 0.5/dt	
		b, a = butter(4, .1)#1.5/nyf)
		fl = filtfilt(b, a, temp_data[:,1])
    	
		sp = UnivariateSpline(temp_data[:,0], temp_data[:,1])

		wi = wiener(temp_data[:,1], mysize=40, noise=10)

		fig1 = plt.figure( facecolor = 'w',figsize = (10,10))
		ax = fig1.add_subplot(1,1,1)		
		try:
			line = ax.plot(temp_data[:,0], temp_data[:,2],'g:', linewidth = 3,label = 'Data2')
		except:
			pass
		line = ax.plot(temp_data[:,0], temp_data[:,1],'g', linewidth = 3,label = 'Data')
		line2 = ax.plot(temp_data[:,0], ga, 'y', linewidth = 3, label = 'Gaussian Conv') # Gaussian Convolution
		line3 = ax.plot(temp_data[:,0], fl, 'c', linewidth = 3, label = 'Butterworth') # butterworth 
		line4 = ax.plot(temp_data[:,0], sp(temp_data[:,0]), 'k', linewidth = 3, label = 'Spline') # bspline
		line5 = ax.plot(temp_data[:,0], wi, 'r', linewidth = 3, label = 'Weiner') # weiner

		ax.grid(b=True, which='major', color='b', linestyle='-')
		ax.grid(b=True, which='minor', color='b', linestyle='--')
		ax.set_title('Heater Voltages = {}'.format(Voltages), fontsize=12)
		ax.set_ylabel('Temperature [Kelvin]')
		ax.set_xlabel('Seconds')
		ax.legend(loc = 'best', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
		plt.show()
		

class sweep:
	'''This class accesses resonance data and fits resonances'''

	
	
	def __init__(self):
		self.loop = loop()
		self.metadata = metadata()
		
		self.data_set_contents = np.dtype([
			("Run"				, 'S10'),
			("Time_Created"		, 'S40'), # 'S40' for format December 23, 2012 12:34:65.675 PM;  'S12' for format '%Y%m%d%H%M' 
			("Num_Ranges"		, np.uint8), # uint8 is Unsigned integer (0 to 255)
			("Num_Powers"		, np.uint8),
			("Num_Temperatures"	, np.uint8), 
			("Sensor"			, 'S20'),
			("Ground_Plane"		, 'S20'),
			("Path"				, 'S100'),
			])

	def _read_scandata_from_file(self,filename_or_path):
		
		# index = filename_or_path.rfind(os.sep)
		# if index > -1: # filename_or_path is a path
		# 	current_path = os.getcwd()
		# 	os.chdir(filename_or_path[0:index])
		# 	mat = scipy.io.loadmat(filename_or_path[index+1:])
		# 	os.chdir(current_path)
		# else: # filename_or_path is a filename
		# 	mat = scipy.io.loadmat(filename_or_path)

		mat = scipy.io.loadmat(filename_or_path)
		self.data = mat
		self.metadata.Data_Source = filename_or_path

	def _download_data(self, URL):
		''' Authenticats to URL containing data.
		Copies the .mat file licated at URL to a local file in local directory.
		.mat file is a Scan_Data matlab structure.
		returns numpy data structure contauning .mat file.
		deletes local file.'''


		passman = urllib2.HTTPPasswordMgrWithDefaultRealm()
		# this creates a password manager
		passman.add_password(None, URL, username, password)
		# because we have put None at the start it will always
		# use this username/password combination for  urls
		# for which `URL` is a super-url

		authhandler = urllib2.HTTPBasicAuthHandler(passman)
		# create the AuthHandler

		opener = urllib2.build_opener(authhandler)

		urllib2.install_opener(opener)
		# All calls to urllib2.urlopen will now use our handler
		# Make sure not to include the protocol in with the URL, or
		# HTTPPasswordMgrWithDefaultRealm will be very confused.
		# You must (of course) use it when fetching the page though.

		pagehandle = urllib2.urlopen(URL)
		# authentication is now handled automatically for us

		#import tempfile # Attempt to download data into a temp file
		#f = tempfile.NamedTemporaryFile(delete=False)
		#f.write(pagehandle.read())
		#f.close()
		#mat = scipy.io.loadmat(f.name)

		output = open('test.mat','wb')
		print('Download Initiated...')
		output.write(pagehandle.read())
		print('Download Completed...')
		output.close()
		#global mat
		mat = scipy.io.loadmat('test.mat')

		#this id how to tell what variables are stored in test.mat
		#print scipy.io.whosmat('test.mat')



		#html = pagehandle.read()
		#pagehandle.close()

		#soup = BeautifulSoup(html)
		#soup.contents
		os.remove('test.mat')
		self.data = mat
		self.metadata.Data_Source = URL

	def plot_loop(self,  aspect='equal', show = True):
		''' Plots currently selected complex transmission in the I,Q plane. Reutrns a tuple, (fig, ax, line),
		where  fig is the figure object, ax is the axes object and line is the line object for the plotted data.

		aspect='equal' makes circles round, aspect='auto' fills the figure space.

		*Must have a loop picked in order to use this function.*
		'''
		try: 
			z = self.loop.z
		except:
			print("Data not available. You probably forgot to load it.")
			return


		fig = plt.figure( figsize=(8, 6), dpi=100)
		ax = fig.add_subplot(111)#,aspect=aspect)
		line = ax.plot(z.real,z.imag,'bo')
		ax.set_xlabel('I [Volts]')
		ax.set_ylabel('Q [Volts]')
		ax.set_title('Run: {0}; Sensor: {1}; Ground: {2}; Record Date: {3}'.format(self.metadata.Run, self.metadata.Sensor, self.metadata.Ground_Plane, self.metadata.Time_Created),fontsize=10)
		if show == True:
			plt.show()
		return  (fig, ax, line)

	def plot_transmission(self, show = True):
		''' Plots currently selected complex transmission in dB as a function of frequency. Reutrns a tuple, (fig, ax, line),
		where fig is the figure object, ax is the axes object and line is the line object for the plotted data.

		*Must have a loop picked in order to use this function.*
		'''
		try: 
			z = self.loop.z
			freq = self.loop.freq
		except:
			print("Data not available. You probably forgot to load it.")
			return

		plt.rcParams["axes.titlesize"] = 10
		fig = plt.figure( figsize=(8, 6), dpi=100)
		ax = fig.add_subplot(111)
		line = ax.plot(freq,20*np.log10(abs(z)),'b-',)
		ax.set_xlabel('Frequency [Hz]')
		ax.set_ylabel('$20*Log_{10}[|S_{21}|]$ [dB]')

		ax.set_title('Run: {0}; Sensor: {1}; Ground: {2}; Record Date: {3}'.format(self.metadata.Run, self.metadata.Sensor, self.metadata.Ground_Plane, self.metadata.Time_Created))
		if show == True:
			plt.show()
		return  (fig, ax, line)

	def _extract_type(self, obj, return_type = None, field = None):
		'''scanandata object, obj, has a lot of single element arrays of arrays. this function gets the element.
		e.g scandata may have [[[ele]]] instead of callling ele = scandata[0][0][0], use this function to get ele.
		if ele is another structured numpy array with field name 'myfield', using keyword field = 'myfield' will get
		the data at field.
		the function will cast ele to be in the data type return_typer. e.g. return_type = 'str' returns a string. 
		If return_type is None, ele is returned as whatever type it was saved as in [[[ele]]] '''
		
		def cast(_obj):
			if (return_type != None) & (_obj != None) :
				_obj = return_type(_obj)
				#pass#exec("_obj = {0}(_obj)".format(return_type))
			return _obj

		def itemloop(_obj):
			while True:
				try:
					_obj = _obj.item()
				except:
					return cast(_obj)
			return cast(_obj)


		if field == None:
			obj = itemloop(obj)

		else:
			while obj.dtype == np.dtype('O'):
				obj = obj.item()
			try:
				obj = obj[field]
			except:
				obj = None
				print('Field named {0} is not found. Returning None'.format(field))
			obj = itemloop(obj)
		return obj

	def _define_sweep_data_columns(self, fsteps):
		self.metadata.Fsteps = fsteps

		self.sweep_data_columns = np.dtype([
			("Fstart"         	, np.float32), # in Hz
			("Fstop"          	, np.float32), # in Hz
			("Heater_Voltage" 	, np.float32), # in Volts
			("Pinput_dB"      	, np.float32), # in dB
			("Preadout_dB"     	, np.float32), # in dB  - The power at the input of the resonator, not inside the resonator
			("Temperature"    	, np.float32), # in Kelvin
			("S21"            	, np.complex128, (fsteps,)), # in complex numbers
			("Frequencies"    	, np.float32,(fsteps,)), # in Hz
			("Q"				, np.float32),
			("Qc"				, np.float32),
			("fr"				, np.float32), # in Hz

			])

	def _define_sweep_array(self,index,**field_names):
		#for field_name in self.sweep_data_columns.fields.keys():
		for field_name in field_names:
			self.Sweep_Array[field_name][index] = field_names[field_name]
			
	def load_scandata(self, file_location):
		''' file_location is the locaiton of the scandata.mat file. It can be a URL, filename or /path/filename.
		assumes that self.data is in the form of matlab ScanData Structure'''

		#delete previous metadata object
		del(self.metadata)
		self.metadata = metadata()

		if file_location.startswith('http'):
			self._download_data(file_location)
		else:
			self._read_scandata_from_file(file_location)

		ScanData = self.data['ScanData']
		
		# These tags specify the data to pull out of self.data['ScanData']. syntax is 
		# (field of self.data['ScanData'] to extract, self.metadata name to save to ('key:sub-key' ifself.metadata.key is a dict), 
		#			 type of value (arrays are None),optional sub-field of self.data['ScanData'] to extract)
		tags = [('Run','Run', str), ('Start_Date','Fridge_Run_Start_Date',str), ('Location','Test_Location', str), 
				('Sensor','Sensor',str), ('Ground_Plane','Ground_Plane',str), ('Box','Box',str), ('Press','Press',str), 
				('Notes','Notes',str),('Time','Time_Created',str),('Temperature','Fridge_Base_Temp',float),
				('Powers','Powers', None), ('Resolution','Min_Freq_Resolution', np.float), ('IFBW','IFBW', np.float),
				('Heater_Voltage','Heater_Voltage',None), ('Average_Factor','NA_Average_Factor', np.float), 
				('Minimum_Q','Minimum_Q', np.float), ('Range','Range',None), ('Added_Atten','Atten_Added_At_NA', np.float),  
				('Num_Points_Per_Scan','Num_Points_Per_Scan',np.float), ('Freq_Range', 'Freq_Range',None), 
				('Pause','Wait_Time',np.float), ('LNA', 'LNA:LNA', str), ('HEMT', 'LNA:Vg', str,'Vg'),
				('HEMT', 'LNA:Id', str,'Id'),  ('HEMT', 'LNA:Vd', str,'Vd'), ('Atten_4K', 'Atten_At_4K', np.float32),
				('Atten_NA_Output', 'Atten_NA_Output',np.float32), ('Atten_NA_Input','Atten_NA_Input',np.float32),
				('Atten_RTAmp_Input','Atten_RTAmp_Input',np.float32), ('RTAmp_In_Use', 'RTAmp_In_Use', int),
				('Elapsed_Time', 'Meaurement_Duration', np.float)]

		for t in tags:
			try:
				if t[1].find(':')>-1:
					t1 = t[1].split(':')

					#This try/except block is for the case where self.metadata.__dict__['?'] is a dictionary
					try:
						self.metadata.__dict__[t1[0]].update([(t1[1],self._extract_type(ScanData[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None))])
					except:
						self.metadata.__dict__[t1[0]] = dict([(t1[1],self._extract_type(ScanData[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None))])	
				else:
					self.metadata.__dict__[t[1]] = self._extract_type(ScanData[t[0]], return_type = t[2],field = t[3] if len(t) > 3 else None)
			except: 
				#the case that the field does not exist or that its in an unexpected format
				print('Field named {0}{1} is not found. Setting value to None'.format(t[0], (':'+t[3]) if len(t) > 3 else ''))

		try:
			self.metadata.Powers                = self.metadata.Powers.squeeze() #for case there are multiples powers
		except:
			self.metadata.Powers                = np.array([self.metadata.Powers]) # for the case there is only one power		

		# Reshape  Heater_Voltage array and  Remove final Heater voltage from self.Heater_Voltage (The final value does just the heater value at which to leave fridge )
		self.metadata.Heater_Voltage = self.metadata.Heater_Voltage.reshape((self.metadata.Heater_Voltage.shape[1],))
		self.metadata.Heater_Voltage = self.metadata.Heater_Voltage[0:-1]

		# Determine lenght of S21 and Frequencies. Note in the next line how self.Freq_Range is structured. 
		print('There are {0} heater voltage(s), {1} input power(s), and {2} frequecy span(s)'.format(self.metadata.Heater_Voltage.shape[0],self.metadata.Powers.shape[0], self.metadata.Freq_Range.shape[0]))
		heater_voltage_num = 0; power_sweep_num = 0; fsteps = 0;

		# determin fsteps, the length of the freq/S21 array
		if self.metadata.Heater_Voltage.shape[0] == 1:
			fsteps = self.metadata.Freq_Range[heater_voltage_num][1]['PowerSweep'][0][0][power_sweep_num][2].squeeze()[()].size # non temp sweep, single freq_range, powersweep
		else:					
			for freq_range_num in xrange(self.metadata.Freq_Range.shape[0]):			
				steps = self.metadata.Freq_Range[freq_range_num][1]['Temp'][0][0][heater_voltage_num][1]['PowerSweep'][0][0][power_sweep_num][2].squeeze()[()].size
				fsteps = max(steps,fsteps)

		self._define_sweep_data_columns(fsteps)

		self.metadata.Num_Powers= self.metadata.Powers.size
		self.metadata.Num_Heater_Voltages = self.metadata.Heater_Voltage.size
		self.metadata.Num_Ranges = self.metadata.Range.shape[0]
		try:
			self.metadata.Cable_Calibration = self._Cable_Calibration
			print('Cable Calibraiton data found and saved in Sweep_Array metadata.')
		except:
			pass

		try:
			self.metadata.Temperature_Calibration = self._Temperature_Calibration
			print('Temperature Calibraiton data found and saved in Sweep_Array metadata.')
		except:
			pass
		### Examples of dealing with Freq_Range Data structure imported from Matlab .mat file			
		#    k.Freq_Range[heater_voltage_num][1]['PowerSweep']
		# 					  k.Freq_Range[0][1]['PowerSweep']
		# j.Freq_Range[0][1]['Temp'][0][0][0][1]['PowerSweep']
		# dt = np.dtype(('O', (2,3)))
		# entended = np.zeros(0,dtype = dt)
		# dt = np.dtype(('O',('O',[('Temp',('O',('O')))])))
		# dt = np.dtype([('O',[('O',[('Temp',[('O',('O'))])])])])
		# #this is the closest I can come to replecating the structure of a Temp Power Sweep
		# dt = np.dtype([('O',[('Temp','O',(1,1))],(1,2))])

		
		i=0
		self.Sweep_Array = np.zeros(self.metadata.Heater_Voltage.shape[0]*self.metadata.Powers.shape[0]*self.metadata.Freq_Range.shape[0], dtype = self.sweep_data_columns)
		for freq_range_num in xrange(self.metadata.Freq_Range.shape[0]):
				if self.metadata.Heater_Voltage.shape[0] == 1:
					heater_voltages = self.metadata.Freq_Range # non temp sweep, single freq_range, powersweep
				else:					
					heater_voltages = self._extract_type(self.metadata.Freq_Range[freq_range_num,1]['Temp'])
				#start here for single res powersweep
				for heater_voltage_num in xrange(heater_voltages.shape[0]):
					sweep_powers = self._extract_type(heater_voltages[heater_voltage_num,1], field = 'PowerSweep')
					for sweep in sweep_powers[:,0:sweep_powers.shape[1]]:
						self._define_sweep_array(i, Fstart = self.metadata.Range[freq_range_num,0],
													Fstop = self.metadata.Range[freq_range_num,1],
													Heater_Voltage = self.metadata.Heater_Voltage[heater_voltage_num],
													Pinput_dB = sweep[0].squeeze()[()],
													S21 = sweep[1].squeeze()[()],
													Frequencies =  sweep[2].squeeze()[()])
						i = i + 1


		del(self.metadata.__dict__['Powers'])
		del(self.metadata.__dict__['Heater_Voltage'])
		del(self.metadata.__dict__['Range'])
		del(self.metadata.__dict__['Freq_Range'])

		# if self.metadata.Atten_NA_Output == None: #redundant to have both
		# 	del(self.metadata.__dict__['Atten_NA_Output'])
		# else:
		# 	del(self.metadata.__dict__['Atten_Added_At_NA'])

	def load_touchstone(self,filename, pick_loop = True):
		''' The function loads S21 and Freq from  Sonnet .s2p or .s3p files into the Sweep_Array structured np array
		All Sij are extracted, but only  S21 is saved into Sweep_Array. Future editions of this code might  find need 
		to load otherSij becuase S21.

		The function only loads one transmission array (S21).  pick_loop = True immediatly selectes this loop as the 
		current loop.
		'''

		import tempfile
		import io

		#delete previous metadata object
		del(self.metadata)
		self.metadata = metadata()

		dt_s2p = [('Freq', np.float64), ('S11r', np.float64), ('S11i', np.float64), ('S12r', np.float64), ('S12i', np.float64), 
										('S21r', np.float64), ('S21i', np.float64), ('S22r', np.float64), ('S22i', np.float64)]
		
		dt_s3p = [('Freq', np.float64), ('S11r', np.float64), ('S11i', np.float64), ('S12r', np.float64), ('S12i', np.float64), ('S13r', np.float64), ('S13i', np.float64),
										('S21r', np.float64), ('S21i', np.float64), ('S22r', np.float64), ('S22i', np.float64), ('S23r', np.float64), ('S23i', np.float64),
										('S31r', np.float64), ('S31i', np.float64), ('S32r', np.float64), ('S32i', np.float64), ('S33r', np.float64), ('S33i', np.float64)] 


		with tempfile.TemporaryFile() as tmp:
			with io.open(filename, mode='r') as f:
				# The following while loop copies the .sNp file into a temp file, which is destroyed when closed,
				# such that the tmp file is formated in the way np.loadtxt can read the data.
				indented = False
				prev_line = ''
				m = 1. # for frequency base conversion
				while 1: 
					line  = f.readline().replace('\n','')

					pos = f.tell()
					if line == '': # End of file is reached
						break
					elif line.startswith('! Data File Written:'): # Save as Metadata
						self.metadata.Time_Created = str(line.split('! Data File Written:')[1].strip())
						tmp.write(line + '\n')
					elif line.startswith('! From Project:') | line.startswith('! From Emgraph Data:'): # Save as Metadata
						self.metadata.Run = str(line.split(':')[1].strip())
						#self.metadata.IS_Sonnet_Simulation = True
						tmp.write(line + '\n')
					elif line[0] == '#':
						line  = line.replace('#','!#')
						if line.find('GHZ') >=-1:
							m = 1.0e9
						freq_convert = lambda s: s*m #Convert to Hertz
						tmp.write(line + '\n')	
					
					elif line[0] == ' ': # in S matrix definition block
						prev_line = prev_line + ' ' + line.strip() + ' '
						next_line = f.readline()
						# if next line is NOT indented date, then S matrix definition block is finished 
						# and we write it to tmp on a single line.
						# for .s2p files the S matrix is fully defined on one line of f
						# for .s3p files, the S matrix is defined in three lines. second two are indented.
						if not ((next_line[0] == '') | (next_line[0] == ' ')):
							tmp.write(prev_line)
							tmp.write('\n')
							prev_line = ''
						f.seek(pos,0)
		
					elif line[0] == '!':
						tmp.write(line + '\n')

					else:
						tmp.write(line)
						next_line = f.readline()
						# add \n to line if it does not begin a S matrix definition block
						if not ((next_line[0] == '') | (next_line[0] == ' ')):
							tmp.write('\n')
						f.seek(pos,0)

			tmp.seek(0,0)
			if filename.endswith('.s2p'):
				dt = dt_s2p
			elif filename.endswith('.s3p'):
				dt = dt_s3p	
			Touchstone_Data = np.loadtxt(tmp, dtype=dt, comments='!', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)
		
		self._define_sweep_data_columns(Touchstone_Data.size)
		j = np.complex(0,1)

		self.Sweep_Array = np.zeros(1, dtype = self.sweep_data_columns)
		
		self._define_sweep_array(0, Fstart = freq_convert(Touchstone_Data['Freq'].min()), #Hz
									Fstop = freq_convert(Touchstone_Data['Freq'].max()), #Hz
									S21 = Touchstone_Data['S21r']+j*Touchstone_Data['S21i'],
									Frequencies = freq_convert(Touchstone_Data['Freq']) #Hz
									)


		self.metadata.Data_Source = filename
		#self.metadata.Min_Freq_Resolution = np.abs(Touchstone_Data['Freq'][:-1]-Touchstone_Data['Freq'][1:]).min()
		self.metadata.Min_Freq_Resolution = np.abs(Touchstone_Data['Freq'][0] - Touchstone_Data['Freq'][-1])/Touchstone_Data['Freq'].size #use average freq resolution
	
		if pick_loop == True: #since there is only one loop in Sweep_Array, we might as well pick it as the current loop
			self.pick_loop(0)

	def downsample_loop(self,N):
		''' Reduce number of loop/freq data point by every Nth point and discarding all others'''
		self.loop.z = self.loop.z[0:-1:N]
		self.loop.freq = self.loop.freq[0:-1:N]

	def save_hf5(self, filename = database_location, overwrite = False):
		'''Saves current self.Sweep_Array into table contained in the hdf5 file speficied by filename.
		If overwite = True, self.Sweep_Array will overwright whatever is previous table data there is.
		'''
		
		if not os.path.isfile(filename):
			print('Speficied h5 database does not exist. Creating new one.')
			pos = filename.find('/')
			if pos >= 0:
				try:
					os.makedirs(filename[0:pos+1])
				except OSError:
					print('{0} exists...'.format(filename[0:pos+1]))
			wmode = 'w'
		else:
			print('Speficied h5 database exists and will be updated.')
			wmode = 'a'
			
		db_title = 'Aggregation of Selected Data Sets'
		group_name = 'Run' + self.metadata.Run
		group_title = self.metadata.Test_Location
		try:
			# case for scan data date
			d = datetime.datetime.strptime(self.metadata.Time_Created, '%B %d, %Y  %I:%M:%S.%f %p') # slightly wrong %f is microseconds. whereas we want milliseconds.
		except:
			pass
		try:
			#Case for sonnet date
			d = datetime.datetime.strptime(self.metadata.Time_Created, '%m/%d/%Y %H:%M:%S')
		except:
			pass
		sweep_data_table_name = 'T' + d.strftime('%Y%m%d%H%M')

		

		with tables.open_file(filename, mode = wmode, title = db_title ) as fileh:
			try:
				table_path = '/' + group_name + '/' + sweep_data_table_name
				sweep_data_table = fileh.get_node(table_path)

				if overwrite == True:
					print('Table {0} exists. Overwriting...'.format(table_path))
					sweep_data_table.remove()
					sweep_data_table = fileh.create_table('/'+ group_name,sweep_data_table_name,description=self.sweep_data_columns,title = 'Sweep Data Table',filters=tables.Filters(0), createparents=True)
				else:
					print('Table {0} exists. Aborting...'.format(table_path))
					return
			except:
				print('Creating table {0}'.format('/'+ group_name+'/'+sweep_data_table_name))
				sweep_data_table = fileh.create_table('/'+ group_name,sweep_data_table_name,description=self.sweep_data_columns,title = 'Sweep Data Table',filters=tables.Filters(0), createparents=True)
			
			# copy Sweep_Array to sweep_data_table
			sweep_data_table.append(self.Sweep_Array)

			# Save metadata
			for data in self.metadata.__dict__.keys():
				exec('sweep_data_table.attrs.{0} = self.metadata.{0}'.format(data))
				if self.metadata.__dict__[data] == None:
					print('table metadata {0} not defined and is set to None'.format(data))	

			sweep_data_table.flush()	

			# try:
			# 	TOC = fileh.get_node('/Contents') # is a table
			# except:
			# 	print('Creating h5 data set table of contents')
			# 	TOC = fileh.create_table('/', 'Contents', self.data_set_contents, "Table listing all tables contained in h5 file", tables.Filters(0)) #tables.Filters(0) means there is no data compression

			# TOC.append()

		# title = 'Data from Run ' + self.metadata.Run + ', Sensor: ' + self.metadata.Sensor + ', Ground Plane: ' + self.metadata.Ground_Plane

		# #determine type of  measurement...
		# if  (self.Sweep_Array.size == 1) | (np.abs(self.Sweep_Array['Fstop'] - self.Sweep_Array['Fstart']).max() >= 100e6):
		# 	groupname = 'Survey'
		# elif (np.unique(swp.Sweep_Array['Heater_Voltage']).size > 1) && (np.unique(swp.Sweep_Array['Pinput_dB']).size == 1):
		# 	groupname = 'T_Sweep'
		# elif (np.unique(swp.Sweep_Array['Heater_Voltage']).size == 1) && (np.unique(swp.Sweep_Array['Pinput_dB']).size > 1):
		# 	groupname = 'P_Sweep'
		# elif (np.unique(swp.Sweep_Array['Heater_Voltage']).size > 1) && (np.unique(swp.Sweep_Array['Pinput_dB']).size > 1):
		# 	groupname = 'TP_Sweep'
		# else:
		# 	groupname = 'Sweep'

		# 	groupname = 'T' + str(np.unique(swp.Sweep_Array['Heater_Voltage']).size) + 'P' +  str(np.unique(swp.Sweep_Array['Pinput_dB']).size)	

	def decompress_gain(self, Compression_Calibration_Index = -1, Show_Plot = True, Verbose = True):
		''' Assumes the two lowest input powers of the power sweep are not gain compressed, thus
		cannot be used if the two lowest powers are gain compressed. '''
		from scipy.interpolate import interp1d
		from scipy.optimize import curve_fit
		from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator

		Sweep_Array_Record_Index = self.loop.index 
		V = self.Sweep_Array['Heater_Voltage'][Sweep_Array_Record_Index]
		Fs = self.Sweep_Array['Fstart'][Sweep_Array_Record_Index]
		P = self.Sweep_Array['Pinput_dB'][Sweep_Array_Record_Index]

		Sweep_Array = np.extract((self.Sweep_Array['Heater_Voltage'] == V) & ( self.Sweep_Array['Fstart']==Fs) , self.Sweep_Array)


		num_sweep_powers = Sweep_Array['Pinput_dB'].shape[0]

		if num_sweep_powers <= 4:
			print('Number of sweep powers, {0}, is insufficient to perform gain decompression.'.format(num_sweep_powers))
			return
		#else:
		#	print('Performing gain decompression on {0} sweep powers.'.format(num_sweep_powers))

		Pin = np.power(10, Sweep_Array['Pinput_dB']/10.0) #mW, Probe Power

		#ChooseCompression calobration data from Power Sweep Data. 
		#It is the S21(Compression_Calibration_Index) for every sweep power 
		compression_calibration_data = np.power(np.abs(Sweep_Array['S21'][:,Compression_Calibration_Index]),2) #Pout/Pin,  

		Pout = compression_calibration_data*Pin

		#calculated_power_gain is power gain calculated from the slope of the two smallest input powers in Pin
		min_index = np.where(Pin == Pin.min())[0][0] # Index of the min values of Pin, unpacked from tuple
		dif = Pin - Pin.min() 
		min_plus = dif[np.nonzero(dif)].min() + Pin.min() # Second lowest value of Pin
		min_plus_index = np.where(np.isclose(Pin,min_plus))[0][0] # index of the second lowest Pin value, unpacked from tuple
		# When Pin = 0, 0 != Pout = Pin*gaain. There is an offset, i.e. a y-intercept, b, such at y = m*x+b. Next, we find m.  
		calculated_power_gain = (Pout[min_plus_index] - Pout[min_index])/(Pin[min_plus_index ]-Pin[min_index]) 

		#Pout_ideal is the output power assuming linear gain
		Pout_ideal = lambda p_in: calculated_power_gain*(p_in-Pin[0]) + Pout[0]

		Probe_Power_Mag = np.power(10,self.Sweep_Array[Sweep_Array_Record_Index]['Pinput_dB']/10)
		S21 = self.Sweep_Array[Sweep_Array_Record_Index]['S21']
		S21_Pout = np.power(np.abs(S21),2)*Probe_Power_Mag

		# create funcation to what Pin would be at an arbitrary Pout
		decompression_function = interp1d(Pout,Pin,kind = 'linear')

		# for polynomial to Pout vs Pin curve and use this to extrapolate values where Pout in not in interpolation domain
		def decompression_function_fit(pout, a,b,c):
			return a*np.power(pout,2)+b*pout+c
		popt,pcov = curve_fit(decompression_function_fit, Pout, Pin)
		decompression_function_extrap = lambda pout : decompression_function_fit(pout,popt[0],popt[1],popt[2])

		
		def decompress_element(z):
			z_Pout = np.power(np.abs(z),2)*Probe_Power_Mag
			if z_Pout <= Pout.min(): #Do nothinge when z_Pout is less than the interpolation range, Pout.min() to Pout.max()
				return z
			elif Pout.min() < z_Pout < Pout.max(): # Interpolate to find ideal Pout (assuming linear gain) when z_Pout is in interpolation domain 
				return z*np.sqrt(Pout_ideal(decompression_function(z_Pout))/Probe_Power_Mag)/np.abs(z)
			else: # Pout.max() <= z_Pout --  Extrapolate to find ideal Pout when z_Pout is above interpolation domain
				return z*np.sqrt(Pout_ideal(decompression_function_extrap(z_Pout))/Probe_Power_Mag)/np.abs(z)

		decompress_array = np.vectorize(decompress_element) # Vectorize for speed

		self.loop.z = S21_Decompressed = decompress_array(S21)

		if Verbose == True:
			print('Gain decompression calculation is based on {0} sweep powers.'.format(num_sweep_powers))
			print('Power out at zero input power is {0} mW'.format(calculated_power_gain*(0-Pin[0]) + Pout[0]))

		if Show_Plot:
			fig1           = plt.figure(figsize = (15,5))
			Freq           = self.Sweep_Array[Sweep_Array_Record_Index]['Frequencies']
			#majorFormatter = FormatStrFormatter('%d')
			majormaxnlocator    = MaxNLocator(nbins = 5)
			minormaxnlocator    = MaxNLocator(nbins = 5*5)
			#minorLocator   = MultipleLocator((Freq.max() - Freq.min())/25)
			

			ax1 = fig1.add_subplot(131)
			ax1.set_xlabel('Power In [mW]')
			line1 = ax1.plot(Pin,Pout, 'b-', label = 'Measured')
			line2 = ax1.plot(Pin,Pout_ideal(Pin), 'r-', label = 'Ideal')
			ax1.set_ylabel('Power Out [mW]', color='b')
			ax1.set_title('Gain Compression', fontsize=9)
			ax1.legend(loc = 'best', fontsize=9)
			plt.setp(ax1.get_xticklabels(),rotation = 45, fontsize=9)
			ax1.grid()
			#fig1.canvas.manager.resize(800,800)

			
			ax2 = fig1.add_subplot(132, aspect='equal')
			line2 = ax2.plot(S21.real,S21.imag, color='blue', linestyle='solid', linewidth = 3, label = 'Measured') 
			line1 = ax2.plot(S21_Decompressed.real, S21_Decompressed.imag, 'g-',linewidth = 3, label = 'Corrected')
			ax2.grid()
			ax2.set_title('Resonance Loop', fontsize=9)
			plt.setp(ax2.get_xticklabels(),rotation = 45)
			#ax2.legend(loc = 'best')

			
			ax3 = fig1.add_subplot(133)
			ax3.set_xlabel('Freq [Hz]')
			line1 = ax3.plot(Freq,10*np.log10(np.abs(S21)), 'b-',label = 'Measured',linewidth = 3)
			line2 = ax3.plot(Freq,10*np.log10(np.abs(S21_Decompressed)), 'g-', label = 'Corrected',linewidth = 3)
			ax3.set_ylabel('$|S_{21}|$ [dB]', color='k')
			ax3.legend(loc = 'best', fontsize=9)
			ax3.xaxis.set_major_locator(majormaxnlocator)
			#ax3.tick_params( axis='both', labelsize=9)
			plt.setp(ax3.get_xticklabels(),rotation = 45, fontsize=9)
			#ax3.xaxis.set_major_formatter(majorFormatter)
			ax3.xaxis.set_minor_locator(minormaxnlocator)
			ax3.set_title('Resonance Dip', fontsize=9)
			ax3.grid()

			fig1.subplots_adjust(wspace = 0.6,bottom = 0.09, top = 0.1)
			fig1.suptitle('Run: {0}, Sensor: {1}, Ground Plane: {2}, Readout Power: {3} dBm, Date: {4}'.format(self.metadata.Run, self.metadata.Sensor,self.metadata.Ground_Plane,self.Sweep_Array[Sweep_Array_Record_Index]['Pinput_dB'],self.metadata.Time_Created), fontsize=10)
			#plt.tight_layout()
			plt.setp(fig1, tight_layout = True)
			plt.show()

	def sweep_array_info(self):
		''' prints information about the Sweep_Array currently loaded'''
		Input_Powers = np.unique(self.Sweep_Array['Pinput_dB'])
		Heater_Voltages = np.unique(self.Sweep_Array['Heater_Voltage'])
		Number_of_Freq_Ranges = max(np.unique(self.Sweep_Array['Fstart']),np.unique(self.Sweep_Array['Fstop']))
		print('{0:03.0f} - Total number of sweeps.\n{1:03.0f} - Number of readout powers.\n{2:03.0f} - Number of readout temperatures.\n{3:03.0f} - Number of frequency bands.'.format(
			self.Sweep_Array.shape[0],
			Input_Powers.shape[0],
			Heater_Voltages.shape[0],
			Number_of_Freq_Ranges.shape[0]))

	def construct_hf5_toc(self,filename = database_location):
		''' Creates a table of contents (toc) of the hf5 database storing all the sweep_data.
		very useful for finding the name and locaiton of a table in the database'''
		if not os.path.isfile(filename):
			print('Speficied h5 database does not exist. Aborting...')
			return 
		
		wmode = 'a'

		# use "with" context manage to ensure file is always closed. no need for fileh.close()
		with tables.open_file(filename, mode = wmode) as fileh:
			table_list = [g for g in fileh.walk_nodes(classname = 'Table')]
			num_tables = len(table_list)
			TOC = np.zeros(num_tables, dtype = self.data_set_contents)
			index = 0
			for table in table_list:
				TOC['Run'][index] 				= table.get_attr('Run') 
				TOC['Time_Created'][index] 		= table.get_attr('Time_Created')
				#TOC['Num_Ranges'][index] 		= table.get_attr('Num_Ranges') if 'Num_Ranges' in table.attrs._v_attrnames else 1
				TOC['Num_Ranges'][index] 		= table.get_attr('Num_Ranges') if table.get_attr('Num_Ranges') !=None else 0
				TOC['Num_Powers'][index] 		= table.get_attr('Num_Powers') if table.get_attr('Num_Powers') !=None else 0
				TOC['Num_Temperatures'][index] 	= table.get_attr('Num_Heater_Voltages') if table.get_attr('Num_Heater_Voltages') !=None else 0
				TOC['Sensor'][index] 			= table.get_attr('Sensor') if table.get_attr('Sensor') !=None else ''
				TOC['Ground_Plane'][index] 		= table.get_attr('Ground_Plane') if table.get_attr('Ground_Plane') !=None  else ''
				TOC['Path'][index] 				= table._v_pathname
				index += 1 

		self.TOC = TOC
		print(TOC)

	def load_hf5(self, tablepath, filename = database_location):
		''' table path is path to the database to be loaded starting from root. e.g. swp.load_hf5('/Run44b/T201312102229')
		filename is the name of the hf5 database to be accessed for the  table informaiton'''

		if not os.path.isfile(filename):
			print('Speficied h5 database does not exist. Aborting...')
			return 
		
		wmode = 'a'
		
		#delete previous metadata object
		del(self.metadata)
		self.metadata = metadata()
		del(self.loop)
		self.loop = loop()

		# use "with" context manage to ensure file is always closed. no need for fileh.close()
		with tables.open_file(filename, mode = wmode) as fileh:
			table = fileh.get_node(tablepath)	
			self.Sweep_Array = table.read()
			for data in self.metadata.__dict__.keys():
				try:
					exec('self.metadata.{0} = table.attrs.{0}'.format(data))
				except:
					print('Table metadata is missing {0}. Setting to None'.format(data))
					exec('self.metadata.{0} = None'.format(data))
		self.sweep_data_columns = self.Sweep_Array.dtype

	def pick_loop(self,index):
		'''Use this function to pick the current loop/transmission data from withing the Sweep_Array. 
		Index is the indes number of sweep/loop to be slected as the current loop.'''
		self.loop.index = index 
		self.loop.z = self.Sweep_Array[index]['S21']
		self.loop.freq = self.Sweep_Array[index]['Frequencies']
		
	def remove_cable_delay(self, Show_Plot = True, Verbose = True):
		

		S21 = self.loop.z
		Freq = self.loop.freq

		j = np.complex(0,1)
		

		# 'base' and 'offset' determine the indices of S21 to use for minimization.
		# Basically, we minimize the distance between complex vectors points S21[base] and S21[base + offset]
		# rotation angle and rescaling factor computed from 'base'
		base = 0
		offset = -1
		normalization = np.abs(S21[base])

		def obj(t):
			''' Build objective function to minimize
			this function outputs the complex vector difference
			cable_delay_term(base)S21(base) - cable_delay_term(base+offset)S21(base+offset)
			'''
			S21a = np.exp(2*np.pi*Freq[base]*j*t)*S21[base]
			S21b = np.exp(2*np.pi*Freq[base+offset]*j*t)*S21[base+offset]
			return np.abs(S21a-S21b)

		cable_delay_guess = 100*10**-9 # Seconds
		out = minimize(obj,cable_delay_guess, method='Nelder-Mead')
		
		cable_delay = out.x[0] # in seconds
		#cable_delay = 80e-9


		# rotation angel to align circle with x-axis
		angle = 0 # angle = 0 -- do not aligh with x-axis

		rescale_factor = np.exp(j*-1*angle)/normalization 
		#For Later: Need to make sure the noise in S21 (perhaps standard dev) is  much much less than distance between S21[base] and S21[base + offset]

		S21_Corrected = rescale_factor*np.exp(2*np.pi*Freq*j*cable_delay)*S21

		if Verbose == True:
			print('cable delay is {} seconds'.format(cable_delay ))

		if Show_Plot:
			from matplotlib.ticker import MaxNLocator
			fig1           = plt.figure(figsize = (6,6))
			majormaxnlocator    = MaxNLocator(nbins = 5)
			minormaxnlocator    = MaxNLocator(nbins = 5*5)
			ax2 = fig1.add_subplot(111, aspect='equal')
			line2 = ax2.plot(S21.real,S21.imag, color='blue', linestyle='solid', linewidth = 3, label = 'Measured') 

			
			line1 = ax2.plot(S21_Corrected.real, S21_Corrected.imag, 'g-',linewidth = 3, label = 'Corrected')
			ax2.grid()
			ax2.set_title('Resonance Loop', fontsize=9)
			plt.setp(ax2.get_xticklabels(),rotation = 45)
			ax2.legend(loc = 'best')

			fig1.subplots_adjust(wspace = 0.6,bottom = 0.09, top = 0.1)
			#plt.tight_layout()
			plt.setp(fig1, tight_layout = True)
			plt.show()

		self.loop.z = S21_Corrected

	def trim_loop(self,N = 20,Verbose = True,):
		import numpy.ma as ma
		f = f1 = ma.array(self.loop.freq)
		z = z1 = ma.array(self.loop.z)
		# estimate resonant freq using resonance dip
		zr_mag_est = np.abs(z).min()
		zr_est_index = np.where(np.abs(z)==zr_mag_est)[0][0]

		# estimate max transmission mag using max valuse of abs(z)
		z_max_mag = np.abs(z).max()

		#Depth of resonance in dB
		depth_est = 20.0*np.log10(zr_mag_est/z_max_mag)

		#Magnitude of resonance dip at half max
		res_half_max_mag = (z_max_mag+zr_mag_est)/2

		#find the indices of the closest points to this magnitude along the loop, one below zr_mag_est and one above zr_mag_est
		a = np.square(np.abs(z[:zr_est_index+1]) - res_half_max_mag)
		lower_index = np.argmin(a)
		a = np.square(np.abs(z[zr_est_index:]) - res_half_max_mag)
		upper_index = np.argmin(a) + zr_est_index

		#estimate the FWHM bandwidth of the resonance
		f_upper_FWHM = f[upper_index]
		f_lower_FWHM = f[lower_index]
		FWHM_est = np.abs(f_upper_FWHM - f_lower_FWHM)
		fr_est = f[zr_est_index]

		#Bandwidth Cut: cut data that is more than N * FWHM_est away from zr_mag_est
		z = z2 = ma.masked_where((f > fr_est + N*FWHM_est) | (fr_est - N*FWHM_est > f),z)
		f = f2 = ma.array(f,mask = z.mask)

		self.loop.z = ma.compressed(z)
		self.loop.freq = ma.compressed(f)

		if Verbose: 
			print('Bandwidth cut:\n\t{1} points outside of fr_est +/- {0}*FWHM_est removed, {2} remaining data points'.format(N, *self._points_removed(z1,z2)))

	def _points_removed(self,initial, final):
		''' Compute and return the number of point removed from inital due to a cut. 
		return this number and the number of points in final'''
		try:
			initial_number = initial.size - initial.mask.sum()
		except:
			initial_number = initial.size

		try: 	
			final_number = final.size - final.mask.sum()
		except:
			final_number = final.size
		return (initial_number - final_number), final_number

	def circle_fit(self, Show_Plot = True):

		S21 = self.loop.z
		Freq = self.loop.freq

		LargeCircle = 10
		def pythag(m,n):
			'''compute pythagorean distance
			   sqrt(m*m + n*n)'''
			return np.sqrt(np.square(m) + np.square(n))

		def eigen2x2(a,b,c):
			'''a,b,c - matrix components 	[[a c]
											 [c d]]	
			   d1,d2 - eigen values where |d1| >= |d2|
			   (Vx,Vy) - unit eigen vector of d1,  Note: (-Vy,Vx) is eigen vector for d2
			'''
			disc = pythag(a-b,2*c) # discriminant
			d1 = max(a+b + disc, a+b - disc)/2
			d2 = (a*b-c*c)/d1

			if np.abs(a-d1) > np.abs(b-d1):
				f = pythag(c,d1-a)
				if f == 0.0:
					Vx = 1.
					Vy = 0.
				else:
					Vx = c/f
					Vy = (d1-a)/f
			else:
				f = pythag(c,d1-b)
				if f == 0.0:
					Vx = 1.
					Vy = 0.
				else:
					Vx = (d1-b)/f
					Vy = c/f					
			return d1,d2,Vx,Vy

		def F(x,y,a,b):
			''' computes and returns the value of the objective fuction.
			do this for the case of a large circle and a small circle  '''
			
			if (np.abs(a) < LargeCircle) and (np.abs(b) < LargeCircle): # Case of Small circle
				xx = x - a
				yy = y - b
				D = pythag(xx,yy)

				r = D.mean()

				return (np.square(D - r)).mean()
			else:	# Case of Large circle
				a0 = a - x.mean()
				b0 = b - y.mean()
				d = 1.0/pythag(a0,b0)
				dd = d*d
				s = b0*d
				c = a0*d

				xx = x - x.mean()
				yy = y - y.mean()
				z = np.square(xx) + np.square(yy)
				p = xx*c + yy*s
				t = d*z - 2.0*p
				g = t/(1+np.sqrt(1+d*t))
				W = (z+p*g)/(2.0+d*g)
				Z = z

				return Z.mean() - W.mean()*(2.0+d*d*W.mean())

		def GradHessF(x,y,a,b):
			'''Compute gradient of F, GradF = [F1,F2] and Hessian of F, HessF = [[A11 A12]
																				  A12 A22]]
			at point p = [a,b].
			Note Hessian is symmetric. 
			'''
			if (np.abs(a) < LargeCircle) and (np.abs(b) < LargeCircle): # Case of Small circle
				xx = x - a
				yy = y - b
				r = pythag(xx,yy)
				u = xx/r
				v = yy/r

				Mr = r.mean()
				Mu = u.mean()
				Mv = v.mean()
				Muu = (u*u).mean()
				Mvv = (v*v).mean()
				Muv = (u*v).mean()
				Muur = (u*u/r).mean()
				Mvvr = (v*v/r).mean()
				Muvr = (u*v/r).mean()
			


				F1 = a + Mu * Mr - x.mean()
				F2 = b + Mv * Mr - y.mean()

				A11 = 1.0 - Mu * Mu - Mr * Mvvr
				A22 = 1.0 - Mv * Mv - Mr * Muur
				A12 = -1.0 * Mu * Mv + Mr * Muvr

			else:	# Case of Large circle
				a0 = a - x.mean()
				b0 = b - y.mean()
				d = 1.0/pythag(a0,b0)
				dd = d*d
				s = b0*d
				c = a0*d

				xx = x - x.mean()
				yy = y - y.mean()
				z = np.square(xx) + np.square(yy)
				p = xx*c + yy*s
				t = 2.0*p - d*z 
				w = np.sqrt(1.0-d*t)				
				g = -1*t/(1.0+w)
				g1 = 1.0/(1.0+d*g)
				gg1 = g*g1
				gg2 = g/(2.0 + d * g)
				aa = (xx+g*c)/w
				bb = (yy+g*s)/w	

				X = (xx*gg1).mean()
				Y = (yy*gg1).mean()
				R = (z+t*gg2).mean()
				T = (t*gg1).mean()
				W = (t*gg1*gg2).mean()	
				AA = (aa*aa*g1).mean()
				BB = (bb*bb*g1).mean()
				AB = (aa*bb*g1).mean()
				AG = (aa*gg1).mean()
				BG = (bb*gg1).mean()
				GG = (g*gg1).mean()	

				U = (T-b*W)*c*0.5 - X + R*c*0.5
				V = (T-b*W)*s*0.5 - Y + R*s*0.5

				F1 = d * ((dd*R*U - d*W*c + T*c)*0.5 - X)
				F2 = d * ((dd*R*V - d*W*s + T*s)*0.5 - Y)

				UUR = ((GG-R*0.5)*c + 2.0*(AG-U))*c + AA
				VVR = ((GG-R*0.5)*s + 2.0*(BG-V))*s + BB
				UVR = ((GG-R*0.5)*c + (AG-U))*s + (BG-V)*c + AB

				A11 = dd*(U*(2.0*c-dd*U) - R*s*s*0.5 - VVR*(1.0+dd*R*0.5))
				A22 = dd*(V*(2.0*s-dd*V) - R*c*c*0.5 - UUR*(1.0+dd*R*0.5))
				A12 = dd*(U*s + V*c + R*s*c*0.5 - dd*U*V + UVR*(1.0 + dd*R*0.5))
			return F1,F2,A11,A22,A12
		
		def sigma(x,y,loop):
			'''estimate of Sigma = square root of RSS divided by N
			gives the root-mean-square error of the geometric circle fit'''
			dx = x-loop.a
			dy = x-loop.b
			loop.sigma = (pythag(dx,dy)-loop.r).mean()
			return loop

		def CircleFitByChernovHoussam(x,y, init, lambda_init):
			import copy
			import sys


			REAL_EPSILON = sys.float_info.epsilon
			REAL_MAX = sys.float_info.max


			IterMAX=200
			check_line= True
			#dmin = 1.0

			ParLimit2 = 100.
			epsilon = 1.e+7*REAL_EPSILON
			factor1 = 32.
			factor2 = 32.
			ccc = 0.4
			factorUp = 10.
			factorDown = 0.1

			new = copy.copy(init)
			#new = sys.modules[__name__].loop() #This is how to access the loop class from inside this function
			#old = loop()

			new.s = F(x,y,init.a,init.b) # compute root mean square error
			F1,F2,A11,A22,A12 = GradHessF(x,y,init.a,init.b) # compute gradient vector and Hessian matrix
			new.Gx = F1
			new.Gy = F2
			new.g = pythag(F1,F2) # The gradient vector and its norm
			lambda_ = lambda_init
			sBest = gBest = progess = REAL_MAX

			enough = False
			i = 0
			ii = 0
			while not enough:
				if i > 0:
					# evaluate the progress made during the previous iteration
					progress = (np.abs(new.a - old.a)+np.abs(new.b - old.b))/(np.square(old.a) + np.square(old.b) + 1.0)
				old = copy.copy(new)

				i = i+1
				if i > IterMAX: #termination due to going over the limit
					enough = True
					break
				d1,d2,Vx,Vy = eigen2x2(A11,A22,A12) #eigendecomposition of the Hessian matrix
				dmin = min(d1,d2) #recording the smaller e-value
				AB = pythag(old.a,old.b) + 1.0 # approximation to the circle size
				# main stopping rule: terminate iterations if 
				# the gradient vector is small enough and the 
				# progress is not substantial 
				if (old.g < factor1*REAL_EPSILON) and (progress<epsilon):
					#print('primary stopping rule')
					enough = True
					break
				# secondary stopping rule (prevents some stupid cycling)
				if (old.s >= sBest) and (old.g >= gBest):
					print(old.s, sBest, old.g, gBest)
					#print('secondary stopping rule')
					enough = True
					break

				if (sBest > old.s):
					sBest = old.s  # updating the smallest value of the objective function found so far
				if (gBest > old.g): 
					gBest = old.g  # updating the smallest length of the gradient vector found so far

				G1 = Vx*F1 + Vy*F2  # rotating the gradient vector
				G2 = Vx*F2 - Vy*F1  # (expressing it in the eigensystem of the Hessian matrix)

				while not enough: # starting point of an "inner" iteration (adjusting lambda)
					# enforcing a lower bound on lambda that guarantees that
					# (i)  the augmented Hessian matrix is positive definite
					# (ii) the step is not too big (does not exceed a certain 
					# fraction of the circle size) the fraction is defined by 
					# the factor "ccc")
					if lambda_ < (np.abs(G1)/AB/ccc) - d1:
						lambda_ = np.abs(G1)/AB/ccc - d1
					if lambda_ < (np.abs(G2)/AB/ccc) - d2: 
						lambda_ = np.abs(G2)/AB/ccc - d2

					# compute the step (dX,dY) by using the current value of lambda
					dX = old.Gx*(Vx*Vx/(d1+lambda_)+Vy*Vy/(d2+lambda_)) + old.Gy*Vx*Vy*(1.0/(d1+lambda_)-1.0/(d2+lambda_))
					dY = old.Gx*Vx*Vy*(1.0/(d1+lambda_)-1.0/(d2+lambda_)) + old.Gy*(Vx*Vx/(d2+lambda_)+Vy*Vy/(d1+lambda_))

					# updating the loop parameter
					new.a = old.a - dX
					new.b = old.b - dY

					if (new.a==old.a) and (new.b==old.b): #if no change, terminate iterations
						enough  = True
						break

					#check if the circle is very large
					if np.abs(new.a)>ParLimit2 or np.abs(new.b)>ParLimit2:
						#when the circle is very large for the first time, check if 
						#the best fitting line gives the best fit
						if check_line:   # initially, check_line= True, then it is set to zero

							#compute scatter matrix
							dx = x - x.mean()
							dy = y - y.mean()
							Mxx = (dx*dx).sum()
							Myy = (dy*dy).sum()
							Mxy = (dy*dx).sum()
							dL1,dL2,VLx,VLy = eigen2x2(Mxx,Myy,Mxy)  # eigendecomposition of scatter matrix

							#compute the third mixed moment (after rotation of coordinates)
							dx = (x - x.mean())*VLx + (y - y.mean())*VLy
							dy = (y - y.mean())*VLx - (x - x.mean())*VLy
							Mxxy = (dx*dx*dy).sum()

							#rough estimate of the center to be used later to recoved from the wrong valley
							if Mxxy > 0.0:
								R = ParLimit2
							else:
								R = -ParLimit2

							aL = -VLy*R
							bL =  VLx*R                 
							check_line = False

						# check if the circle is in the wrong valley
						if (new.a*VLy - new.b*VLx)*R>0.0: 
							# switch to the rough circle estimate (precomupted earlier)
							new.a = aL;                 
							new.b = bL;                 
							new.s = F(x,y,new.a,new.b)    # compute the root-mean-square error
							
							# compute the gradient vector and Hessian matrix
							F1,F2,A11,A22,A12 = GradHessF(x,y,new.a,new.b)  

							# the gradient vector and its norm 
							new.Gx = F1;  
							new.Gy = F2;   
							new.g = pythag(F1,F2)  
							lambda_ = LambdaIni        #reset lambda
							sBest = gBest = REAL_MAX  #reset best circle characteristics 
							break
					
					# compute the root-mean-square error
					new.s = F(x,y,new.a,new.b) 
					# compute the gradient vector and Hessian matrix
					F1,F2,A11,A22,A12 = GradHessF(x,y,new.a,new.b)

					# the gradient vector and its norm  
					new.Gx = F1  
					new.Gy = F2   
					new.g = pythag(F1,F2) 

					# check if improvement is gained
					if new.s < sBest*(1.0+factor2*REAL_EPSILON):  #yes, improvement
						lambda_ *= factorDown     # reduce lambda
						break 
					else:
						ii += 1
						if ii > IterMAX: #termination due to going over the limit
							enough = True
							break
						lambda_ *= factorUp #increace lambda
						continue
			

			old.r = pythag(x - old.a, y - old.b).mean() 
			old.outer_iterations = i
			old.inner_iterations = ii
			loop = old
			exit_code = 0
			if old.outer_iterations  > IterMAX:
				exit_code  = 1

			if old.inner_iterations  > IterMAX:
				exit_code = 3

			if (dmin <= 0.0) and (exit_code==0):
				exit_code  = 3

			loop.circle_fit_exit_code = exit_code
			loop = sigma(x,y,loop)

			return loop
		


		#initial guess
		self.loop.a =  0
		self.loop.b =  0
		lambda_init = 0.001

		x = S21.real
		y = S21.imag
	
		self.loop = CircleFitByChernovHoussam(x,y, self.loop, lambda_init)

		if Show_Plot:
			fig, ax = self.plot_loop(show = False)[:2]		
			t = np.linspace(0, 2.0*np.pi, num=50, endpoint=True)
			j = np.complex(0,1); zc = self.loop.a + j*self.loop.b;  r = self.loop.r
			line = ax.plot(zc.real + r*np.cos(t),zc.imag + r*np.sin(t),'y-', label = 'Circle Fit')
			line = ax.plot([zc.real],[zc.imag],'yx', markersize = 10, markeredgewidth = 4, label = 'Center')

			plt.show()	

	def phase_fit(self, Fit_Method = 'Multiple', Verbose = True, Show_Plot = True):
		from scipy.stats import chisquare
		def angle(z, deg = 0):
			''' If z is a masked array. angle(z) returns the angle of the elements of z
			within the branch [0,360] instead of [-180, 180], which is the branch used
			in np.angle(). The mask of angle(z) is set to be the mask of the input, z.

			If z is not a masked array, then angle(z) is the same as np.angle except 
			that range is [0,360] instead of [-180, 180]

			If z is a vector, then an angle shift is added to z  so the z[0] is 0 degrees
			If z is a number, then dont shift angle'''
			a = np.angle(z, deg = deg)
			try:
				a = a - a[0] #if a is not a vector, then a[0] will throw an error
			except:
				pass
			p = np.where(a<=0,1,0)
			n = 2
			units = n*np.pi if deg == 0 else n*180
			try:
				a = ma.array(a + p*units,mask =z.mask) 
			except:
				a = a + p*units #if z is not a masked array do this
			return a

		j = np.complex(0,1)
		try:
			zc = self.loop.a + j*self.loop.b
			r = self.loop.r
		except:
			print('Phase fit needs loop center and radius, which are not currently defined. Aborting phase fit.')
			return
		f = f0 = self.loop.freq
		z = z0 = self.loop.z
		

		# Remove duplicate frequency elements in z and f, e.g. places where f[n] = f[n+1]
		f_adjacent_distance = np.abs(f[:-1]-f[1:])
		z = z1 = ma.masked_where(f_adjacent_distance==0.0, z[:-1])
		f = f1 = ma.array(f[:-1],mask = z.mask) #Syncronize mask of f to match mask of z


		#Estimate Resonance frequency using minimum Dip or max adjacent distance
		Use_Dip = 1 
		if Use_Dip: #better for non-linear resonances with point near loop center
			zr_mag_est = np.abs(z).min()
			zr_est_index = np.where(np.abs(z)==zr_mag_est)[0][0]
		else:
			z_adjacent_distance = np.abs(z[:-1]-z[1:])
			zr_est_index = np.argmax(z_adjacent_distance) 
			zr_mag_est = np.abs(z[zr_est_index])


		#Transmission magnitude off resonance 
		Use_Fit = 1
		if Use_Fit:
			z_max_mag = np.abs(zc)+r
		else: #suspected to be better for non-linear resonances
			z_max_mag = np.abs(z).max()

		#Depth of resonance in dB
		depth_est = 20.0*np.log10(zr_mag_est/z_max_mag)

		#Magnitude of resonance dip at half max
		res_half_max_mag = (z_max_mag+zr_mag_est)/2

		#find the indices of the closest points to this magnitude along the loop, one below zr_mag_est and one above zr_mag_est
		a = np.square(np.abs(z[:zr_est_index+1]) - res_half_max_mag)
		lower_index = np.argmin(a)#np.where(a == a.min())[0][0]
		a = np.square(np.abs(z[zr_est_index:]) - res_half_max_mag)
		upper_index = np.argmin(a) + zr_est_index

		#estimate the FWHM bandwidth of the resonance
		f_upper_FWHM = f[upper_index]
		f_lower_FWHM = f[lower_index]
		FWHM_est = np.abs(f_upper_FWHM - f_lower_FWHM)
		fr_est = f[zr_est_index]
		

		
		#translate circle to origin, and rotate so that z[zr_est_index] has angle 0
		z = z2 = ma.array((z.data-zc)*np.exp(-j*(angle(zc)-np.pi)), mask = z.mask)
		
		
		
		#remove points that occur within r_cutoff of the origin of the centered data. 
		#(For non-linear resonances that have spurious point close to loop center)	
		r_fraction = 0.75
		r_cutoff  = r_fraction*r
		z = z3 = ma.masked_where(np.abs(z)<r_cutoff,z)
		f = f3 = ma.array(f,mask = z.mask)
		

		#Bandwidth Cut: cut data that is more than N * FWHM_est away from zr_mag_est
		N = 10
		z = z4 = ma.masked_where((f > fr_est + N*FWHM_est) | (fr_est - N*FWHM_est > f),z)
		f = f4 = ma.array(f,mask = z.mask)
		z_theta = angle(z)


		#Angle jump cut : masks points where angle jumps to next branch, +/- theta_cutoff
		theta_cutoff = 345 #degrees
		mask = (f > fr_est + 0.5*FWHM_est) | (f < fr_est + -0.5*FWHM_est)
		f_in_FWHM = ma.masked_where(mask,f) # or alternatively: f_in_FWHM = f; f_in_FWHM[mask] = ma.masked 
		edge1,edge2 = ma.flatnotmasked_edges(f_in_FWHM)
		angle_slope = (z_theta[edge2]-z_theta[edge1])/(f[edge2]-f[edge1]) # angle is decreasing if negative slope
		upper_cond = ((f > fr_est +  0.5*FWHM_est) & ((z_theta[edge2]<z_theta) if (angle_slope<0) else (z_theta[edge2]>z_theta))) 
		lower_cond = ((f < fr_est + -0.5*FWHM_est) & ((z_theta[edge1]>z_theta) if (angle_slope<0) else (z_theta[edge1]<z_theta))) 
		z = z5 = ma.masked_where(lower_cond|upper_cond,z)
		f = f5 = ma.array(f,mask = z.mask)
		z_theta = z_theta5 = ma.array(z_theta,mask = z.mask)
		



		theta_est = np.extract(f==fr_est,z_theta)[0]
		Q_est = fr_est/FWHM_est


		#consider reducing computation by extracting only the unmasked values of z,f, and z_theta of the minimization
		#These commands return a masked array where all the masked elements are removed.
		#z = z[~z.mask]
		#f = f[~f.mask]
		#z_theta = z_theta[~z_theta.mask]

		# These commands return np array
		#z = ma.compressed(z)
		#f = ma.compressed(f)
		#z_theta  = ma.compressed(z_theta)
		

		if mysys.startswith('Windows'):
			dt = np.float64
		else:	
			dt = np.float128

		def hess(x, z_theta,f): #to avoid overflow try to re write hessian so that all numbers are of order 1
			theta,fr,Q = x	
			H = np.zeros((3,3), dtype = dt)
			ff = (1-(f/fr))
			denom = (1+4.0*np.square(ff*Q))
			numer = (theta+z_theta-2.0*np.arctan(2.0*ff*Q))
			H[0,0] = (2.0*np.ones_like(z_theta)).sum()
			H[0,1] = ((-8.0*f*Q)/(np.square(fr)*denom)).sum()
			H[0,2] = ((8.0*ff)/denom).sum()
			H[1,0] = H[0,1] #((8.0*f*Q)/(np.square(fr)*denom)).sum()
			H[1,1] = ((32.0*np.square(f*Q/(np.square(fr)*denom)))  +   (64.0*np.square(f/(np.square(fr)*denom))*ff*np.power(Q,3)*numer)   +  ((16.0*f*Q/np.power(fr,3))*(numer/denom))).sum()
			H[1,2] = (((32.0*f*Q*ff)/np.square(fr*denom))  +  ((64.0*f*np.square(ff*Q)*numer)/(np.square(fr*denom)))  - ((8.0*f*numer)/(np.square(fr)*denom))).sum()
			H[2,0] = H[0,2] #((8.0*ff)/denom).sum()
			H[2,1] = H[1,2] #(((32.0*f*ff*Q)/np.square(fr*denom))  +  ((64.0*f*np.square(ff*Q)*numer)/(np.square(fr*denom)))  -  ((8.0*f*numer)/(np.square(fr)*denom))).sum()
			H[2,2] = (((32.0*np.square(ff))/np.square(denom))  +  ((64.0*np.power(ff,3)*Q*numer)/np.square(denom))).sum()				
			return H

		def jac(x,z_theta,f):
			theta,fr,Q = x
			J = np.zeros((3,),dtype = dt)    #np.zeros_like(x)
			ff = (1-(f/fr))
			denom = (1+4.0*np.square(ff*Q))
			numer = (theta+z_theta-2.0*np.arctan(2.0*ff*Q))	
			J[0] = np.sum(2.0*numer)
			J[1] = np.sum(-8.0*f*Q*numer/(np.square(fr)*denom))
			J[2] = np.sum(-8.0*ff*numer/denom)
			return J


		def obj(x,z_theta,f):
			theta,fr,Q = x
			return np.square(z_theta - theta - 2.0*np.arctan(2.0*Q*(1-f/fr))).sum()	 #<--- Need hessian of this


		def obj_ls(x,z_theta,f):
			'''object fuctinon for least squares fit'''
			theta,fr,Q = x
			residual  = z_theta - theta - 2.0*np.arctan(2.0*Q*(1-f/fr))	
			return residual

		p0 = np.array([theta_est,fr_est ,Q_est])
		
		fit_func = {}
		fit_func['Powell'] =  lambda : minimize(obj, p0, args=(z_theta,f), method='Powell', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-15, callback=None, options={'disp':False})
		fit_func['Nelder-Mead']  = lambda : minimize(obj, p0, args=(z_theta,f), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-15, callback=None, options={'disp':False, 'xtol' : 1e-6,'maxfev':1000})
		fit_func['Newton-CG'] = lambda : minimize(obj, p0, args=(z_theta,f), method='Newton-CG', jac=jac, hess=hess, hessp=None, bounds=None, constraints=(),tol=1e-15, callback=None, options={'maxiter' : 50,'xtol': 1e-4,'disp':False})

		fit = {}

		if Fit_Method == 'Multiple':
			for method in fit_func.keys():
				fit[method] = fit_func[method]()
		elif Fit_Method in fit_func.keys():
			fit[Fit_Method] = fit_func[Fit_Method]()
		else:
			print("Unrecognized fit method. Aborting fit. \n\t Must choose one of {0} or 'Multiple'".format(fit_func.keys()))
			return
		



		#Does not work if the objective function is re-arranged as in the following
		# print('Nelder-Mead 2 ################# ')
		# def obj(x,z_theta,f):
		# 	theta,fr,Q = x
		# 	return np.square(np.tan((z_theta - theta)/2) - (2.0*Q*(1-f/fr))).sum()
		# res = minimize(obj, p0, args=(z_theta,f), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-20, callback=None, options={'disp':True})
		# print(res)
	
		# Least square method does not find a good Q fit and the sum of the squares for solution is fairly high
		# print('Least Square ################# ')
		# print(fit['Least-Squares'])
		# print(np.square(fit['Least-Squares'][2]['fvec']).sum()) # this is the value of the sum of the squares for the solution
		# x = fit['Least-Squares'][0] 
		
		#x = res.x 
		bestfit = list(fit)[0]
		lowest = fit[bestfit].fun
		for key in fit.keys(): 
			if fit[key].fun < lowest:
				bestfit = key

		
		self.loop.Phase_Fit_Method = bestfit
		self.loop.Q = Q = fit[bestfit].x[2]
		self.loop.Qc = Qc = Q/(2*r)
		self.loop.Qi = Q*Qc/(Qc-Q)
		self.loop.fr = fr = fit[bestfit].x[1]
		self.loop.FWHM = fr/Q
		self.loop.phi = (fit[bestfit].x[0]-1*np.pi)*180/np.pi
		self.loop.chisquare, self.loop.pvalue = chisquare( z_theta,f_exp=fit[bestfit].x[0] + 2.0*np.arctan(2.0*Q*(1-f/fr)))
		
		#estimated quantities from MAG S21 
		self.loop.fr_est = fr_est
		self.loop.FWHM_est = FWHM_est
		self.loop.depth_est = depth_est
		self.loop.Q_est = Q_est

		if Verbose: 
			print('Duplicates cuts:\n\t{0} duplicate frequencies removed from loop data, {1} remaining data points'.format(*self._points_removed(z0,z1)))
			print('Radius cut:\n\t{1} points < r_loop*{0} found and removed, {2} remaining data points'.format(r_fraction,*self._points_removed(z2,z3)))
			print('Bandwidth cut:\n\t{1} points outside of fr_est +/- {0}*FWHM_est removed, {2} remaining data points'.format(N, *self._points_removed(z3,z4)))
			print('Angle jump cut:\n\t{1} points with loop angle step > {0} deg removed, {2} remaining data points'.format(theta_cutoff, *self._points_removed(z4,z5)))
			print('Initial Guess:\n\tLoop rotation {0}, fr {1}, Q {2}'.format(*p0))

			for method in fit.keys():
				print('\n{0} Minimzation Result:\n{1}\n'.format(method,fit[method]))



		if Show_Plot:
			total_removed, total_used_in_fit = self._points_removed(z1,z5)
			fig1 = plt.figure( facecolor = 'w',figsize = (10,10))
			ax = fig1.add_subplot(6,1,1)
			ax.set_title('Number of points used in fit = '+str(total_used_in_fit)+', Number of points removed = ' + str(total_removed) )
			line = ax.plot(f1[~f5.mask], np.abs(z1[~z5.mask]),'g-', label = 'Used for Fit')
			line = ax.plot(f1[f5.mask], np.abs(z1[z5.mask]),'r.',markersize = 2,  alpha = 0.2, label = 'Excluded Data')
			line = ax.plot([f1[zr_est_index],f1[zr_est_index]] , [np.abs(z1[zr_est_index]),np.abs(zc)+r] ,'k.', label = 'Magitude Min and Max')
			line = ax.plot([f1[lower_index], f1[upper_index], f1[upper_index]], np.abs([z1[lower_index],z1[lower_index],z1[upper_index]]),'yo-', label = 'FWHM Estimate')
			ax.set_ylabel('Magnitude')
			ax.legend(loc = 'best', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
			
			
			ax = fig1.add_subplot(6,1,(2,4), aspect='equal')
			t = np.linspace(0, 2.0*np.pi, num=50, endpoint=True)
			line = ax.plot([0,zc.real],[0, zc.imag],'y*-', label = 'Center Vector')	
			line = ax.plot(zc.real + r*np.cos(t),zc.imag + r*np.sin(t),'y-', label = 'Circle Fit')		
			line = ax.plot(z1.real, z1.imag,'r:', label = 'Initial Location')
			line = ax.plot(z3.real, z3.imag,'r-', label = 'Aligned w/ Origin')
			line = ax.plot(z4.real, z4.imag,'g:', linewidth = 3,label = 'Bandwidth Cut')
			pt = ax.plot([z1[0].real,z[~z.mask][0].real], [z1[0].imag,z[~z.mask][0].imag],'ko', label = 'First Point')
			pt = ax.plot(z2[zr_est_index].real, z2[zr_est_index].imag,'k*', label = 'Magnitude Min')
			line = ax.plot(z4[z4.mask].data.real, z4[z4.mask].data.imag,'r.', alpha = 0.2, label = 'Excluded Data')
			ax.legend(loc = 'best', fontsize=10, scatterpoints =1, numpoints = 1, labelspacing = .1)#,numpoints)
			
			text = ('$*Resonator Properties*$\n' + '$Q =$ ' + '{0:.2f}'.format(self.loop.Q) +'\nf$_0$ = ' + '{0:.6f}'.format(self.loop.fr/1e6) 
				+  ' MHz\n$Q_c$ = ' + '{0:.2f}'.format(self.loop.Qc) + '\n$Q_i$ = ' + '{0:.2f}'.format(self.loop.Qi) + '\n|S$_{21}$|$_{min}$ = ' 
				+ '{0:.2f}'.format(self.loop.depth_est) + ' dB' + '\nBW$_{FWHM}$ = ' + '{0:.3f}'.format(self.loop.FWHM/1e3) +  ' kHz' 
				+ '\n$\chi^{2}$ = ' + '{0:.4f}'.format(self.loop.chisquare) + '\n$\phi$ = ' + '{0:.3f}'.format(self.loop.phi) +' deg')
			bbox_args = dict(boxstyle="round", fc="0.8")        
			fig1.text(0.10,0.7,text,
					ha="center", va="top", visible = True,
					bbox=bbox_args, backgroundcolor = 'w')


			ax = fig1.add_subplot(6,1,5)
			hline = ax.axhline(y = fit[bestfit].x[0],linewidth=2, color='y', linestyle = '-.',   label = r'$\theta_{r}$')
			vline = ax.axvline(x = fit[bestfit].x[1],linewidth=2, color='y', linestyle = ':',   label = r'$f_{r}$')
			line = ax.plot(f,z_theta,'g-',linewidth = 3,label = 'Data')
			line = ax.plot(f,(fit[bestfit].x[0] + 2.0*np.arctan(2.0*fit[bestfit].x[2]*(1-f/fit[bestfit].x[1]))),'g:', linewidth = 1, label = 'Fit ')
			line = ax.plot(f5[~f5.mask][0],z_theta5[~z_theta5.mask][0],'ko',linewidth = 3,label = 'First Point')
			ax.set_ylabel('Angle [rad]')
			ax.legend(loc = 'best', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
			
			ax = fig1.add_subplot(6,1,6)
			vline = ax.axvline(x = fit[bestfit].x[1],linewidth=2, color='y', linestyle = ':',   label = r'$f_{r}$')
			style  = ['-','--',':','-.','+','x']; s = 0 #Cyclic iterable?
			for key in fit.keys():
				line = ax.plot(f,(z_theta - fit[key].x[0] - 2.0*np.arctan(2.0*fit[key].x[2]*(1-f/fit[key].x[1]))),'b'+style[s], linewidth = 3, label = 'Data - Fit ' + key)
				s += 1
			ax.set_ylabel('Angle [rad]')
			ax.set_xlabel('Freq [Hz]')
			ax.legend(loc = 'best', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
			plt.show()

	def fill_sweep_array(self, Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False ):
		
		if Compute_Preadout == True:
			needed = ('Atten_NA_Output', 'Atten_At_4K','Cable_Calibration')

			for quantities  in needed:				
				if  self.metadata.__dict__[quantities] == None:
					print('{0} metadate missing. Unable to compute Preadout. Setting to 0.'.format(quantities))
					Compute_Preadout = False

				Atten_NA_Output = self.metadata.Atten_NA_Output
				Atten_At_4K = self.metadata.Atten_At_4K
				k = self.metadata.Cable_Calibration

			if Fit_Resonances == False:
				print('Resonance fit not selected. Computation of Preadout_dB requires knowledge of resonance frequency and may not work.')


			if Compute_Preadout == True:
				Preadout = lambda f: k[0]*np.sqrt(f)+k[1]*f+k[2] - Atten_NA_Output - Atten_At_4K

		if Add_Temperatures == True:
			Temperature_Calibration = self.metadata.Temperature_Calibration
			if type(Temperature_Calibration) == list: 
				Temperature_Calibration = np.array(Temperature_Calibration)
				# Temperature_Calibration[:,0] is heater voltages
				# Temperature_Calibration[:,1] is temperatures voltages
				tol =  0.0005
			else:
				print('Temperature_Calibration metadata is not found or not of the correct type. Unable to add temperatures.')
				Add_Temperatures = False

			
		num_records = self.Sweep_Array.size
		for index in xrange(num_records): 
			sys.stdout.write('\r {0} of {1} '.format(index+1, num_records))
			sys.stdout.flush()

			#set current loop
			self.pick_loop(index)

			if Fit_Resonances == True:
				# Remove Gain Compression
				self.decompress_gain(Compression_Calibration_Index = -1, Show_Plot = False, Verbose = False)

				# Remove Cable Delay
				self.remove_cable_delay(Show_Plot = False, Verbose = False)	

				# Fit loop to circle
				self.circle_fit(Show_Plot = False)

				# Fit resonance parameters
				self.phase_fit(Fit_Method = 'Multiple',Verbose = False, Show_Plot = False)
				

				self._define_sweep_array(index, Q = self.loop.Q,
												Qc = self.loop.Qc,
												fr = self.loop.fr)

			if Compute_Preadout == True:
				if self.loop.fr != None:
					self._define_sweep_array(index, Preadout_dB = self.Sweep_Array['Pinput_dB'][index] + Preadout(self.loop.fr))
				elif np.abs(self.loop.freq[-1]-self.loop.freq[0]) > 1e9:
					print('Sweep bandwidth is {0} Hz. Sweep looks more like a survey. Preadout_dB is meaningless for a survey. Aborting Preadout computation... '.format(np.abs(self.loop.freq[-1]-self.loop.freq[0])))
					
				else:
					print('No resonance frquency (fr) on record for selected resonance. Estimating fr using sweep minimum.')
					fr = np.extract(np.abs(self.loop.z).min() == np.abs(self.loop.z),self.loop.freq)[0]
					self._define_sweep_array(index, Preadout_dB = self.Sweep_Array['Pinput_dB'][index] + fr)

			if Add_Temperatures == True:
				condition = (self.Sweep_Array['Heater_Voltage'][index] + tol > Temperature_Calibration[:,0]) & (self.Sweep_Array['Heater_Voltage'][index] - tol < Temperature_Calibration[:,0])
				if condition.sum() == 1:
					self.Sweep_Array['Temperature'][index] = Temperature_Calibration[condition,1][0]
				else:
					print('Unable to match unique temperature to heater voltage value for Sweep_Array[{1}]. {} matches found.'.format(index,condition.sum() ))

			# Clear out loop
			del(self.loop)
			self.loop = loop() 

	def fit_cable_loss(self, freq_range = [500e6, 1e9], Verbose = True, Show_Plot = True):
		'''produces fit to cable loss in the functional form:
		term1 + term2 + term3 = a * sqrt(f) + b * f + c
		term1 is the sum of inner and outer coaxial cable conductor losses
		term2 is due to coaxial cable dielectric loss
		term3 is a constant fudge factor
		The loss evaluates to units of dB.

		Two used this function load transmission for complete cable loop only (not amps or attens).
		Then call this function on that transmission data. This funciton creats the tuple (a,b,c,run) in 
		metadata, where run is the name of the calibration run.

		Create a function from a,b,c and it to the effect of attenuators on the input side of the cable loop.

		set freq_range = None to use full freq range	
		'''

		f   = self.loop.freq
		s21 = self.loop.z

		if freq_range == None:
			condition = f == f
		else:
			condition = (f>freq_range[0]) & (f<freq_range[1])
		
		f = np.extract(condition,f)
		s21 = np.extract(condition,s21)
		


		def obj(x,s21,f):
			a,b,c = x
			return np.square(20*np.log10(np.abs(s21)) - a*np.sqrt(f) - b*f - c).sum()	
		
		p0 = np.array([-3.0e-4,-1.0e-9 ,0.5])

		res = minimize(obj, p0, args=(s21,f), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=1e-15, callback=None, options={'disp':False, 'xtol' : 1e-6,'maxfev':1000})
		
		k = list(res.x/2.0) #devide by 2 to get one way loss
		k.append(self.metadata.Run)

		self.metadata.Cable_Calibration = self._Cable_Calibration = tuple(k)

		if Verbose == True:
			print(res)

		if Show_Plot == True:
			(fig,ax,) = self.plot_transmission(show = False)[:2]
			Cal  = lambda f: k[0]*np.sqrt(f)+k[1]*f+k[2]
			line = ax.plot(f, Cal(f)*2.0, 'r--', linewidth=3, label = 'fit - round trip')
			line = ax.plot(f, Cal(f), 'g-', linewidth=3, label = 'fit - one way')
			ax.set_xlim([freq_range[0]*0.75, freq_range[1]*1.25])
			leftvline = ax.axvline(x = freq_range[0],linewidth=2, color='k', linestyle = ':')
			rightvline = ax.axvline(x = freq_range[1],linewidth=2, color='k', linestyle = ':')
			ax.legend(loc = 'best', fontsize=10,scatterpoints =1, numpoints = 1, labelspacing = .1)
			plt.show()



		
		






























		 		













	