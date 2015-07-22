# run this code from terminal useing 'run -i Create_Library'

import KAM
reload(KAM)

temp = KAM.thermometry()
swp = KAM.sweep()

path = '/Users/miguel_daal/Desktop/Some_Data/'

Find_Temperatures = 0
Cable_Calibration = 1
Run52b = 0
Run51a = 0
Run51b = 0
Run49a = 0
Run48a = 0
Run48b = 0
Run47a = 0
Run46a = 0
Run45b = 0
Run45a = 0
Run44b = 1
Run44a = 0

sensor_ids = np.array(['S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16'])
Z3 = np.array([85.7779239,	67.8250245,	50.1762772,	33.824545,	21.2469277,	12.7273743,	7.18826586,	3.5793348,	78.4335348,	61.6228064,	45.6375047,	31.5852207,	20.407557,	12.372514,	7.10001542,	3.89603542])
Z3_dict_1 = dict(zip(sensor_ids,Z3))

Z1 = np.array([52.1559515,52.1558281,52.1560869,52.1561139,52.1483422,52.144045,52.1437542,52.0581491,51.2646322,51.2662686,51.2663971,51.2663751,51.2663306,51.266291,51.2661294, 51.2656037])
Z1_dict_1 = dict(zip(sensor_ids,Z1))

Eeff = np.array([4.19652819,	3.81480738,	3.25172475,	2.61741572,	2.09048564,	1.70180979,	1.44820567,	1.27737833,	4.0067732,	3.60877767,	3.11707476,	2.590019,	2.10964846,	1.73077079,	1.46394603,	1.29009339])
Eeff_dict_1 = dict(zip(sensor_ids,Eeff))

ground_plane_wafer = np.array(['DSAT3','DSAT4','DSAT5','DSAT6','DSAT7','GSAT1', 'HSA2'])
ground_plane_wafer_thickness = np.array([348.0,	309.0,	309.0,	164.0,	85.0,	305.5,	265.33])
gp_thickness_dict_1 = dict(zip(ground_plane_wafer,ground_plane_wafer_thickness))

def _set_metadata(Z1_dict,Z3_dict,Eeff_dict,gp_thickness_dict):
	def _extract_wafer_label(string_input):
		'''returns first word in string for which all cased characters are upper case are AND for which the word is entirely 
		alpha numeric (no '-', '(', etc) AND for which word is >= 4 characters. returns warning if there are two matches or no matches.
		'''
		
		matches = 0
		wafer_label = None
		split_string = string_input.replace('(','|').replace(')','|').replace(' ','|').split('|')[::-1] # reversed so that first occurance comes last
		for word in split_string:
			if word.isalnum() and word.isupper() and len(word) >= 4:
				wafer_label = word
				matches = matches + 1
		if matches > 1:
			warnings.warn('More than one wafer label found. Using first match.')
		if matches == 0:
			warnings.warn('No wafer label matches found. Returning None')

		return wafer_label

	def _extract_sensor_id():
		sensor = swp.metadata.Sensor
		wafer_label  = _extract_wafer_label(sensor)

		split_string  = sensor.replace(wafer_label ,'').split()
		matches = 0
		sensor_id = None
		for word in split_string:
			if word.isalnum() and word.isupper() and len(word) <= 4:
				sensor_id = word
				matches = matches + 1
		if matches > 1:
			warnings.warn('More than one sensor id found. Using first match.')
		if matches == 0:
			warnings.warn('No sensor id matches found. Returning None')

		return sensor_id

		

	def _extract_width():
		''' extracts resonator width from swp.metadata.Sensor,'''
		sensor = swp.metadata.Sensor
		e = sensor.rfind('um')
		s = e - 1
		while True:
			if sensor[s-1].isdigit() is True:
				s = s -1;
			else:
				break
		width  = int(sensor[s:e])
		return width	




	gid = _extract_wafer_label(swp.metadata.Ground_Plane)
	sid = _extract_sensor_id()
	swp.metadata.Resonator_Width = _extract_width()*1e-6 # in meters
	swp.metadata.Resonator_Thickness =  200*1e-9 #list if more than one
	swp.metadata.Resonator_Impedance = Z3_dict[sid] #in ohms
	swp.metadata.Resonator_Eeff = Eeff_dict[sid] # Resonator Dielectric Constant
	swp.metadata.Feedline_Impedance = Z1_dict[sid] #in ohms
	swp.metadata.Ground_Plane_Thickness = gp_thickness_dict[gid]*1e-9 #in meters

if Cable_Calibration:
	filename = path + 'Calibrations/Cables_140520/Channel_1_Coax_2_50mK_KIDs_Run_50a_ScanData_37mK_2014520183338_Fixed20140625.mat'
	swp.load_scandata(filename)
	index = 1 
	swp.pick_loop(index)
	swp.fit_cable_loss(freq_range = [400e6, 1e9], Verbose = False, Show_Plot = False)



if Run52b: #256 um resonator (S8) on hex ground plane
	if 1:	#Survey Sweep  - Low Freq
		filename = path + 'Run52b/Survey/Low_Frequency_Survey/52b_ScanData_37mK_201531123136.mat'
		swp.load_scandata(filename)
		swp.metadata.Num_Temperatures  = 0 #Measured temps are no good bc forgot to set thermometer bias
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)	
	if 1:	#Survey Sweep  - High Freq - cannot see fundamental, only first harmonic
		filename = path + 'Run52b/Survey/52b_ScanData_42mK_201522705010.mat'
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)		
	if 1:	#TP Sweep - 1st Harmonic Low Power, decreasing temp 
		filename = path + 'Run52b/TP_Sweep/844MHz_Resonance/Low_Power_Sweep/52b_ScanData_40mK_20153405551.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1:	#TP Sweep - 1st Harmonic high Power, decreasing temp
		filename = path + 'Run52b/TP_Sweep/844MHz_Resonance/High_Power_Sweep/52b_ScanData_40mK_20153518447.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1:	#TP Sweep - 1st Harmonic downward, low & high power, increasing temp
		filename = path + 'Run52b/TP_Sweep/844MHz_Resonance/Upward_Temp_Sweep/52b_ScanData_39mK_2015311192843.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1:	#TP Sweep - Fundamental Harmonic downward, low & high power, decreasing temp
		filename = path + 'Run52b/TP_Sweep/423MHz_Resonance/52b_ScanData_Merged_40mK_20150315.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)
	if 1: 	#P Sweep - Fundamental
		filename = path + 'Run52b/P_Sweep/423MHz_Resonance/52b_ScanData_39mK_2015312102455.mat'	
		swp.load_scandata(filename)
		swp.metadata.Num_Temperatures  = 0 #Measured temps are no good bc forgot to set thermometer bias
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)
	if 1: 	#P Sweep - first harmonic
		filename = path + 'Run52b/P_Sweep/844MHz_Resonance/52b_ScanData_39mK_20153417229.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)		



if Run51b: #listed at 32um wide should be 16um, listed at sensor wafer FHN2, should be FHN1
	if 1:	#Survey Sweep 
		filename = path + 'Run51b/Survey/51b_ScanData_42mK_2014111721582.mat'
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S12 FHN1 (16um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)		
	if 1:	#TP Sweep
		filename = path + 'Run51b/TP_Sweep/51b_ScanData_Merged_42mK_20141128.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S12 FHN1 (16um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1: 	#P Sweep
		filename = path + 'Run51b/P_Sweep/51b_ScanData_41mK_2014111851960.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S12 FHN1 (16um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	

if Run51a: #isted at sensor wafer FHN2, should be FHN1
	if 1:	#Survey Sweep 
		filename = path + 'Run51a/Survey/51a_ScanData_45mK_20141116201531.mat'
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)		
	if 1:	#TP Sweep
		filename = path + 'Run51a/TP_Sweep/51a_ScanData_41mK_2014112663253.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1: 	#P Sweep
		filename = path + 'Run51a/P_Sweep/51a_ScanData_44mK_2014111702433.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	

if Run49a: #Listed at FHN4 change to FHN2!
	Run49a_Temps =  [(0.7, 0.602), (0.65,0.559), (0.6,0.502), (0.55,0.4456), (0.5,0.375), (0.45,0.3289), (0.4,0.290), (0.35,0.2524), (0.3, 0.2167), (0.25,0.184), (0.225,0.165), (0.2,0.1484), (0.175,0.1324), (0.15,0.1158), (0.1,0.0802), (0.0,0.0372)]
	if Find_Temperatures:
		filename = path + 'Run49a/TP_Sweep/Run49a_Temperatures'
		temp.load_MonitoringVI_file(filename,temp_list = Run49a_Temps)
	swp._Temperature_Calibration = Run49a_Temps

	if 1:	#Survey Sweep 
		filename = path + 'Run49a/Survey/49a_ScanData_45mK_2014417175147.mat'
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S2 FHN2 (4um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1:	#TP Sweep
		filename = path + 'Run49a/TP_Sweep/49a_ScanData_45mK_2014420142743.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S2 FHN2 (4um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep
		filename = path + 'Run49a/P_Sweep/49a_ScanData_45mK_201442443323.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S2 FHN2 (4um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

if Run48a:
	Run48a_Temps = [(0.65, 0.653), (0.6, 0.63), (0.55,0.575), (0.5, 0.535), (0.45, 0.498), (0.4,0.455), (0.35,0.395), (0.3,0.325), (0.25, 0.258), (0.2,.206), (0.15, 0.154), (0.1,0.106), (0.0, 0.051)]
	if Find_Temperatures:
		filename = path + 'Run48a/TP_Sweep/Run48a'
		temp.load_MonitoringVI_file(filename, temp_list = Run48a_Temps)
	swp._Temperature_Calibration = Run48a_Temps
	
	if 1:	#Survey Sweep 
		filename = path + 'Run48a/Survey/48a_ScanData_50mK_2014331143139.mat'
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S3 FHN1 (8um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1:	#TP Sweep
		filename = path + 'Run48a/TP_Sweep/48a_ScanData_70mK_201433101142.mat'	
		swp.load_scandata(filename)
		swp.Sweep_Array['Pinput_dB'] = np.round(swp.Sweep_Array['Pinput_dB']) #NA Rounding Error not fixed yet
		swp.metadata.Sensor = 'S3 FHN1 (8um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep
		filename = path + 'Run48a/P_Sweep/48a_ScanData_45mK_201433112722.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S3 FHN1 (8um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

if Run48b:
	# Temperature data incomplete, No date for Heater_Voltage: 0.25, 0.15, 0.1
	Run48b_Temps = [(0.65,0.553), (0.55, 0.52), (0.45, 0.468), (0.35,0.390),(0.2,0.261), (0.25,0.0), (0.2,0.261), (0.15, 0.0), (0.1,0.0), (0.0, 0.045)]
	if Find_Temperatures:
		filename = path + 'Run48b/TP_Sweep/Run48b_Crashed_Midway_Use_48a' # Run48a temp data is inconsisten with Run48b temp data e.g. Run48a has (0.2,.206), while Run48b has (0.2,0.261)
		temp.load_MonitoringVI_file(filename, temp_list = Run48b_Temps)
	swp._Temperature_Calibration = Run48b_Temps

	if 1:	#Survey Sweep 
		filename = path + 'Run48b/Survey/48b_ScanData_45mK_201441123255.mat'
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1:	#TP Sweep
		filename = path + 'Run48b/TP_Sweep/48b_ScanData_45mK_201441101927.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep
		filename = path + 'Run48b/P_Sweep/48b_ScanData_45mK_201441132117.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

if Run47a: #this run failed
	if 1:	#Survey Sweep 
		filename = path + 'Run47a/Survey/47a_ScanData_50mK_20143234307.mat'
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

if Run46a:
	#Temps noted in Run logbook, no temperature file recorded this run
	Run46a_Temps = [(0.9,0.629), (0.75,0.485),(0.6,0.385),(0.5,0.328),(0.4,0.2745),(0.3,0.2085),(0.2,0.152),(0.1,0.0935), (0.0,0.044)]
	swp._Temperature_Calibration = Run46a_Temps

	if 1:	#Survey Sweep 
		filename = path + 'Run46a/Survey/46a_ScanData_49mK_2014221194941.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep - Low Power
		filename = path + 'Run46a/P_Sweep/Low_Power_Sweep/46a_ScanData_49mK_201422114464.mat'	
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep - High Power - Signicifant Amp Saturation
		filename = path + 'Run46a/P_Sweep/High_Power_Sweep_w_20dB_Atten_Added/46a_ScanData_49mK_2014221141745.mat'	
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 20
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #TP Sweep
		# Due to persistent crashes of MonitoringVI, we piece together temps from three different files
		filename = path + 'Run46a/TP_Sweep/Prev_Attempt_500mV_to_200mV/ScanData_Backup.mat'
		swp.load_scandata(filename)
		Sweep_Array_500mV_to_200mV = swp.Sweep_Array
		filename = path + 'Run46a/TP_Sweep/Prev_Attempt_100mV_to_000mV/46a_ScanData_39mK_2014228232848.mat'	
		swp.load_scandata(filename)
		Sweep_Array_100mV_to_000mV = swp.Sweep_Array		
		filename = path + 'Run46a/TP_Sweep/Prev_Attempt_900mV_to_600mV/ScanData_Backup.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		ind = swp.Sweep_Array.size- Sweep_Array_500mV_to_200mV.size
		swp.Sweep_Array[ind:] = Sweep_Array_500mV_to_200mV
		ind = swp.Sweep_Array.size- Sweep_Array_100mV_to_000mV.size
		swp.Sweep_Array[ind:] = Sweep_Array_100mV_to_000mV
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

if Run45b: # Had gross Microphonics
	if 1:	#Survey Sweep 
		filename = path + 'Run45b/Survey/45b_ScanData_40mK_20142622445.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'InP #2'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay(Show_Plot = False, Verbose = True, center_freq = None, Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep  - Feb 6 - like Feb 4 but extended to higher power 
		filename = path + 'Run45b/P_Sweep/Feb6/45b_ScanData_47mK_201426213523.mat'	
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'InP #2'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 10
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep - Feb 4
		filename = path + 'Run45b/P_Sweep/Feb4/45b_ScanData_44mK_201424175218.mat'	
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'InP #2'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #TP Sweep 
		# No TP taken
		pass	

if Run45a:
	#Temp for heat voltage 0.5 is very badly defined due to large downward slip (0.433 to 0.309mK) of temp during scan. should remove that from data
	# Due to thermometry crashes, no data for heater voltage 0.2, 0.1 and 0.0 for TP sweep
	# Using base temperature from  power sweep  for heater voltage 0.0
	Run45a_Temps = [(0.9,0.61), (0.8, 0.590), (0.7, 0.524), (0.6,0.467),(0.5,0.395), (0.4,0.248), (0.3,0.210), (0.2, 0), (0.1,0), (0.0,0.043)]
	if Find_Temperatures:
		filename = path + 'Run45a/TP_Sweep/Run45a_Temp_140205' 
		temp.load_MonitoringVI_file(filename, process_therm = 2,temp_list = Run45a_Temps)
	swp._Temperature_Calibration = Run45a_Temps

	if 1:	#Survey Sweep 
		filename = path + 'Run45a/Survey/45a_ScanData_48mK_201426222214.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay(Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1: #TP Sweep 
		# Data not valid for heater voltage 0.5, 0.2 and 0.1
		# Data missing for heater voltages 0.2, 0.1 and 0.0 for TP sweep	
		filename = path + 'Run45a/TP_Sweep/ScanData_Backup.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		swp.Sweep_Array = swp.Sweep_Array[:5*21] # remove zeros from Sweep_Array which correspond to 0.2 and 0.1 and 0.0
		swp.metadata.Num_Heater_Voltages = 7
		for rec in swp.Sweep_Array:
			if rec['Heater_Voltage'] == 0.5:
				rec['Is_Valid'] = False
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)
	
	if 1: #P Sweep  - Feb 6 - again, like Feb 4 but extended to higher power 
		filename = path + 'Run45a/P_Sweep/Feb6/45a_ScanData_50mK_201426175813.mat'	
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 10
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep - Feb 4
		filename = path + 'Run45a/P_Sweep/Feb4/45a_ScanData_44mK_201424185419.mat'	
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 40		
		swp.metadata.Atten_NA_Output = 0
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

if Run44b:
	Run44b_Temps = [(0.0,0.0495), (0.1, 0.0827), (0.2, 0.156), (0.3,0.225), (0.4,0.280), (0.5,0.346), (0.6, 0.412), (0.7,0.476), (0.8,0.544)]
	if Find_Temperatures:
		filename = path + 'Run44b/TP_Sweep/Run44b_131218' 
		temp.load_MonitoringVI_file(filename, process_therm = 1,temp_list = Run44b_Temps)
	swp._Temperature_Calibration = Run44b_Temps

	if 1:	#Survey Sweep 
		filename = path + 'Run44b/Survey/44b_ScanData_40mK_20131218112528.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.LNA['Vd'] = 2.5
		swp.metadata.LNA['Id'] = 0.010
		swp.metadata.LNA['Vg'] = 0.0
		swp.metadata.Atten_At_4K = 30		
		swp.metadata.Atten_NA_Output = 20
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1: #TP Sweep 	
		filename = path + 'Run44b/TP_Sweep/44b_ScanData_40mK_201312196158.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.Atten_At_4K = 30		
		swp.metadata.Atten_NA_Output = 20
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 0: #P Sweep 
		filename = path + 'Run44b/P_Sweep/44b_ScanData_60mK_2013121105756.mat'	
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'SiGe #1'
		swp.metadata.LNA['Vd'] = 2.5
		swp.metadata.LNA['Id'] = 0.010
		swp.metadata.LNA['Vg'] = 0.0
		swp.metadata.Atten_At_4K = 30		
		swp.metadata.Atten_NA_Output = 20
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)


if Run44a: 
	Run44a_Temps = [(1.1,0.769), (1.0,0.705), (0.9,0.628), (0.8,0.555), (0.7,0.481),  (0.6,0.414 ),  (0.5,0.347),  (0.4,0.280), (0.3,0.227), (0.2,0.1677),(0.1,0.0933 ),(0.0,0.0495)]
	if Find_Temperatures:
		filename = path + 'Run44a/TP_Sweep/Run44a_131219_Fix140627' 
		temp.load_MonitoringVI_file(filename, process_therm = 1,temp_list = Run44a_Temps)
	swp._Temperature_Calibration = Run44a_Temps

	if 1:	#Survey Sweep 
		filename = path + 'Run44a/Survey/44a_ScanData_40mK_2013121819025.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'InP #2'
		swp.metadata.LNA['Vd'] = 1.0
		swp.metadata.LNA['Id'] = 0.015
		swp.metadata.LNA['Vg'] = 1.9
		swp.metadata.Atten_At_4K = 30		
		swp.metadata.Atten_NA_Output = 20
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1: #TP Sweep 	
		filename = path + 'Run44a/TP_Sweep/44a_ScanData_40mK_20131220154159.mat'
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'InP #2'
		swp.metadata.Atten_At_4K = 30		
		swp.metadata.Atten_NA_Output = 20
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep 
		filename = path + 'Run44a/P_Sweep/44a_ScanData_60mK_20131210213718.mat'	
		swp.load_scandata(filename)
		swp.metadata.LNA['LNA'] = 'InP #2'
		swp.metadata.LNA['Vd'] = 1.0
		swp.metadata.LNA['Id'] = 0.015
		swp.metadata.LNA['Vg'] = 1.9
		swp.metadata.Atten_At_4K = 30		
		swp.metadata.Atten_NA_Output = 20
		swp.metadata.Atten_NA_Input = 0
		swp.metadata.Atten_RTAmp_Input = 0
		swp.metadata.RTAmp_In_Use = 1
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)


