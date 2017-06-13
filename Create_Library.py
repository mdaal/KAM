# run this code from terminal useing 'run -i Create_Library'
import time
import KAM
import numpy as np
reload(KAM)

start = time.time()

temp = KAM.thermometry()
swp = KAM.sweep()

path = '/Users/miguel_daal/Desktop/Some_Data/'

Find_Temperatures = 0
Cable_Calibration = 1
Run53ab = 1
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
Run44b = 0
Run44a = 0
Save_Sonnet_Sims = 0
Mock_Data = 0

#############
#
#
#       Calibration information for redout devices incl. noise and gain of amplifiers and NA
#
#
#############

System_Calibration = {}
# Units:
# F = Hz
# G = dB
# Nm = Kelvin amplitude noise
# Np = Kelvin phase noise
# P1dB = dBm at input

F = [409600000., 509500000, 609400000, 709300000, 809200000, 909100000, 1009000000]
G = [48.2976, 48.1631, 48.0342, 47.684, 47.3322, 46.9015, 46.4087]
Nm = [5.48601,5.72677,5.77245,6.0146,5.94466,5.95258,6.13591]
Np = Nm
P1dB = -48.
System_Calibration['SiGe #2'] = dict(freq =F,g= G,Tn_m = Nm,Tn_p = Np,P1dB= P1dB)

F = [400000000., 500000000, 600000000, 700000000, 800000000, 900000000, 1000000000]
G = [47.5, 47.3, 47.0, 46.6, 46.0, 45.6, 45.0]
Nm = [5.7, 6.3, 6.1, 6.3, 6.6, 6.6, 7.0]
Np = Nm
P1dB = -48.
System_Calibration['SiGe #1'] = dict(freq =F,g= G,Tn_m = Nm,Tn_p = Np,P1dB= P1dB)

F = [400000000., 500000000, 600000000, 700000000, 800000000, 900000000, 1000000000]
G = [31.6958, 31.7332, 33.3491, 34.582, 36.4723, 35.4597, 35.4288]
Nm = [10.0454, 9.13519, 8.59369, 8.28292, 8.11164, 7.86398, 7.5084]
Np = Nm
P1dB = -40.
System_Calibration['InP #2'] = dict(freq =F,g= G,Tn_m = Nm,Tn_p = Np,P1dB= P1dB)

F = [395000000.,424500000, 454000000, 483500000, 513000000, 542500000, 572000000, 601500000, 631000000, 660500000, 690000000, 719500000, 749000000, 778500000, 808000000, 837500000, 867000000, 896500000, 926000000, 955500000, 985000000]
Nm = [488.6053076, 489.2697104, 474.3872914, 473.1617407, 467.3813032, 470.0790808, 460.311823, 446.0557503, 439.5120869, 439.2718329, 429.4232912, 423.771945, 414.857644, 419.2513872, 410.1854824, 401.0360685, 398.6567012, 400.5274794, 389.6199922, 390.1977107, 411.6113207]
G = [ 39.1,  39.1,  39.1,  39. ,  39. ,  39. ,  39. ,  38.9,  38.9, 38.9,  38.8,  38.7,  38.7,  38.6,  38.6,  38.5,  38.4,  38.4, 38.3,  38.2,  38.1]
Np = Nm
P1dB = -11.2 # at input
System_Calibration['AML016P3411'] = dict(freq =F,g= G,Tn_m = Nm,Tn_p = Np,P1dB= P1dB)

F = [3000000., 3000000000]
G = [1.,1]
Nm = [1448594313,1448594313] ### np.square(2e-6)/(4*k*50)   <--  is Tn (Hz)^2, must devide by BW^2 to get Tn
Np =  [1752799118,1752799118]### np.square(2.2e-6)/(4*k*50)
P1dB = +20. # Fault power level. 
System_Calibration['NA']  = dict(freq =F,g= G,Tn_m = Nm,Tn_p = Np,P1dB= P1dB)

#############
#
#
#       Impedance, width, and Eeff data  for resonators
#
#
#############
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


#############
#
#
#       Code for connecting impedance, width, and Eeff  and readout calibration data to data sets
#
#
#############

def _set_metadata(Z1_dict,Z3_dict,Eeff_dict,gp_thickness_dict,System_Calibration):
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
	swp.metadata.System_Calibration = System_Calibration
	swp.fit_system_calibration()
	swp.metadata.Digitizer = 'NA'
	swp.metadata.RTAmp = 'AML016P3411'

if Cable_Calibration:
	filename = path + 'Calibrations/Cables_140520/Channel_1_Coax_2_50mK_KIDs_Run_50a_ScanData_37mK_2014520183338_Fixed20140625.mat'
	swp.load_scandata(filename)
	index = 1 
	swp.pick_loop(index)
	swp.fit_cable_loss('One_Way_40mK',freq_range = [400e6, 1e9], Verbose = False, Show_Plot = False)
	filename = path + 'Calibrations/Cables_111016/Channel_2_4K/Channel 2 4K_ScanData_300000mK_2011101616135.mat'
	swp.load_scandata(filename)
	index = 0 
	swp.pick_loop(index)
	swp.fit_cable_loss('One_Way_4K',freq_range = [400e6, 1e9], Verbose = False, Show_Plot = False)
	filename = path + 'Calibrations/Cables_110908/Flexible_Coax_Loop_300K/Full Loop_ScanData_300000mK_20119812370.mat'
	swp.load_scandata(filename)
	index = 0 
	swp.pick_loop(index)
	swp.fit_cable_loss('One_Way_300K',freq_range = [400e6, 1e9], Verbose = False, Show_Plot = False)
	swp.metadata.Cable_Calibration['300K_to_4K'] = (swp.metadata.Cable_Calibration['One_Way_4K'][0]  - swp.metadata.Cable_Calibration['One_Way_300K'][0],
													swp.metadata.Cable_Calibration['One_Way_4K'][1]  - swp.metadata.Cable_Calibration['One_Way_300K'][1], 
													swp.metadata.Cable_Calibration['One_Way_4K'][2]  - swp.metadata.Cable_Calibration['One_Way_300K'][2],
													'One_Way_4K - One_Way_300K', 400995077.27409261, 998491719.05006254)
	swp.metadata.Cable_Calibration['4K_to_40mK'] = (swp.metadata.Cable_Calibration['One_Way_40mK'][0]  - swp.metadata.Cable_Calibration['One_Way_4K'][0],
													swp.metadata.Cable_Calibration['One_Way_40mK'][1]  - swp.metadata.Cable_Calibration['One_Way_4K'][1],
													swp.metadata.Cable_Calibration['One_Way_40mK'][2]  - swp.metadata.Cable_Calibration['One_Way_4K'][2],
													'One_Way_40mK - One_Way_4K', 400995077.27409261, 998491719.05006254)

if Run53ab: #256 um resonator (S8) in mu-metal sheild
	if 1: #Survey Sweep  - Low Freq
		filename = path + 'Run53ab/Run_53ab_Data.h5'
		table_path = '/SweepP/T201603112022_SweepP'
		swp.load_hf5_2(filename, table_path)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)	
	if 0: 	#P Sweep 
		filename = path + 'Run53ab/Run_53ab_Data.h5'
		table_path = '/SweepNP/T201603130517_SweepNP'
		swp.load_hf5_2(filename, table_path)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)
		Power_Calibration  = swp.metadata.Power_Calibration # forgot to  include in TP sweep below
	if 1:	#TP Sweep 
		filename = path + 'Run53ab/Run_53ab_Data.h5'
		table_path = '/SweepNPT/T201603120024_SweepNPT'
		swp.load_hf5_2(filename, table_path)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.metadata.Power_Calibration = Power_Calibration 
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)



if Run52b: #256 um resonator (S8) on hex ground plane
	if 1:	#Survey Sweep  - Low Freq
		filename = path + 'Run52b/Survey/Low_Frequency_Survey/52b_ScanData_37mK_201531123136.mat'
		swp.load_scandata(filename)
		swp.metadata.Num_Temperatures  = 0 #Measured temps are no good bc forgot to set thermometer bias
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)	
	if 1:	#Survey Sweep  - High Freq - cannot see fundamental, only first harmonic
		filename = path + 'Run52b/Survey/52b_ScanData_42mK_201522705010.mat'
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.System_Calibration = System_Calibration
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)		
	if 1:	#TP Sweep - 1st Harmonic Low Power, decreasing temp 
		filename = path + 'Run52b/TP_Sweep/844MHz_Resonance/Low_Power_Sweep/52b_ScanData_40mK_20153405551.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1:	#TP Sweep - 1st Harmonic high Power, decreasing temp
		filename = path + 'Run52b/TP_Sweep/844MHz_Resonance/High_Power_Sweep/52b_ScanData_40mK_20153518447.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1:	#TP Sweep - 1st Harmonic downward, low & high power, increasing temp
		filename = path + 'Run52b/TP_Sweep/844MHz_Resonance/Upward_Temp_Sweep/52b_ScanData_39mK_2015311192843.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1:	#TP Sweep - Fundamental Harmonic downward, low & high power, decreasing temp
		filename = path + 'Run52b/TP_Sweep/423MHz_Resonance/52b_ScanData_Merged_40mK_20150315.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)
	if 1: 	#P Sweep - Fundamental
		filename = path + 'Run52b/P_Sweep/423MHz_Resonance/52b_ScanData_39mK_2015312102455.mat'	
		swp.load_scandata(filename)
		swp.metadata.Num_Temperatures  = 0 #Measured temps are no good bc forgot to set thermometer bias
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)
	if 1: 	#P Sweep - first harmonic
		filename = path + 'Run52b/P_Sweep/844MHz_Resonance/52b_ScanData_39mK_20153417229.mat'	
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)		



if Run51b: #listed at 32um wide should be 16um, listed at sensor wafer FHN2, should be FHN1
	if 1:	#Survey Sweep 
		filename = path + 'Run51b/Survey/51b_ScanData_42mK_2014111721582.mat'
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S12 FHN1 (16um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)		
	if 1:	#TP Sweep
		filename = path + 'Run51b/TP_Sweep/51b_ScanData_Merged_42mK_20141128.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S12 FHN1 (16um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1: 	#P Sweep
		filename = path + 'Run51b/P_Sweep/51b_ScanData_41mK_2014111851960.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S12 FHN1 (16um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	

if Run51a: #isted at sensor wafer FHN2, should be FHN1
	if 1:	#Survey Sweep 
		filename = path + 'Run51a/Survey/51a_ScanData_45mK_20141116201531.mat'
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)		
	if 1:	#TP Sweep
		filename = path + 'Run51a/TP_Sweep/51a_ScanData_41mK_2014112663253.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)	
	if 1: 	#P Sweep
		filename = path + 'Run51a/P_Sweep/51a_ScanData_44mK_2014111702433.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		swp.metadata.LNA['LNA'] =  'SiGe #1'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1:	#TP Sweep
		filename = path + 'Run49a/TP_Sweep/49a_ScanData_45mK_2014420142743.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S2 FHN2 (4um width)'
		swp.metadata.LNA['LNA'] =  'SiGe #1'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep
		filename = path + 'Run49a/P_Sweep/49a_ScanData_45mK_201442443323.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S2 FHN2 (4um width)'
		swp.metadata.LNA['LNA'] =  'SiGe #1'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		swp.metadata.LNA['LNA'] =  'SiGe #1'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 0:	#TP Sweep
		filename = path + 'Run48a/TP_Sweep/48a_ScanData_70mK_201433101142.mat'	
		swp.load_scandata(filename)
		swp.Sweep_Array['Pinput_dB'] = np.round(swp.Sweep_Array['Pinput_dB']) #NA Rounding Error not fixed yet
		swp.metadata.Sensor = 'S3 FHN1 (8um width)'
		swp.metadata.LNA['LNA'] =  'SiGe #1'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep
		filename = path + 'Run48a/P_Sweep/48a_ScanData_45mK_201433112722.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S3 FHN1 (8um width)'
		swp.metadata.LNA['LNA'] =  'SiGe #1'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.pick_loop(0)
		Electrical_Delay = swp.remove_cable_delay( Show_Plot = False, Verbose = True, center_freq = None,Force_Recalculate = True); 
		swp.fill_sweep_array(Fit_Resonances = False, Compute_Preadout = False, Add_Temperatures = True) #Dont Compute Preadout for Surverys. That would be meaningless.
		swp.save_hf5(overwrite = True)

	if 1:	#TP Sweep
		filename = path + 'Run48b/TP_Sweep/48b_ScanData_45mK_201441101927.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep
		filename = path + 'Run48b/P_Sweep/48b_ScanData_45mK_201441132117.mat'	
		swp.load_scandata(filename)
		swp.metadata.Sensor = 'S5 FHN1 (32um width)'
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

if Run47a: #this run failed
	if 1:	#Survey Sweep 
		filename = path + 'Run47a/Survey/47a_ScanData_50mK_20143234307.mat'
		swp.load_scandata(filename)
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)

	if 1: #P Sweep 
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
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
		_set_metadata(Z1_dict_1,Z3_dict_1,Eeff_dict_1,gp_thickness_dict_1,System_Calibration)
		swp.metadata.Electrical_Delay = Electrical_Delay
		swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = True, Add_Temperatures = True)
		swp.save_hf5(overwrite = True)




if Save_Sonnet_Sims:
	# Add Simulation data to Database
	path = '/Users/miguel_daal/Desktop/Some_Data/'
	swp = KAM.sweep()

	filename = path + 'Resonator_Simulations/S3.s2p'
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False)
	swp.metadata.Resonator_Width = 2.**3 * 1e-6
	swp.metadata.Run = 'S3_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S3_no_teflon.s2p' #This simulation replaces the teflon in the dielectric with sapphire.
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False)
	swp.metadata.Resonator_Width = 2.**3 * 1e-6
	swp.metadata.Run = 'S3_sapphire_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S3_no_teflon_v16.s2p' #This simulation replaces the teflon in the dielectric with sapphire.
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False)
	swp.metadata.Resonator_Width = 2.**3 * 1e-6
	swp.metadata.Run = 'S3_sapphire_v16_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S4_no_teflon.s2p' #This simulation replaces the teflon in the dielectric with sapphire.
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False)
	swp.metadata.Resonator_Width = 2.**4 * 1e-6
	swp.metadata.Run = 'S4_sapphire_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S4.s2p'
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False) 
	swp.metadata.Resonator_Width = 2.**4 * 1e-6
	swp.metadata.Run = 'S4_Sim'
	swp.save_hf5(overwrite = True)


	filename = path + 'Resonator_Simulations/S5_no_teflon.s2p' #This simulation replaces the teflon in the dielectric with sapphire.
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False)
	swp.metadata.Resonator_Width = 2.**5 * 1e-6
	swp.metadata.Run = 'S5_sapphire_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S5.s2p'
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False) 
	swp.metadata.Resonator_Width = 2.**5 * 1e-6
	swp.metadata.Run = 'S5_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S6.s2p'
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False) 
	swp.metadata.Resonator_Width = 2.**6 * 1e-6
	swp.metadata.Run = 'S6_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S7.s2p'
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False) 
	swp.metadata.Resonator_Width = 2.**7 * 1e-6
	swp.metadata.Run = 'S7_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S5_no_teflon.s2p' #This simulation replaces the teflon in the dielectric with sapphire.
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False)
	swp.metadata.Resonator_Width = 2.**8 * 1e-6
	swp.metadata.Run = 'S5_sapphire_Sim'
	swp.save_hf5(overwrite = True)

	filename = path + 'Resonator_Simulations/S8.s2p'
	swp.load_touchstone(filename)
	swp.metadata.Electrical_Delay = 0
	swp.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False, Complete_Fit = False) 
	swp.metadata.Resonator_Width = 2.**8 * 1e-6
	swp.metadata.Run = 'S8_Sim'
	swp.pick_loop(0)
	swp.trim_loop(N =2)
	swp.circle_fit(Show_Plot = False)
	swp.phase_fit(Fit_Method = 'Multiple',Verbose = False, Show_Plot = False)
	swp._define_sweep_array(0 , Q = swp.loop.Q,
								Qc = swp.loop.Qc,
								Fr = swp.loop.fr,
								#Mask = swp.loop.phase_fit_mask,
								Chi_Squared = swp.loop.chisquare,
								R = swp.loop.R,
								r = swp.loop.r,
								a = swp.loop.a,
								b = swp.loop.b,
								Theta = swp.loop.theta,
								Phi = swp.loop.phi)
	swp.save_hf5(overwrite = True)

if Mock_Data:
	import KAM
	reload(KAM)	
	database_location = 'Data/My_Data_Library.h5'
	Run45aP     = KAM.sweep(); Run45aP.load_hf5('/Run45a/T201402061732', filename = database_location);
	Run48aP     = KAM.sweep(); Run48aP.load_hf5('/Run48a/T201403311119', filename = database_location);

	def Gen_Mock_Data_Sweep(Sweep, Indexing = (None,None,None)):
		index = 0 
		Sweep.pick_loop(index)
		fit, fig, ax = Sweep.nonlinear_fit(Save_Fig = True, Indexing = Indexing)
		Sweep.pick_loop(index)
		

		#Construct Simulated data for run  Sweep, then non-linear fit it, then construct sweep array including Concurrent and Stepwise fits and save
		bestfit = 'Powell'
		Zfl = Sweep.metadata.Feedline_Impedance
		Zres = Sweep.metadata.Resonator_Impedance
		V30V30 = fit['V30V30']
		eta = fit[bestfit].x[4]
		delta =  fit[bestfit].x[5]
		f0 = fit[bestfit].x[0]
		Qi = fit[bestfit].x[1]
		Qc = fit[bestfit].x[2]
		phi31 = fit[bestfit].x[3]
		# use noise calculated from complete fit
		#Amplitude_Noise_Variance = cfit['sigma_squared_m'][0]
		#Phase_Noise_Variance = cfit['sigma_squared_p'][0]

	
		sweep  = KAM.sweep();
		sweep.generate_nonlinear_data(Show_Plot = True, Phase_Noise_Variance = None, Amplitude_Noise_Variance = None, Like = Sweep, Save_Fig = True,
		curve_parameter_dict = {'f_0':f0, 'Qtl':Qi, 'Qc':Qc, 'eta':eta, 'delta':delta, 'Zfl':Zfl, 'Zres':Zres, 'phi31': phi31, 'phiV1':0, 'V30V30':V30V30 },
		sweep_parameter_dict = {'Run': 'Mock_' + Sweep.metadata.Run, 'Pprobe_dBm_Start' :-100.0,'Pprobe_dBm_Stop': -64.0, 'Pprobe_Num_Points':17, 'numBW':40,'num': 2000, 'Up_or_Down': 'Up', 'Freq_Spacing':'Linear'})

		sweep.nonlinear_fit(Save_Fig = True,Indexing = (None,None,None))
		sweep.fill_sweep_array(Fit_Resonances = True, Compute_Preadout = False, Add_Temperatures = False,Complete_Fit = True )
		sweep.save_hf5(overwrite = True)
		return sweep 
	#sweep = 	Gen_Mock_Data_Sweep(Run45aP, Indexing = (None,-1,None))#(None,-1,None) *********'Pprobe_dBm_Stop': -54.0,

	sweep = 	Gen_Mock_Data_Sweep(Run48aP, Indexing = (None,-18,2))

finished = time.time()
elapsed = (finished - start )/60.0 #minutes
print 'Library Creation took {:.2f} minutes'.format(elapsed)
