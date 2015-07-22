import KAM
reload(KAM)

#syntax .....metadata.custom = [show on plot?, add this text to legend, plot linestyle]

database_location = 'Data/My_Data_Library.h5'

if Ver == 1: #power sweeps
	Run52b1st	= KAM.sweep(); Run52b1st.load_hf5('/Run52b/T201503041056',filename = database_location); Run52b1st.metadata.width = 256; Run52b1st.metadata.custom = [True, '; hex gnd', '--',None]
	Run52b2nd	= KAM.sweep(); Run52b2nd.load_hf5('/Run52b/T201503112352',filename = database_location); Run52b2nd.metadata.width = 256; Run52b2nd.metadata.custom = [False,'; hex gnd; $f_1$', ':',None]
	Run51aP	    = KAM.sweep(); Run51aP.load_hf5('/Run51a/T201411162028',filename = database_location); Run51aP.metadata.width = 32;	Run51aP.metadata.custom = [True,'', '-',None]
	Run51bP     = KAM.sweep(); Run51bP.load_hf5('/Run51b/T201411180125',filename = database_location); Run51bP.metadata.width = 16; Run51bP.metadata.custom = [True,'', '-',None]
	Run49aP     = KAM.sweep(); Run49aP.load_hf5('/Run49a/T201404232144', filename = database_location); Run49aP.metadata.width = 4; Run49aP.metadata.custom = [True,'', '-',None]
	Run48bP     = KAM.sweep(); Run48bP.load_hf5('/Run48b/T201404011305', filename = database_location); Run48bP.metadata.width = 32; Run48bP.metadata.custom = [True,'', '-',None]
	Run48aP     = KAM.sweep(); Run48aP.load_hf5('/Run48a/T201403311119', filename = database_location); Run48aP.metadata.width = 8; Run48aP.metadata.custom = [True,'', '-',None]
	Run46aPl    = KAM.sweep(); Run46aPl.load_hf5('/Run46a/T201402211440', filename = database_location); Run46aPl.metadata.width = 256; Run46aPl.metadata.custom = [True,'', '-',None]
	Run46aPh    = KAM.sweep(); Run46aPh.load_hf5('/Run46a/T201402211413', filename = database_location); Run46aPh.metadata.width = 256; Run46aPh.metadata.custom = [True,'', '-',None]# high power
	Run45bP     = KAM.sweep(); Run45bP.load_hf5('/Run45b/T201402062055', filename = database_location); Run45bP.metadata.width = 128; Run45bP.metadata.custom = [True,';70nm gnd; loose torq', '-',None]# to high power
	Run45aP     = KAM.sweep(); Run45aP.load_hf5('/Run45a/T201402061732', filename = database_location); Run45aP.metadata.width =64; Run45aP.metadata.custom = [True,'', '-',None]# to high power
	Run44bP     = KAM.sweep(); Run44bP.load_hf5('/Run44b/T201312102229', filename = database_location); Run44bP.metadata.width = 32; Run44bP.metadata.custom = [True,'; Cu box', '-', None]
	Run44aP     = KAM.sweep(); Run44aP.load_hf5('/Run44a/T201312101928', filename = database_location); Run44aP.metadata.width = 128; Run44aP.metadata.custom = [True,';70nm gnd', '-', None]
	
	#Run49aP.metadata.Sensor = 'S2 FHN2 (4um width)'
	#Run51bP.metadata.Sensor = 'S12 FHN1 (16um width)'
	#Run51aP.metadata.Sensor = 'S5 FHN1 (32um width)'

	Power_Sweeps = [Run51aP,Run51bP,Run49aP,Run48bP,Run48aP,Run46aPl,Run46aPh,Run45bP,Run45aP,Run44bP,Run44aP,Run52b1st,Run52b2nd]

	#needs to be incorporated in to  load_scandata or fill_sweep_array
	# for sweep in Power_Sweeps:
	# 	if sweep.metadata.Atten_Added_At_NA != 0.0:
	# 		sweep.Sweep_Array['Preadout_dB']  = sweep.Sweep_Array['Preadout_dB']-sweep.metadata.Atten_Added_At_NA

if Ver == 1: #temperature sweeps
	Run52bTP2ndh = KAM.sweep(); Run52bTP2ndh.load_hf5('/Run52b/T201503041901',filename = database_location); Run52bTP2ndh.metadata.width = 256
	Run52bTP2ndl = KAM.sweep(); Run52bTP2ndl.load_hf5('/Run52b/T201503021724',filename = database_location); Run52bTP2ndl.metadata.width = 256
	Run52bTP1st = KAM.sweep(); Run52bTP1st.load_hf5('/Run52b/T201503121941',filename = database_location); Run52bTP1st.metadata.width = 256
	Run51aTP = KAM.sweep(); Run51aTP.load_hf5('/Run51a/T201411241956',filename = database_location); Run51aTP.metadata.width = 32	
	Run51bTP = KAM.sweep(); Run51bTP.load_hf5('/Run51b/T201411212056',filename = database_location); Run51bTP.metadata.width = 16
	Run49aTP = KAM.sweep(); Run49aTP.load_hf5('/Run49a/T201404181957', filename = database_location); Run49aTP.metadata.width = 4
	Run48bTP = KAM.sweep(); Run48bTP.load_hf5('/Run48b/T201403312220', filename = database_location); Run48bTP.metadata.width = 32
	Run48aTP = KAM.sweep(); Run48aTP.load_hf5('/Run48a/T201403292355', filename = database_location); Run48aTP.metadata.width = 8
	Run46aTP = KAM.sweep(); Run46aTP.load_hf5('/Run46a/T201402262113', filename = database_location); Run46aTP.metadata.width = 256
	Run45aTP = KAM.sweep(); Run45aTP.load_hf5('/Run45a/T201402052038', filename = database_location); Run45aTP.metadata.width =64# to high power #has invalid data
	Run44bTP = KAM.sweep(); Run44bTP.load_hf5('/Run44b/T201312182102', filename = database_location); Run44bTP.metadata.width = 32
	Run44aTP = KAM.sweep(); Run44aTP.load_hf5('/Run44a/T201312191950', filename = database_location); Run44aTP.metadata.width = 128 

	#Run49aTP.metadata.Sensor = 'S2 FHN2 (4um width)'
	#Run51bTP.metadata.Sensor = 'S12 FHN1 (16um width)'
	#Run51aTP.metadata.Sensor = 'S5 FHN1 (32um width)'	


	TP_Sweeps = [Run51aTP,Run51bTP,Run48bTP,Run49aTP,Run48aTP,Run46aTP,Run45aTP,Run44bTP,Run44aTP,Run52bTP2ndh,Run52bTP2ndl,Run52bTP1st]




#######################################################
##
##              +  VERSION 2  +    
##
##
#######################################################

if Ver == 2: #power sweeps
	Run52b1st	= KAM.sweep(); Run52b1st.load_hf5('/Run52b/T201503041056', filename = database_location); Run52b1st.metadata.custom = [True, '; HX', '--', None]
	Run52b2nd	= KAM.sweep(); Run52b2nd.load_hf5('/Run52b/T201503112352', filename = database_location);  Run52b2nd.metadata.custom = [False,'; HX; $f_1$', ':', None]
	Run51aP	    = KAM.sweep();   Run51aP.load_hf5('/Run51a/T201411162028', filename = database_location); 	Run51aP.metadata.custom = [True,'', '-', None]
	Run51bP     = KAM.sweep();   Run51bP.load_hf5('/Run51b/T201411180125', filename = database_location);  Run51bP.metadata.custom = [True,'', '-', None]
	Run49aP     = KAM.sweep();   Run49aP.load_hf5('/Run49a/T201404232144', filename = database_location); Run49aP.metadata.custom = [True,'', '-', None]
	Run48bP     = KAM.sweep();   Run48bP.load_hf5('/Run48b/T201404011305', filename = database_location);  Run48bP.metadata.custom = [True,'', '-', None]
	Run48aP     = KAM.sweep(); Run48aP.load_hf5('/Run48a/T201403311119', filename = database_location);  Run48aP.metadata.custom = [True,'', '-', None]
	Run46aPl    = KAM.sweep(); Run46aPl.load_hf5('/Run46a/T201402211440', filename = database_location);  Run46aPl.metadata.custom = [True,'', '-', None]
	Run46aPh    = KAM.sweep(); Run46aPh.load_hf5('/Run46a/T201402211413', filename = database_location);  Run46aPh.metadata.custom = [True,'', '-',"Run46aPl.Sweep_Array['Fr'][0]"]# high power
	Run45bP     = KAM.sweep(); Run45bP.load_hf5('/Run45b/T201402062055', filename = database_location);  Run45bP.metadata.custom = [True,'; TQ', '-',None]# has loose torque on sensor; to high power
	Run45aP     = KAM.sweep(); Run45aP.load_hf5('/Run45a/T201402061732', filename = database_location);  Run45aP.metadata.custom = [True,'', '-',None]# to high power
	Run44bP     = KAM.sweep(); Run44bP.load_hf5('/Run44b/T201312102229', filename = database_location); Run44bP.metadata.custom = [True,'; CuBx', '-',None]
	Run44aP     = KAM.sweep(); Run44aP.load_hf5('/Run44a/T201312101928', filename = database_location);  Run44aP.metadata.custom = [True,'', '-', None]
	

	#Run49aP.metadata.Sensor = 'S2 FHN2 (4um width)'
	#Run51bP.metadata.Sensor = 'S12 FHN1 (16um width)'
	#Run51aP.metadata.Sensor = 'S5 FHN1 (32um width)'

	Power_Sweeps = [Run51aP,Run51bP,Run49aP,Run48bP,Run48aP,Run46aPl,Run46aPh,Run45bP,Run45aP,Run44bP,Run44aP,Run52b1st,Run52b2nd]

	#needs to be incorporated in to  load_scandata or fill_sweep_array
	# for sweep in Power_Sweeps:
	# 	if sweep.metadata.Atten_Added_At_NA != 0.0:
	# 		sweep.Sweep_Array['Preadout_dB']  = sweep.Sweep_Array['Preadout_dB']-sweep.metadata.Atten_Added_At_NA

if Ver == 2: #temperature sweeps
	Run52bTP2ndh = KAM.sweep(); Run52bTP2ndh.load_hf5('/Run52b/T201503041901', filename = database_location); Run52bTP2ndh.metadata.custom = [False,'; HX; $f_1$', ':', None]
	Run52bTP2ndl = KAM.sweep(); Run52bTP2ndl.load_hf5('/Run52b/T201503021724', filename = database_location); Run52bTP2ndl.metadata.custom = [False,'; HX; $f_1$', ':', None]
	Run52bTP1st =  KAM.sweep();  Run52bTP1st.load_hf5('/Run52b/T201503121941', filename = database_location); Run52bTP1st.metadata.custom  = [True ,'; HX'       ,'--', None]
	Run51aTP = KAM.sweep();         Run51aTP.load_hf5('/Run51a/T201411241956', filename = database_location); Run51aTP.metadata.custom     = [True ,''           , '-', None]
	Run51bTP = KAM.sweep();         Run51bTP.load_hf5('/Run51b/T201411212056', filename = database_location); Run51bTP.metadata.custom     = [True ,''           , '-', None]
	Run49aTP = KAM.sweep();         Run49aTP.load_hf5('/Run49a/T201404181957', filename = database_location); Run49aTP.metadata.custom     = [True ,''           , '-', None]
	Run48bTP = KAM.sweep(); Run48bTP.load_hf5('/Run48b/T201403312220', filename = database_location); Run48bTP.metadata.custom             = [True ,''           , '-', None]
	Run48aTP = KAM.sweep(); Run48aTP.load_hf5('/Run48a/T201403292355', filename = database_location); Run48aTP.metadata.custom             = [True ,''           , '-', None]
	Run46aTP = KAM.sweep(); Run46aTP.load_hf5('/Run46a/T201402262113', filename = database_location); Run46aTP.metadata.custom = [True,'', '-',"Run46aPl.Sweep_Array['Fr'][0]"]# high power
	Run45aTP = KAM.sweep(); Run45aTP.load_hf5('/Run45a/T201402052038', filename = database_location); Run45aTP.metadata.custom = [True,'', '-',None]# to high power #has invalid data
	Run44bTP = KAM.sweep(); Run44bTP.load_hf5('/Run44b/T201312182102', filename = database_location); Run44bTP.metadata.custom = [True,'; CuBx', '-',None] 
	Run44aTP = KAM.sweep(); Run44aTP.load_hf5('/Run44a/T201312191950', filename = database_location); Run44aTP.metadata.custom = [True,'', '-', None]
	TP_Sweeps = [Run51aTP,Run51bTP,Run48bTP,Run49aTP,Run48aTP,Run46aTP,Run45aTP,Run44bTP,Run44aTP,Run52bTP2ndh,Run52bTP2ndl,Run52bTP1st]

	#Run49aTP.metadata.Sensor = 'S2 FHN2 (4um width)'
	#Run51bTP.metadata.Sensor = 'S12 FHN1 (16um width)'
	#Run51aTP.metadata.Sensor = 'S5 FHN1 (32um width)'	

	TP_Sweeps = [Run51aTP,Run51bTP,Run48bTP,Run49aTP,Run48aTP,Run46aTP,Run45aTP,Run44bTP,Run44aTP,Run52bTP2ndh,Run52bTP2ndl,Run52bTP1st]

