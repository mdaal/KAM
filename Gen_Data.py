#Generate nonlinear resonance test data.
import datetime
import numpy as np
from matplotlib import pylab, mlab, pyplot
plt = pyplot

import KAM
reload(KAM)
import sys
MHz = np.power(10,6)


######################### Choose Variable Values
f_0 = 700*MHz
Qtl = 80e3
Qc = 30e3
phi31  = np.pi/2.03 #np.pi/2.015
eta = .01 #Q nonlinearity
delta =0.000008 #freq nonlinearity


Z3  = 30.0
Z1  = 50.0


phiV1 = 1*np.pi/10

V30 = 0.1
V30V30 = np.square(np.abs(V30))

Q = 1.0/ ((1.0/Qtl) + (1.0/Qc))

Pprobe_dBm_Start = -65.0
Pprobe_dBm_Stop = -25.0
Pprobe_Num_Points = 10

numBW  = 24
num = 1000

#Up or Down

#Show_Plot

############################## Make Pprobe Array
Pprobe_dBm = np.linspace(Pprobe_dBm_Start,Pprobe_dBm_Stop, Pprobe_Num_Points)
Pprobe = 0.001* np.power(10.0,Pprobe_dBm/10.0)
V1V1 = Pprobe *2*Z1
V1 = np.sqrt(V1V1) * np.exp(np.complex(0,1)*phiV1)


################################# Create f array making sure it contains f_0

BW = numBW*f_0/Q # 2*(f_0 * 0.25) 


if 0: #Triangular numbers - Denser around f_0
	T = np.linspace(1, num,  num=num, endpoint=True, retstep=False, dtype=None)
	T = T*(T+1.0)/2.0
	f_plus = (T*(BW/2)/T[-1]) + f_0
	f_minus = (-T[::-1]/T[-1])*(BW/2) + f_0
	f = np.hstack((f_minus,f_0,f_plus))

if 1: #linear
	f_plus = np.linspace(f_0, f_0 + BW/2,  num=num, endpoint=True, retstep=False, dtype=None)
	f_minus = np.linspace(f_0 - BW/2,f_0,   num=num-1, endpoint=False, retstep=False, dtype=None)
	f = np.hstack((f_minus,f_plus))


if 0: #logerithmic - Denser around f_0, for wide band sweeps
	f_plus = np.logspace(np.log10(f_0), np.log10(f_0 + BW/2),  num=num, endpoint=True, dtype=None)
	f_minus = -f_plus[:0:-1] + 2*f_0
	f = np.hstack((f_minus,f_plus))
################################# 

#################### Declare Arrays
Number_of_Roots = 3
V3V3 = np.ma.empty((f.shape[0],Number_of_Roots), dtype = np.complex128)

V3V3_up = np.empty_like(f)
V3_up = np.empty_like(f)
V2_out_up = np.empty(f.shape, dtype = np.complex128)
S21_up = np.empty(f.shape, dtype = np.complex128)

V3V3_down = np.empty(f.shape)
V3_down = np.empty_like(f)
V2_out_down = np.empty_like(f,dtype = np.complex128)
S21_down = np.empty_like(f,dtype = np.complex128)
################ 


#################### Creat KAM object and fill it.
swp    = KAM.sweep() # Create KAM object so we can fill it with fake data

tpoints = 0
swp._define_sweep_data_columns(f.size,tpoints)
swp.Sweep_Array = np.zeros(Pprobe_dBm.size, dtype = swp.sweep_data_columns) #Sweep_Array holdes all sweep data. Its length is the number of sweeps

fig = plt.figure( figsize=(5, 5), dpi=150)
ax = {}
ax[1] = fig.add_subplot(1,1,1)
dff = (f - f_0)/f_0



swp.metadata.Electrical_Delay = 0.0
swp.metadata.Through_Line_Impedance = Z1
swp.metadata.Resonator_Impedance = Z3
swp.metadata.Time_Created =     '05/01/2015 12:00:00'
swp.metadata.Run = 'F1'

print 'Generating False Data...'
for index in xrange(Pprobe_dBm.size):
	sys.stdout.write('\r {0} of {1} '.format(index+1, Pprobe_dBm.size))
	sys.stdout.flush()
	for n in xrange(f.shape[0]):
		coefs = np.array([np.square(delta * f[n]/V30V30 )+ np.square(eta*f_0/(2*Qtl*V30V30)), 2*(delta*(f[n]-f_0)*f[n]/V30V30 + eta*f_0*f_0/(4*Qtl*V30V30*Q) ),np.square(f[n]-f_0) +  np.square(f_0/(2*Q)),  -1.0*f_0*f_0*Z3*V1V1[index]/(4*np.pi*Qc*Z1)])
		V3V3[n] =np.ma.array(np.roots(coefs),mask= np.iscomplex(np.roots(coefs)),fill_value = 1)

		# calculate observed for upsweep and down sweep
		# min |--> up sweep (like at UCB),  
		# max |--> down sweep
	 	V3V3_down[n]    = np.extract(~V3V3[n].mask,V3V3[n]).max().real
	 	V3_down[n]      = np.sqrt(V3V3_down[n])
		V2_out_down[n]  = V1[index]*((1-np.exp(np.complex(0,2.0)*phi31))/2 +( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*V3V3_down[n]/V30V30) + np.complex(0,2)* (((f[n]-f_0)/f_0) + delta*(V3V3_down[n]/V30V30)*(f[n]/f_0))))*np.exp(np.complex(0,2.0)*phi31))
		S21_down[n]     = (V2_out_down[n]/V1[index])

		V3V3_up[n]      = np.extract(~V3V3[n].mask,V3V3[n]).min().real
		V3_up[n]        = np.sqrt(V3V3_up[n])
		V2_out_up[n]    = V1[index]*((1-np.exp(np.complex(0,2.0)*phi31))/2 +( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*V3V3_up[n]/V30V30) + np.complex(0,2)* (((f[n]-f_0)/f_0) + delta*(V3V3_up[n]/V30V30)*(f[n]/f_0))))*np.exp(np.complex(0,2.0)*phi31))
		S21_up[n]       = V2_out_up[n]/V1[index] 

		
			
	swp._define_sweep_array(index, Fstart = f[0], #Hz
								Fstop = f[-1], #Hz
								S21 = S21_up,
								Frequencies = f, #Hz
								Preadout_dB = Pprobe_dBm[index],
								Pinput_dB = Pprobe_dBm[index] - 50.0,
								Is_Valid = True,
								Mask = np.zeros(f.shape, dtype=np.bool),
								Fr = f_0, #Note! we are using the resonance freq of the lowest power S21 for all 
								Q = Q,
								Qc = Qc,
								Heater_Voltage = 0.0
								)

	curve_up = ax[1].plot(dff,np.abs(S21_up), linestyle = '-')
	#curve_down = ax[1].plot(dff,np.abs(S21_down), linestyle = '--', color = 'r')



################ Configure Plot
ax[1].set_title('Mag Transmission')
ax[1].set_xlabel(r'$\delta f_0 / f_0$', color='k')
ax[1].set_ylabel(r'Mag[$S_{21}$]', color='k')
ax[1].yaxis.labelpad = -4
ax[1].ticklabel_format(axis='x', style='sci',scilimits = (0,0), useOffset=True)

for k in ax.keys():
	ax[k].tick_params(axis='y', labelsize=9)
	ax[k].tick_params(axis='x', labelsize=5)


#plt.subplots_adjust(left=.1, bottom=.1, right=None, top=.95 ,wspace=.4, hspace=.4)
plt.suptitle('Nonlinear Resonator Plots')
plt.show()
