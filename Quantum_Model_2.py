import numpy as np
from matplotlib import pylab, mlab, pyplot
plt = pyplot


MHz = np.power(10,6)


# w_0 = 2*np.pi*f_0
# psi_1 = 0
# phi_r = np.pi / 2
# K =  -1.0 * np.power(10,-4.0)*w_0
# gamma_nl = 0.01*K/np.sqrt(3)
# gamma_r = 0.01 * w_0
# gamma_l = 1.1 * gamma_r 
# sigmap_12 = 0.5
# sigmap_13 = np.power(0.5,0.5)
# b1c_in = np.sqrt((4/(3*np.sqrt(3))) * np.power(sigmap_13*gamma_r,-1) * np.power((gamma_r+gamma_l)/(np.abs(K)-np.sqrt(3)*gamma_nl), 3) * (K*K + gamma_nl*gamma_nl)) #Threshold of bi-stability
# b1_in = .8* b1c_in
#computed values:
# psi_2, V3, phi_B, b2_out 

Manual = False
Use_Run_45a = 1

if Use_Run_45a:
	#fit, fig, ax = Run45aP.nonlinear_fit(Save_Fig = True, Indexing = (None,-1,None))
	bestfit = 'Powell'
	eta = fit[bestfit].x[4] #
	delta =  fit[bestfit].x[5]#
	f_0 = fit[bestfit].x[0] #
	Qtl = fit[bestfit].x[1]#
	Qc = fit[bestfit].x[2]#
	phi31 = fit[bestfit].x[3]
	Z1 = Run45aP.metadata.Feedline_Impedance
	Z3 = Run45aP.metadata.Resonator_Impedance
	V30 = np.sqrt(fit['V30V30'])
	Pprobe_dBm = -54
	 #dBm
	phiV1 = 0

if Manual:
	f_0 = 700*MHz
	delta =0.009 #freq nonlinearity
	eta = 10#5.9 #Q nonlinearity
	V30 = 0.1
	Qtl = 100e3
	Qc = 30e3
	Z3  = 50.0
	Z1  = 50.0
	Pprobe_dBm = -65 
	phi31  = np.pi/2.05 #np.pi/2.015
	phiV1 = 1*np.pi/10

Z2 = Z1
Pprobe = 0.001* np.power(10.0,Pprobe_dBm/10.0)
V1V1 = Pprobe *2*Z1
V1 = np.sqrt(V1V1) * np.exp(np.complex(0,1)*phiV1)
V30V30 = np.square(np.abs(V30))
Q = 1.0/ ((1.0/Qtl) + (1.0/Qc))
#tau = 50e-9 #cable delay - only affect V3, but not V2

################################# Create f array making sure it contains f_0
numBW  = 20
BW = numBW*f_0/Q # 2*(f_0 * 0.25) 
num = 3000

if 1: #Triangular numbers
	T = np.linspace(1, num,  num=num, endpoint=True, retstep=False, dtype=None)
	T = T*(T+1.0)/2.0
	f_plus = (T*(BW/2)/T[-1]) + f_0
	f_minus = (-T[::-1]/T[-1])*(BW/2) + f_0
	f = np.hstack((f_minus,f_0,f_plus))

if 0: #linear
	f_plus = np.linspace(f_0, f_0 + BW/2,  num=num, endpoint=True, retstep=False, dtype=None)
	f_minus = np.linspace(f_0 - BW/2,f_0,   num=num-1, endpoint=False, retstep=False, dtype=None)
	f = np.hstack((f_minus,f_plus))


if 0: #logerithmic
	f_plus = np.logspace(np.log10(f_0), np.log10(f_0 + BW/2),  num=num, endpoint=True, dtype=None)
	f_minus = -f_plus[:0:-1] + 2*f_0
	f = np.hstack((f_minus,f_plus))
################################# 

Number_of_Roots = 3

V3V3 = np.ma.empty((f.shape[0],Number_of_Roots), dtype = np.complex128)
V3 = np.ma.empty_like(V3V3)
exp_phi_V3 = np.ma.empty_like(V3V3)
V2_out = np.ma.empty_like(V3V3)

V3V3_up = np.empty_like(f)
V3_up = np.empty_like(f)
V2_out_up = np.empty(f.shape, dtype = np.complex128)

V3V3_down = np.empty(f.shape)
V3_down = np.empty_like(f)
V2_out_down = np.empty_like(f,dtype = np.complex128)
################ 3rd Version

for n in xrange(f.shape[0]):
	coefs = np.array([np.square(delta * f[n]/V30V30 )+ np.square(eta*f_0/(2*Qtl*V30V30)), 2*(delta*(f[n]-f_0)*f[n]/V30V30 + eta*f_0*f_0/(4*Qtl*V30V30*Q) ),np.square(f[n]-f_0) +  np.square(f_0/(2*Q)),  -1.0*f_0*f_0*Z3*V1V1/(4*np.pi*Qc*Z1)])
	V3V3[n] =np.ma.array(np.roots(coefs),mask= np.iscomplex(np.roots(coefs)),fill_value = 1)
	V3[n] =  np.ma.sqrt(V3V3[n])
	
	# exp_phi_V3[n] is e^{i phi_V3} - no minus sign
	exp_phi_V3[n] =  f_0*np.exp(np.complex(0,1.0)*phi31)*V1*np.sqrt(Z3)/(2*np.sqrt(np.pi*Qc*Z1)) *np.power(  ( ((f_0/(2*Q)) +  np.complex(0,1)*(f[n]-f_0)) * V3[n]) + ((eta*f_0/(2*Qtl*V30V30))  +  np.complex(0,1)*(delta * f[n]/V30V30))* V3V3[n]*V3[n],-1.0   )
	
	# V2_out[n] is V_2^out * e^(-i phi_2)
	V2_out[n] = V1*((1-np.exp(np.complex(0,2.0)*phi31))/2 +( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*V3V3[n]/V30V30) + np.complex(0,2)* (((f[n]-f_0)/f_0) + delta*(V3V3[n]/V30V30)*(f[n]/f_0))))*np.exp(np.complex(0,2.0)*phi31))

	# calculate observed for upsweep and down sweep
	# min |--> up sweep (like at UCB),  
	# max |--> down sweep
 	V3V3_down[n]  = np.extract(~V3V3[n].mask,V3V3[n]).max().real
 	V3_down[n]    = np.sqrt(V3V3_down[n])
	V2_out_down[n] = V1*((1-np.exp(np.complex(0,2.0)*phi31))/2 +( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*V3V3_down[n]/V30V30) + np.complex(0,2)* (((f[n]-f_0)/f_0) + delta*(V3V3_down[n]/V30V30)*(f[n]/f_0))))*np.exp(np.complex(0,2.0)*phi31))

	V3V3_up[n] = np.extract(~V3V3[n].mask,V3V3[n]).min().real
	V3_up[n] = np.sqrt(V3V3_up[n])
	V2_out_up[n] = V1*((1-np.exp(np.complex(0,2.0)*phi31))/2 +( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*V3V3_up[n]/V30V30) + np.complex(0,2)* (((f[n]-f_0)/f_0) + delta*(V3V3_up[n]/V30V30)*(f[n]/f_0))))*np.exp(np.complex(0,2.0)*phi31))
####################

# sV3V3_up = ((2.0*V1/Qc) * np.power(2*V2_out_up*np.exp(np.complex(0,-2.0)*phi31)+V1*(1-np.exp(np.complex(0,-2.0)*phi31)),-1) - (1.0/Qc) - (1.0/Qtl) - (np.complex(0,2.0) * (f-f_0)/f_0)) * V30V30  *np.power((eta/Qtl) + np.complex(0,2.0)*(delta*f/f_0),-1)
# sV3V3_up = sV3V3_up.real
# sV3_up =  np.sqrt(sV3V3_up)
# residual1 = np.empty(sV3V3_up.shape,dtype = np.complex128)
# residual2 = np.empty(sV3V3_up.shape,dtype = np.complex128)
# residual3 = np.empty(sV3V3_up.shape,dtype = np.complex128)
# residual4 = np.empty(sV3V3_up.shape,dtype = np.complex128)

# for n in xrange(f.shape[0]):
# 	coefs = np.array([np.square(delta * f[n]/V30V30 )+ np.square(eta*f_0/(2*Qtl*V30V30)), 2*(delta*(f[n]-f_0)*f[n]/V30V30 + eta*f_0*f_0/(4*Qtl*V30V30*Q) ),np.square(f[n]-f_0) +  np.square(f_0/(2*Q)),  -1.0*f_0*f_0*Z3*V1V1/(4*np.pi*Qc*Z1)])
# 	residual1[n] = f_0*np.exp(np.complex(0,1)*phi31)*V1*np.sqrt(Z3)/(2*np.sqrt(np.pi*Qc*Z1))    -    (( ((f_0/(2.0*Q)) +  np.complex(0,1.0)*(f[n]-f_0)) * sV3_up[n]) + ((eta*f_0/(2.0*Qtl*V30V30))  +  np.complex(0,1)*(delta * f[n]/V30V30))* sV3V3_up[n]*sV3_up[n])
# 	residual2[n] = f_0*np.exp(np.complex(0,1)*phi31)*V1*np.sqrt(Z3)/(2*np.sqrt(np.pi*Qc*Z1))    -    (( ((f_0/(2.0*Q)) +  np.complex(0,1.0)*(f[n]-f_0)) *  V3_up[n]) + ((eta*f_0/(2.0*Qtl*V30V30))  +  np.complex(0,1)*(delta * f[n]/V30V30))*  V3V3_up[n]* V3_up[n])
# 	residual3[n] = np.polyval(coefs,sV3V3_up[n] ) # Exaluate the V3V3 qubic using the  sV3V3_up synthesized from S21 
# 	residual4[n] = np.polyval(coefs,V3V3_up[n] ) # Exaluate the V3V3 qubic using the  V3V3_up computed from  polynomial roots
# 	#if residual2 - residual3 = 0 then V3V3_up = sV3V3_up to high enough accuracy
# sumsq = np.square(residual2).sum()

#We use the solution to the cubic for for one scan direction to construct the other two solutions
V2cubic = V2_out_down
S21 = V2cubic/V1
V3__ = np.empty_like(f)
V3__  = (S21 + (np.exp(np.complex(0,2.0)*phi31)-1)/2.)*V1*np.sqrt(Z3*Qc/(Z1*np.pi))* np.exp(np.complex(0,-1.0)*phi31)
z1 = eta/(Qtl*V30V30)+ np.complex(0,1.0)*(2*delta*f)/(V30V30*f_0)
z2 = (1.0/Qc) + (1/Qtl) + np.complex(0,2.0) *(f-f_0)/f_0
z1z2c = z1*z2.conjugate()
z1z1 = z1*z1.conjugate()
z2z2 = z2*z2.conjugate()
v1 = V3__*V3__.conjugate()
term1 = -(z1z2c.real/z1z1) - v1/2.0
term2 = np.complex(0,1)*np.sqrt(4*z1z2c.imag*z1z2c.imag + 3*v1*v1*z1z1*z1z1 + 4*z1z1*z1z2c.real*v1)/(2*z1z1)
v2  = term1 + term2
v3  = term1 - term2

V3p__ = np.sqrt(v2)
V3m__ = np.sqrt(v3)

S21p= ((1-np.exp(np.complex(0,2.0)*phi31))/2 +( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*v2/V30V30) + np.complex(0,2)* (((f-f_0)/f_0) + delta*(v2/V30V30)*(f/f_0))))*np.exp(np.complex(0,2.0)*phi31))
S21m = ((1-np.exp(np.complex(0,2.0)*phi31))/2 +( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*v3/V30V30) + np.complex(0,2)* (((f-f_0)/f_0) + delta*(v3/V30V30)*(f/f_0))))*np.exp(np.complex(0,2.0)*phi31))



#V3c__ = V3__.conjugate()
#f_0*np.exp(np.complex(0,1)*phi31)*V1*np.sqrt(Z3)/(2*np.sqrt(np.pi*Qc*Z1))    -    (( ((f_0/(2.0*Q)) +  np.complex(0,1.0)*(f[n]-f_0)) * sV3_up[n]) + ((eta*f_0/(2.0*Qtl*V30V30))  +  np.complex(0,1)*(delta * f[n]/V30V30))* sV3V3_up[n]*sV3_up[n])







############### 2nd Version
# def roots(freq):
# 	coefs = [np.square(delta * freq/V30V30 )+ np.square(eta*f_0/(4*Qtl*V30V30)), 2*(delta*(freq-f_0)*freq/V30V30 + eta*f_0*f_0/(4*Qtl*V30V30*Q) ),np.square(freq-f_0) +  np.square(f_0/(2*Q)),  -1.0*f_0*f_0*Z3*V1V1/(4*np.pi*Qc*Z1)]
# 	return np.roots(coefs)

# for n in xrange(f.shape[0]):
# 	V3V3[n] = roots(f[n])
# 	V3[n] =  np.sqrt(V3V3[n])
# 	test[n] = V3[n]
# 	# exp_phi_V3[n] is e^{i phi_V3} - no minus sign
# 	exp_phi_V3[n] =  f_0*np.exp(np.complex(0,1)*phi31)*V1*np.sqrt(Z3)/(2*np.sqrt(np.pi*Qc*Z1)) / (  ((f_0/(2*Q)) +  np.complex(0,1)*(f[n]-f_0) *V3[n]) + ((eta*f_0/(2*Qtl*V30V30))  +  np.complex(0,1)*(delta * f[n]/V30V30))* V3V3[n]*V3[n]   )
# 	# V2_out[n] is V_2^out * e^(-i phi_2)
# 	V2_out[n] = V1*(1 -( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*V3V3[n]/V30V30) + 2*np.complex(0,1)* (((f[n]-f_0)/f_0) + delta*(V3V3[n]/V30V30)*(f[n]/f_0)))))

# V3V3 = np.ma.masked_where(np.iscomplex(V3V3), V3V3, copy=True)
# V3V3.fill_value = 1
# V3.mask = V3V3.mask
# exp_phi_V3.mask = V3V3.mask
# V2_out.mask = V3V3.mask

# for n in xrange(f.shape[0]):
# 	#calculate observed upsweep values
# 	V3V3_up[n] = np.max(np.abs(V3V3[n].compressed()))
# 	V3_up[n] = np.sqrt(V3V3_up[n])
# 	V2_out_up[n] = V1*(1 -( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*V3V3_up[n]/V30V30) + 2*np.complex(0,1)* (((f[n]-f_0)/f_0) + delta*(V3V3_up[n]/V30V30)*(f[n]/f_0)))))
#################

# ################## 1st Version
# for n in xrange(f.shape[0]):
# 	V3V3[n] = roots(f[n])
	
# 	#Where are there 3 real solutions?
# 	# if  np.isreal(V3V3[n,0]) & np.isreal(V3V3[n,1]) & np.isreal(V3V3[n,2]): 
# 	# 	print(n)
	


# V3V3 = np.ma.masked_where(np.iscomplex(V3V3), V3V3, copy=True)
# V3 = np.sqrt(V3V3)

# exp_phi_V3 = np.ma.empty_like(V3)
# V2_out = np.ma.empty_like(V3)
# for n in xrange(f.shape[0]):

	
	
# 	# exp_phi_V3[n] is e^{i phi_V3} - no minus sign
# 	exp_phi_V3[n] =  f_0*np.exp(np.complex(0,1)*phi31)*V1*np.sqrt(Z3)/(2*np.sqrt(np.pi*Qc*Z1)) / (  ((f_0/(2*Q)) +  np.complex(0,1)*(f[n]-f_0) *V3[n]) + ((eta*f_0/(2*Qtl*V30V30))  +  np.complex(0,1)*(delta * f[n]/V30V30))* V3V3[n]*V3[n]   )
# 	# V2_out_phasor[n] is V_2^out * e^(-i phi_2)
# 	V2_out[n] = V1*(1 -( (1/Qc) / ((1/Qc) + (1/Qtl)*(1+eta*V3V3[n]/V30V30) + 2*np.complex(0,1)* (((f[n]-f_0)/f_0) + delta*(V3V3[n]/V30V30)*(f[n]/f_0)))))
# ##################


fig = plt.figure( figsize=(6, 6), dpi=150)
ax = {}
ax[1] = fig.add_subplot(2,2,1)
dff = (f - f_0)/f_0
trans = (V2_out/V1)
# dfff = np.array([dff,dff,dff]).transpose()
# dff = ma.array(dfff, mask = trans.mask)
# trans2 = trans.compressed()
# dff2 = dff.compressed()
trans_up = V2_out_up/V1 
trans_down = (V2_out_down/V1)
transp=S21p[~np.iscomplex(V3p__)]
transm=S21m[~np.iscomplex(V3m__)]

curve = ax[1].plot(dff,20*np.log10(np.abs(trans)),color = 'g', linestyle = '-',linewidth = 2)
curve_up = ax[1].plot(dff,20*np.log10(np.abs(trans_up)), color = 'k', linestyle = ':', alpha  = .35,linewidth = 1, label = 'Up Sweep')
curve_down = ax[1].plot(dff,20*np.log10(np.abs(trans_down)), color = 'k', linestyle = '--', alpha  = .35, linewidth = 1,label = 'Down Sweep')
ax[1].set_title('Mag Transmission')
ax[1].set_xlabel(r'$\delta f_0 / f_0$', color='k')
ax[1].set_ylabel(r'$20 \cdot \log_{10}|S_{21}|$ [dB]', color='k')
ax[1].yaxis.labelpad = 0
ax[1].ticklabel_format(axis='x', style='sci',scilimits = (0,0), useOffset=True)
#ax[1].xaxis.set_ticks(np.hstack((np.arange(-numBW/2.0,0,f_0/Q),np.arange(0,numBW/2.0,f_0/Q))) )
parameter_dict = {'f_0':f_0, 'Qtl':Qtl, 'Qc':Qc, 'phi31':phi31, 'eta':eta, 'delta':delta, 'Zfl':Z1, 'Zres':Z3,  'phiV1':phiV1, 'V30V30':V30*V30}
note = '$P_{probe}$' + ' {:3.0f} dBm, '.format(Pprobe_dBm)+'\n' +(r'$f_0$ = {f_0:3.2e} Hz,' + '\n' + '$Q_{sub1}$ = {Qtl:3.2e},' +'\n' +' $Q_c$ = {Qc:3.2e},' +
	'\n' + r'$\phi_{sub2}$ = {ang:3.2f}$^\circ$,'+ '\n' + '${l1}$ = {et:3.2e},' + '\n' +'${l2}$ = {de:3.2e}').format(
	nl = '\n', et = parameter_dict['eta']/parameter_dict['V30V30'],
	de = parameter_dict['delta']/parameter_dict['V30V30'], 
	l1 = r'{\eta}/{V_{3,0}^2}',
	l2  = r'{\delta}/{V_{3,0}^2}',
	ang = parameter_dict['phi31']*180/np.pi, 
	sub1 = '{i}', sub2 = '{31}',**parameter_dict)
ax[1].text(0.99, 0.01, note,
	verticalalignment='bottom', horizontalalignment='right',
	transform=ax[1].transAxes,
	color='black', fontsize=4)

ax[2] = fig.add_subplot(2,2,2)
curve = ax[2].plot(dff,np.abs(V3),color = 'g', linestyle = '-',linewidth = 2)# <- V3 has complex values when it shouldn't !! should this be real part or abs?
upcurve = ax[2].plot(dff,np.abs(V3_up),color = 'k', linestyle = ':', alpha  = .35,linewidth = 1, label = 'Up Sweep')
upcurve = ax[2].plot(dff,np.abs(V3_down),color = 'k', linestyle = '--', alpha  = .35, linewidth = 1,label = 'Down Sweep')
#upcurve__ = ax[2].plot(dff,np.abs(V3__),linestyle = '--')
#curve__ = ax[2].plot(dff[~np.iscomplex(V3p__)].real,V3p__[~np.iscomplex(V3p__)].real,linestyle = '--', marker = '1')
#curve__ = ax[2].plot(dff[~np.iscomplex(V3m__)].real,V3m__[~np.iscomplex(V3m__)].real,linestyle = '--', marker = '2')
ax[2].set_title('Cavity Amplitude')
ax[2].set_xlabel(r'$\delta f_0 / f_0$', color='k')
ax[2].set_ylabel(r'Volts', color='k')
ax[2].ticklabel_format(axis='x', style='sci',scilimits = (0,0),useOffset=False)



ax[3] = fig.add_subplot(2,2,3,aspect='equal')
loop  = ax[3].plot(trans.real, trans.imag,color = 'g', linestyle = '-',linewidth = 2)#, label = 'Full Solution')
loop[0].set_label('Full Solution')
loop_down  = ax[3].plot(trans_down.real, trans_down.imag,color = 'k', linestyle = '--', alpha  = .35, linewidth = 1,label = 'Down Sweep')
#loop = ax[3].plot(transp.real,transp.imag,linestyle = '--', marker = '1')
#oop = ax[3].plot(transm.real,transm.imag,linestyle = '--', marker = '2')
#firstpt  = ax[3].plot(trans.real[0:10], trans.imag[0:10], 'ok')
loop_up  = ax[3].plot(trans_up.real, trans_up.imag,color = 'k', linestyle = ':', alpha  = .35,linewidth = 1, label = 'Up Sweep')
ax[3].set_title('Resonance Loop')
ax[3].set_xlabel(r'$\Re$[$S_{21}$]', color='k')
ax[3].set_ylabel(r'$\Im$[$S_{21}$]', color='k')
ax[3].yaxis.labelpad = 0
ax[3].ticklabel_format(axis='x', style='sci',scilimits = (0,0),useOffset=False)
ax[3].legend(loc = 'upper center', fontsize=7, bbox_to_anchor=(1.5, -.15),  ncol=3,scatterpoints =1, numpoints = 1, labelspacing = .02)

ax[4] = fig.add_subplot(2,2,4)
trans_phase = np.ma.array(np.angle(trans),mask = trans.mask)
trans_up_phase = np.angle(trans_up)
trans_down_phase = np.angle(trans_down)
phase_ang_curve = ax[4].plot(dff,trans_phase,color = 'g', linestyle = '-',linewidth = 2)
phase_up_ang_curve = ax[4].plot(dff,trans_up_phase,color = 'k', linestyle = ':', alpha  = .35,linewidth = 1, label = 'Up Sweep')
phase_down_ang_curve = ax[4].plot(dff,trans_down_phase,color = 'k', linestyle = '--', alpha  = .35,linewidth = 1, label = 'Down Sweep')
ax[4].set_title('Transmitted Phase Angle')
ax[4].set_xlabel(r'$\delta f_0 / f_0$', color='k')
ax[4].set_ylabel(r'Ang[$S_{21}$]', color='k')
ax[4].yaxis.labelpad = 0
ax[4].ticklabel_format(axis='x', style='sci',scilimits = (0,0),useOffset=False)



for k in ax.keys():
	ax[k].tick_params(axis='y', labelsize=5)
	ax[k].tick_params(axis='x', labelsize=5)

plt.subplots_adjust(left=.1, bottom=.1, right=None ,wspace=.35, hspace=.3)
#plt.subplots_adjust(left=.1, bottom=.1, right=None, top=.95 ,wspace=.4, hspace=.4)
#plt.suptitle('Nonlinear Resonator Plots')
plt.show()
if Use_Run_45a:
	Title = '45a_Nonlinear_Solition_Pprobe_-54dBm'
	#swp._save_fig_dec(fig, Title.replace('\n','_').replace(' ','_'), Use_Date = Use_Date )

#fig.savefig('Nonlinear_Res',dpi=300, transparency  = True)

