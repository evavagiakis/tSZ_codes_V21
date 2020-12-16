import numpy as np
import matplotlib.pyplot as plt

ln_tau=-6.4 #parameters and errors for hydrosim relationship 
m=0.49
err_lt=0.04
err_m=0.08

pars=[ln_tau,m,err_lt,err_m] #put into array

def tau_calc(y):
	tau = np.e**(ln_tau+m*np.log(y/10**(-5))) #Nick's hydrosim relationship
	return tau 

def tau_err_calc(y, yerr):
	tauerrval = 0.229462  #(m*e**lntau0)/(1e-5)**m, factor for statistical error propagation
	tauerr = tauerrval*yerr/(y**(1.0-m)) #propagate yerr
	return tauerr

def tau_scaling_relation(y,ln_t,m):
        ln_t=ln_tau +m*np.log(y/1e-5) #Nick's hydrosim relationship
        return ln_t

def mc_tau_scaling_relation_error(y,pars,n): #Nick's hydrosim relationship MCMC sampler for sys err estimate
        ln_tau,m,err_lt,err_m=pars
        mc_ln_tau = ln_tau + np.random.normal(0,err_lt,n)
        mc_m = np.random.normal(m,err_m*m,n)
        mc_taus = tau_scaling_relation(y,mc_ln_tau,mc_m)
        ans=np.std(np.exp(mc_taus))
        return ans

#Bins are in the following order: L116, L98, L79, L61, L43, L98D, L79D, L61D, L43D

f150d=np.loadtxt('S18f150Donutraw.txt') #Load raw aperture photometry results for three map analyses 
f090d=np.loadtxt('S18f090raw.txt')
ILCd= np.loadtxt('S16ILCraw.txt')

f150=f150d[0,:] #extract dTs 
f150err=f150d[1,:] #extract errs
f090=f090d[0,:] #extract dTs
f090err=f090d[1,:] #extract errs
ILC=ILCd[0,:] #extract ys 
ILCerr=ILCd[1,:] #extract errs

################################ Dust corrections ###########################################################
f150Dust=[0.043735,0.036873 ,0.028311 ,0.027812 ,0.025492 ,0.045007 ,0.028342 ,0.036278 ,0.026553] #Herschel dust corrections from Stefania
f150Dusterr=[0.059536,0.028298 ,0.01674 ,0.017594 ,0.013509 ,0.031073 ,0.028272 ,0.039795 ,0.02056] #Herschal dust correction upper error bar

f090Dust=[0.018791,0.015356 ,0.011432 ,0.012232 ,0.010424 ,0.019142 ,0.011636 ,0.015865 ,0.01118] #Herschel dust corrections from Stefania
f090Dusterr=[0.03034,0.012307 ,0.007802 ,0.008122 ,0.006161 ,0.01478 ,0.013613 ,0.018879 ,0.010186] #Herschal dust correction upper error bar

f150=[f150[x]-f150Dust[x] for x in range(len(f150))] #Perform Stefania's Herschel dust correction
f150err=[np.sqrt(f150err[x]**2.+f150Dusterr[x]**2.) for x in range(len(f150))] #Propogate Stefania's Herschel dust correction error

f090=[f090[x]-f090Dust[x] for x in range(len(f090))] #Perform Stefania's Herschel dust correction
f090err=[np.sqrt(f090err[x]**2.+f090Dusterr[x]**2.) for x in range(len(f090))] #Propogate Stefania's Herschel dust correction error
#############################################################################################################

################################ Beam corrections ###########################################################
f090BeamCorr=[0.762,0.762 ,0.763 ,0.764 ,0.764 ,0.763 ,0.764 ,0.766 ,0.766] #From Stefania's beam corrections
ILCBeamCorr=[1.051,1.053 ,1.053 ,1.056 ,1.056 ,1.053 ,1.056 ,1.059 ,1.059] #From Stefania's beam corrections 

f090=[f090[x]/f090BeamCorr[x] for x in range(len(f090))] #Perform Stefania's beam correction
f090err=[f090err[x]/f090BeamCorr[x] for x in range(len(f090))] #Perform Stefania's beam correction

ILC=[ILC[x]/ILCBeamCorr[x] for x in range(len(ILC))] #Perform Stefania's beam correction
ILCerr=[ILCerr[x]/ILCBeamCorr[x] for x in range(len(ILC))] #Perform Stefania's beam correction
#############################################################################################################

################################ Conversion to y  ###########################################################
cfac150=0.00000038312 #fSZ for S18 f150
cfac090=0.00000024007 #fSZ for S18 f090

f150=[-x*cfac150 for x in f150] #Converting dT to y
f150err=[x*cfac150 for x in f150err] #Converting dTerr to yerr

f090=[-x*cfac090 for x in f090] #Converting dT to y
f090err=[x*cfac090 for x in f090err] #Converting dTerr to yerr
#############################################################################################################

############################### Calculate tau ###############################################################
tau150=[tau_calc(f150[x]) for x in range(len(f150))] #Tau estimates from hydrosim relationship
tau090=[tau_calc(f090[x]) for x in range(len(f090))]
tauILC=[tau_calc(ILC[x]) for x in range(len(ILC))]

tau150err=[tau_err_calc(f150[x],f150err[x]) for x in range(len(f150err))] #Tau statistical errors from hydrosim relationship
tau090err=[tau_err_calc(f090[x],f090err[x]) for x in range(len(f090err))]
tauILCerr=[tau_err_calc(ILC[x],ILCerr[x]) for x in range(len(ILCerr))]
#############################################################################################################

############################## MCMC Sampler for sys err estimate on tau #####################################
n=1000
tau150sys=[mc_tau_scaling_relation_error(x,pars,10000) for x in f150] #Tau systematic errors from MCMC sampler
tau090sys=[mc_tau_scaling_relation_error(x,pars,10000) for x in f090]
tauILCsys=[mc_tau_scaling_relation_error(x,pars,10000) for x in ILC]
#############################################################################################################

############################## Theory tau comparisons #######################################################
theorytau = [4.439284824479226E-4,3.348820653104621E-4 ,2.4216433534549755E-4 ,1.767778718898415E-4 ,1.3927485788583183E-4 ,2.0889705718306475E-4 ,1.5270659977613805E-4 ,1.058116994456508E-4 ,0.696660489891753E-4]

fc150=[tau150[x]/theorytau[x] for x in range(len(tau150))]
fc090=[tau090[x]/theorytau[x] for x in range(len(tau090))]
fcILC=[tauILC[x]/theorytau[x] for x in range(len(tauILC))]

fc150err=[tau150err[x]/theorytau[x] for x in range(len(tau150))]
fc090err=[tau090err[x]/theorytau[x] for x in range(len(tau090))]
fcILCerr=[tauILCerr[x]/theorytau[x] for x in range(len(tauILC))]

fc150sys=[tau150sys[x]/theorytau[x] for x in range(len(tau150))]
fc090sys=[tau090sys[x]/theorytau[x] for x in range(len(tau090))]
fcILCsys=[tauILCsys[x]/theorytau[x] for x in range(len(tauILC))]

#############################################################################################################

print fc150sys
print fc090sys
print fcILCsys







