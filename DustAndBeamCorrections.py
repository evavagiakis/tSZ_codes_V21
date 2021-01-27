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

f150d=np.loadtxt('DR5f150_nocoreexcised_20201223_AP.txt') #Load raw aperture photometry results for three map analyses 
f090d=np.loadtxt('DR5f090_nocoreexcised_20201223_AP.txt')
ILCd= np.loadtxt('ILC_20201223_AP.txt')

f150=f150d[0,:] #extract dTs 
f150err=f150d[1,:] #extract errs
f090=f090d[0,:] #extract dTs
f090err=f090d[1,:] #extract errs
ILC=ILCd[0,:] #extract ys 
ILCerr=ILCd[1,:] #extract errs

print 'raw s/n:', [f150[x]/f150err[x] for x in range(len(f150))]
print 'raw s/n:', [f090[x]/f090err[x] for x in range(len(f090))]
print 'raw s/n:', [ILC[x]/ILCerr[x] for x in range(len(ILC))]


################################ Dust corrections ###########################################################
d150=np.loadtxt('dust_coadd_f150.txt')
d090=np.loadtxt('dust_coadd_f090.txt')

binselect=[4,3,2,1,0,8,7,6,5] #bins are in a different order in the text files than what we are working with here

f150Dust=[d150[x,2] for x in binselect]#Herschel dust corrections from Stefania
f150Dusterr=[d150[x,3] for x in binselect] #Herschal dust correction upper error bar

f090Dust=[d090[x,2] for x in binselect] #Herschel dust corrections from Stefania
f090Dusterr=[d090[x,3] for x in binselect] #Herschal dust correction upper error bar

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

####Print dT results
print 'dT f150 corrected:', [np.round(x,2) for x in f150]
print 'dT f150 corrected err:', [np.round(x,2) for x in f150err]
print 'dT f090 corrected:', [np.round(x,2) for x in f090]
print 'dT f090 corrected err:', [np.round(x,2) for x in f090err]

################################ Conversion to y  ###########################################################
cfac150=0.00000038312 #fSZ for DR5 f150
cfac090=0.00000024007 #fSZ for DR5 f090

f150=[-x*cfac150 for x in f150] #Converting dT to y
f150err=[x*cfac150 for x in f150err] #Converting dTerr to yerr

f090=[-x*cfac090 for x in f090] #Converting dT to y
f090err=[x*cfac090 for x in f090err] #Converting dTerr to yerr
#############################################################################################################

####Print y results
print 'y f150 corrected:', [np.round(x/1e-7,2) for x in f150]
print 'y f150err corrected:', [np.round(x/1e-7,2) for x in f150err]
print 'y f090 corrected:', [np.round(x/1e-7,2) for x in f090]
print 'y f090err corrected:', [np.round(x/1e-7,2) for x in f090err]
print 'y ILC corrected:', [np.round(x/1e-7,2) for x in ILC]
print 'y ILCerr corrected:', [np.round(x/1e-7,2) for x in ILCerr]



############################### Calculate tau ###############################################################
tau150=[tau_calc(f150[x]) for x in range(len(f150))] #Tau estimates from hydrosim relationship
tau090=[tau_calc(f090[x]) for x in range(len(f090))]
tauILC=[tau_calc(ILC[x]) for x in range(len(ILC))]

tau150err=[tau_err_calc(f150[x],f150err[x]) for x in range(len(f150err))] #Tau statistical errors from hydrosim relationship
tau090err=[tau_err_calc(f090[x],f090err[x]) for x in range(len(f090err))]
tauILCerr=[tau_err_calc(ILC[x],ILCerr[x]) for x in range(len(ILCerr))]
#############################################################################################################

#####Print tau results
print 'tau 150', [np.round(x/1e-4,2) for x in tau150]
print 'tau 150err', [np.round(x/1e-4,2) for x in tau150err]
print 'tau 090', [np.round(x/1e-4,2) for x in tau090]
print 'tau 090err', [np.round(x/1e-4,2) for x in tau090err]
print 'tau ILC', [np.round(x/1e-4,2) for x in tauILC]
print 'tau ILCerr', [np.round(x/1e-4,2) for x in tauILCerr]


############################## MCMC Sampler for sys err estimate on tau #####################################
n=1000
tau150sys=[mc_tau_scaling_relation_error(x,pars,10000) for x in f150] #Tau systematic errors from MCMC sampler
tau090sys=[mc_tau_scaling_relation_error(x,pars,10000) for x in f090]
tauILCsys=[mc_tau_scaling_relation_error(x,pars,10000) for x in ILC]
#############################################################################################################

########Print tau sims
print 'tau 150 sys err', [np.round(x/1e-4,2) for x in tau150sys]
print 'tau 90 sys err', [np.round(x/1e-4,2) for x in tau090sys]
print 'tau ILC sys err', [np.round(x/1e-4,2) for x in tauILCsys]

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

####Print fc results
print 'fc150',[np.round(x,2) for x in fc150]
print 'fc090',[np.round(x,2) for x in fc090]
print 'fcILC',[np.round(x,2) for x in fcILC]

print 'fc150err',[np.round(x,2) for x in fc150err]
print 'fc090err',[np.round(x,2) for x in fc090err]
print 'fcILCerr',[np.round(x,2) for x in fcILCerr]

print 'fc150sys',[np.round(x,2) for x in fc150sys]
print 'fc090sys',[np.round(x,2) for x in fc090sys]
print 'fcILCsys',[np.round(x,2) for x in fcILCsys]





