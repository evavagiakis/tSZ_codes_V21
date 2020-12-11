import numpy as np
import matplotlib.pyplot as plt

ln_tau=-6.4
m=0.49
err_lt=0.04
err_m=0.08

#ln_tau=-8.06
#m=0.29
#err_lt=0.15
#err_m=0.93

pars=[ln_tau,m,err_lt,err_m]

def tau_scaling_relation(y,ln_t,m):
	ln_t=ln_tau +m*np.log(y/1e-5)
	return ln_t

def mc_tau_scaling_relation_error(y,pars,n):
	ln_tau,m,err_lt,err_m=pars
	mc_ln_tau = ln_tau + np.random.normal(0,err_lt,n)
	mc_m = np.random.normal(m,err_m*m,n)
	mc_taus = tau_scaling_relation(y,mc_ln_tau,mc_m)
	ans=np.std(np.exp(mc_taus))/1e-4
	return ans

ybar=[3.26,2.17 ,1.32 ,0.92 ,0.67 ,1.11 ,0.61 ,0.57 ,0.28] #S18 f150
#ybar=[3.04,2.02 ,1.27 ,0.82 ,0.63 ,1.01 ,0.64 ,0.41 ,0.32] #S18 f090
ybar=[x*1e-7 for x in ybar]
#ybar=[3.92E-07,2.52E-07 ,1.57E-07 ,9.45E-08 ,6.30E-08 ,1.08E-07 ,7.55E-08 ,3.33E-08 ,1.03E-08] #S16 ILC
n=1000

for i in range(len(ybar)):
	print mc_tau_scaling_relation_error(ybar[i],pars,10000)



