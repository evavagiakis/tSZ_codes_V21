import argparse
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
import numpy
import pylab
from pylab import *
from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel
from scipy import signal
from scipy.integrate import quad
from scipy.integrate import dblquad
import pandas as pd
from scipy import interpolate

dts18150=(-0.84994731,-0.5656846779 ,-0.3441183991 ,-0.2405159473 ,-0.1761323741 ,-0.2888925015 ,-0.1591828689 ,-0.1477630680 ,-0.0731822341)
dts18150err=(0.08883084,0.0539037406 ,0.0338921927 ,0.0275226831 ,0.0215700681 ,0.0658364357 ,0.0471320325 ,0.0495105052 ,0.0330823968)

dts1890=(-1.265355965,-0.840460965 ,-0.527565256 ,-0.343358492 ,-0.262665673 ,-0.421167194 ,-0.268306625 ,-0.171861893 ,-0.134699160)
dts1890err=(0.093433783,0.060832087 ,0.038352344 ,0.027698164 ,0.021982106 ,0.079441095 ,0.050980324 ,0.041571323 ,0.037002407)

dtILC=(-1.06316309982,-0.68403463443 ,-0.42652672757 ,-0.25640983359 ,-0.17091666932 ,-0.29220881417 ,-0.20494594345 ,-0.09025365414 ,-0.02805127692)
dtILCerr=(0.09027830721,0.06082676724 ,0.04071116034 ,0.02877718669 ,0.02343171171 ,0.07723131179 ,0.05010420328 ,0.03558547460 ,0.03530605876)

maxt = 1.8
if maxt == 1.3:
     lntau0 = 6.34
     tauerrval = 0.243651 #(m*e**lntau0)/(1e-5)**m, factor for statistical error propagation
if maxt == 1.8:
     lntau0 = 6.40
     mval=0.49
     tauerrval = 0.229462
if maxt == 2.6:
     lntau0 = 6.47
     tauerrval = 0.213949

tauar=[]
tauerrar=[]

cval150=0.00000038312
cval90=0.00000024007
cvalILC=0.00000036856 #doesn't matter because we are working in y, this just converts between recorded dT and y 
for i in range(len(dts18150)):
	print 'dt'
	print dts18150[i]
	yval = -cval150*dts18150[i] #y0 value from temperature decrement 
	print yval
	yvalerr = cval150*dts18150err[i] #y0 error from temperature decrement error
	print yvalerr
	#r500 = Da*t500*3.086*10**24
	#sarray = np.arange(-2.0,2.0,0.001)
	#thetaarray = np.arange(3e-6,maxt,1e-2)/60.0*deg2rad
	#pij = np.zeros([len(sarray),len(thetaarray)])
	#for i in range(len(sarray)):
	#    for j in range(len(thetaarray)):
	#        pij[i,j] = Pgnfw(sarray[i],thetaarray[j],r500,t500)
	#ps = np.zeros([len(thetaarray)])
	#for j in range(len(thetaarray)):
	#    ps[j] = np.trapz(pij[:,j],sarray)
	#ps = ps/max(ps)*yval
	#pint = np.trapz(ps*thetaarray,thetaarray)
	#yintd = 2.*pint/(maxt/60.*deg2rad)**2.0
	yintd=yval#/(np.pi*2.1**2.)
	#yintd=(yval*1.8**2.)/(2.1**2.)
	tau = e**(-lntau0+mval*log(yintd/10**(-5)))
	yerr = yvalerr*yintd/yval
	tauerr = tauerrval*yerr/(yintd**0.51)
	tauar.append(tau/10**(-4))
	tauerrar.append(tauerr/10**(-4))
	print 'y',yintd,'yerr',yerr,'tau',tau,'tauerr',tauerr

for i in range(len(tauar)):
	print tauar[i]
print 'errors:'
for i in range(len(tauar)):
	print tauerrar[i]
