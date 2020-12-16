import numpy as np
import math
import matplotlib.pyplot as plt
#from enlib import enmap


def get_group(i,Ngroups,population):
	lenGroup=len(population)/Ngroups
	sel=np.ones(len(population),dtype=bool)
	if i==0:
		sel[0]=False
	else:
		sel[i*int(lenGroup):(i+1)*int(lenGroup)]=False
	return sel

def estimatorFunction(dt,divsmap,sel):
	return np.sum(np.multiply(dt[sel],divsmap[sel]))/np.sum(divsmap[sel])

data=np.loadtxt('S18f150_donut_duplicationcheck.txt')

ra=data[0,:]
dec=data[1,:]
lum=data[2,:]
z=data[3,:]
disk_mean=data[4,:]
disk_std=data[5,:]
annulus_mean=data[6,:]
annulus_std=data[7,:]
divs=data[8,:]

signals=np.subtract(disk_mean,annulus_mean)

lumcuts=[11.6e10,9.8e10,7.9e10,6.1e10,4.3e10]
lumcuts2=[11.6e10,9.8e10,7.9e10,6.1e10,4.3e10]


for m in [1]:
	sn=[]
	dts=[]
	errs=[]
	intsoverbins=[]
	for i in [0,1,2]:#,3,4,10,11,12,13]:
		if i <5:
			ras=[ra[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]
			decs=[dec[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]
			zs=[z[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]
			lums=[lum[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]
			dt=np.array([signals[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]])
			divsbin=[divs[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]
			
		if i==10:
			ras=[ra[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]] 
	                decs=[dec[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]
	                zs=[z[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]
	                lums=[lum[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]
	                dt=np.array([signals[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]])
			divsbin=[divs[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]
			
		if i==11:
			ras=[ra[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]
	                decs=[dec[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]
	                zs=[z[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]
	                lums=[lum[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]
	                dt=np.array([signals[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]])
			divsbin=[divs[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]
		if i==12:
                        ras=[ra[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]
                        decs=[dec[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]
                        zs=[z[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]
                        lums=[lum[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]
                        dt=np.array([signals[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]])
			divsbin=[divs[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]
		if i==13:
                        ras=[ra[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]] 
                        decs=[dec[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]
                        zs=[z[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]
                        lums=[lum[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]
                        dt=np.array([signals[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]])
			divsbin=[divs[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]

		print len(dt)
		print len(divsbin)
		divsbin=np.array(divsbin)
		avgs=[estimatorFunction(dt,divsbin,get_group(j,500,dt)) for j in range(len(ras))]
		err=np.sqrt(((len(avgs)-1.)/len(avgs))*np.sum((avgs-(np.sum(np.multiply(dt,divsbin))/np.sum(divsbin)))**2.0))
		dtl=np.sum(np.multiply(dt,divsbin))/np.sum(divsbin)

		print 'BIN ', i
		print '-----------------------------'
	        print 'N=',len(zs)
	        print '<z>=',np.mean(zs)
		print 'dt',dtl 
		print 'jk err',err
		dts.append(dtl)
		errs.append(err)

np.savetxt('testoutput.txt',np.array([dts,errs]))

