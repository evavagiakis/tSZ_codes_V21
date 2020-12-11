import numpy as np
import math
import matplotlib.pyplot as plt
from enlib import enmap


def get_group(i,Ngroups,population):
	lenGroup=len(population)/Ngroups
	sel=np.ones(len(population),dtype=bool)
	#sel[i*int(lenGroup):(i+1)*int(lenGroup)]=False
	if i==0:
		sel[0]=False
	else:
		sel[i*int(lenGroup):(i+1)*int(lenGroup)]=False
	return sel

def estimatorFunction(dt,divsmap,sel):
	return np.sum(np.multiply(dt[sel],divsmap[sel]))/np.sum(divsmap[sel])

data=np.loadtxt('S18f090_2p1arcmin_V20cat_additionalf090PScutPato_20201209.txt')

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


for m in [1]:#[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]:#,1.0]:
	sn=[]
	dts=[]
	errs=[]
	intsoverbins=[]
	for i in [0,1,2,3,4,10,11,12,13]:
		if i <5:
			ras=[ra[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]# and divs[x]>=(1./(45.**2.))]# and 148.<ra[x]<244. and -4.<dec[x]<21.]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]  
			decs=[dec[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]# and divs[x]>=(1./(45.**2.))]# and 148.<ra[x]<244. and -4.<dec[x]<21.]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
			zs=[z[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]# and divs[x]>=(1./(45.**2.))]# and 148.<ra[x]<244. and -4.<dec[x]<21.]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
			lums=[lum[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]# and divs[x]>=(1./(45.**2.))]# and 148.<ra[x]<244. and -4.<dec[x]<21.]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
			dt=np.array([signals[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]])# and divs[x]>=(1./(45.**2.))])# and 148.<ra[x]<244. and -4.<dec[x]<21.])# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]])
			divsbin=[divs[x] for x in range(len(lum)) if lum[x]>=lumcuts2[i]]# and divs[x]>=(1./(45.**2.))]# and 148.<ra[x]<244. and -4.<dec[x]<21.]
			#pixarray = enmap.sky2pix(enmapmap.shape,enmapmap.wcs,np.vstack((decs,ras))*np.pi/180.,safe=True)
			#raselectpixels=pixarray[1,:]
			#decselectpixels=pixarray[0,:]
			#coords=np.array([raselectpixels,decselectpixels])
			#coords=coords.T
			#np.savetxt('coords_bin7_m='+str(m)+'.txt',coords,delimiter=',',newline='\n')

			#ints=[fgs[x] for x in range(len(lum)) if lum[x]>=lumcuts[i] and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
#			plt.hist(zs,30,alpha=0.5,label='m='+str(m))
		if i==10:
			ras=[ra[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
	                decs=[dec[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]# and divs[x]>=(1./(45.**2.))]#  and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
	                zs=[z[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
	                lums=[lum[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
	                dt=np.array([signals[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]])# and divs[x]>=(1./(45.**2.))])# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]])
			divsbin=[divs[x] for x in range(len(lum)) if lumcuts[1]<=lum[x]<lumcuts[0]]# and divs[x]>=(1./(45.**2.))]
			#ints=[fgs[x] for x in range(len(lum)) if lum[3]>=lumcuts[1] and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
#			plt.hist(zs,30,alpha=0.5,label='m='+str(m))
		if i==11:
			ras=[ra[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
	                decs=[dec[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]# and divs[x]>=(1./(45.**2.))]#  and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
	                zs=[z[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
	                lums=[lum[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
	                dt=np.array([signals[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]])#) and divs[x]>=(1./(45.**2.))])
			divsbin=[divs[x] for x in range(len(lum)) if lumcuts[2]<=lum[x]<lumcuts[1]]# and divs[x]>=(1./(45.**2.))]
			# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]])
			#ints=[fgs[x] for x in range(len(lum)) if lum[7]>=lumcuts[3] and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
		if i==12:
                        ras=[ra[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
                        decs=[dec[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]# and divs[x]>=(1./(45.**2.))]#  and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
                        zs=[z[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
                        lums=[lum[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
                        dt=np.array([signals[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]])# and divs[x]>=(1./(45.**2.))])# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]])
			divsbin=[divs[x] for x in range(len(lum)) if lumcuts[3]<=lum[x]<lumcuts[2]]# and divs[x]>=(1./(45.**2.))]
		if i==13:
                        ras=[ra[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
                        decs=[dec[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]# and divs[x]>=(1./(45.**2.))]#  and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
                        zs=[z[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
                        lums=[lum[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]# and divs[x]>=(1./(45.**2.))]# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]]
                        dt=np.array([signals[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]])# and divs[x]>=(1./(45.**2.))])# and fgs[x]<fgssorted[int(len(fgssorted)*m)-1]])
			divsbin=[divs[x] for x in range(len(lum)) if lumcuts[4]<=lum[x]<lumcuts[3]]# and divs[x]>=(1./(45.**2.))]
	#intsoverbins.append(np.mean(ints))
	#plt.plot(intsoverbins,'--o',label='Bin'+str(i))	
#			plt.hist(zs,30,alpha=0.5,label='m='+str(m))
		#errso=[]	
		#for u in [2,5,10,25,50,75,100,150,250,350,500,650,800,950,1100,1250,1400,1550,1700,1850,2000]:
		print len(dt)
		print len(divsbin)
		divsbin=np.array(divsbin)
		avgs=[estimatorFunction(dt,divsbin,get_group(j,500,dt)) for j in range(len(ras))]
		err=np.sqrt(((len(avgs)-1.)/len(avgs))*np.sum((avgs-(np.sum(np.multiply(dt,divsbin))/np.sum(divsbin)))**2.0))
		#	errso.append(err)
	
		dtl=np.sum(np.multiply(dt,divsbin))/np.sum(divsbin)

		print 'BIN ', i
		print '-----------------------------'
	        print 'N=',len(zs)
	        print '<z>=',np.mean(zs)
		print 'dt',dtl #np.sum(np.multiply(dt,divsbin))/np.sum(divsbin)
		print 'jk err',err
		#print 'div',np.mean(divsbin)
#		print 'S/N',np.mean(dt)/err
#		sn.append(-np.mean(dt)/err)
		dts.append(dtl)
		errs.append(err)

