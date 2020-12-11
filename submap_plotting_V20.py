import h5py
import numpy as np
import math
import matplotlib.pyplot as plt
from enlib import enmap
import pandas as pd

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

cpath='V20_DR15_Catalog_v2.csv'
catalog=pd.read_csv(cpath, comment = '#')
clum=np.array(catalog["lum"])
S18coadd=np.array(catalog["S18coadd"])

submapsLocation='samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f150_night_map.fits_submaps.h5'
divmapsLocation='samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f150_night_map.fits_divmaps.h5'

submaps = h5py.File(submapsLocation, 'r')
divmaps = h5py.File(divmapsLocation,'r')

fmax = np.sqrt(2)
radmax = 2.1 #arcmin 
radius = float(radmax)
mult = float(fmax)
radius_out = radius*mult
r=np.deg2rad(radius/60.)
rout=np.deg2rad(radius_out/60.)
modrmap=submaps['modr']


lupper=[11.6e10,9.8e10,7.9e10,6.1e10,4.3e10]
llower=[11.6e10,9.8e10,7.9e10,6.1e10,4.3e10]


#counter=0
#for l in [11.6e10,9.8e10,7.9e10,6.1e10,4.3e10]:
#	counter=counter+1
#	
#	stacked_maps_tot=np.zeros(np.shape(np.array(submaps['submaps'][0,:,:],dtype='float')))
#	divsarray=[]
#	
#	for i in range(len(clum)): #selecting lums of sources 
#		if (l<=clum[i] and S18coadd[i]==1):# and (zd[i]<0.47354):# and divcut[i]==2 and pscut[i]==1):
#                	#lumselect.append(lumd[i])
#                	#raselect.append(rad[i])
#                	#decselect.append(decd[i])
#                	#zselect.append(zd[i])
#                	#inds.append(i)                     
#			divmapl=np.array(divmaps['divmaps'][i,:,:],dtype='float')
#			weight=divmapl[modrmap<r].mean()
#			stacked_maps_tot+=np.multiply(np.array(submaps['submaps'][i,:,:],dtype='float'),weight)   
#			divsarray.append(weight)
## 
counter=5
for l in range(len(lupper)):
        counter=counter+1

        stacked_maps_tot=np.zeros(np.shape(np.array(submaps['submaps'][0,:,:],dtype='float')))
        divsarray=[]

        for i in range(len(clum)): #selecting lums of sources 
                if (llower[l+1]<=clum[i]<lupper[l] and S18coadd[i]==1):# and (zd[i]<0.47354):# and divcut[i]==2 and pscut[i]==1):
                        #lumselect.append(lumd[i])
                        #raselect.append(rad[i])
                        #decselect.append(decd[i])
                        #zselect.append(zd[i])
                        #inds.append(i)                     
                        divmapl=np.array(divmaps['divmaps'][i,:,:],dtype='float')
                        weight=divmapl[modrmap<r].mean()
                        stacked_maps_tot+=np.multiply(np.array(submaps['submaps'][i,:,:],dtype='float'),weight)
                        divsarray.append(weight)	

#############################





#
#
	stacked_maps_tot=stacked_maps_tot/np.sum(divsarray)
#	
	print stacked_maps_tot[modrmap<r].mean()-stacked_maps_tot[(modrmap>=r) & (modrmap<rout)].mean()
#	
#	#normalize
#	#stacked_maps_tot=np.subtract(stacked_maps_tot,stacked_maps_tot[modrmap<r].mean())
#
	np.savetxt('Bin'+str(counter)+'_S18_150GHz_20201123.txt',stacked_maps_tot)

	#plt.imshow(stacked_maps_tot)
	#cbar=plt.colorbar()
	#cbar.set_label('[uK]',fontsize=16)
	#cbar.ax.tick_params(labelsize=16)
	#circle = np.rad2deg(submaps['modr'])*60.
	#plt.contour(circle, [2.1,2.9698484802])
	#plt.ylabel('arcmin',fontsize=16)
	#plt.xlabel('arcmin',fontsize=16)
	#plt.title('16.3<L',fontsize=16)
	#plt.xticks(fontsize=16)
	#plt.yticks(fontsize=16)
	#plt.xticks(np.arange(0, 200, step=25),np.arange(0, 20.0, step=2.5))
	#plt.yticks(np.arange(0, 200, step=25),np.arange(0, 20.0, step=2.5))
	#plt.gca().invert_yaxis()
	##plt.savefig('Bin'+str(counter)+'_20200504')
	#plt.show()
	










