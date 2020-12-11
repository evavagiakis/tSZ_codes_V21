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

submaps2=h5py.File('samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_cmb_map_v1.0.0_rc_joint.fits_submaps.h5','r')
submaps=h5py.File('samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_cmb_map_v1.0.0_rc_joint.fits_submaps.h5','r')
divmaps2=h5py.File('samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_cmb_map_v1.0.0_rc_joint.fits_divmaps.h5','r')
divmaps=h5py.File('samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_cmb_map_v1.0.0_rc_joint.fits_divmaps.h5','r')


print 'submap shapes:'
print np.shape(submaps2['submaps'])
print np.shape(submaps['submaps'])

# # BN and D56 regions: If incompsepregion==0, in neither region. If ==1, in D56 region. If ==2, in BN region. Taken based on average value of mask submap > 0.9
cpath='V20_DR15_Catalog_v2.csv'
catalog=pd.read_csv(cpath, comment = '#')
clum=np.array(catalog["lum"])
print len(clum)
incompsepregion=np.array(catalog["incompsepregion"])
S16ILC=np.array(catalog["S16ILC"])


deletear=[]
selectarbn=[]
selectard56=[]


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


counter=0
for l in [11.6e10,9.8e10,7.9e10,6.1e10,4.3e10]:
	counter=counter+1
	stacked_maps_tot=np.zeros(np.shape(np.array(submaps['submaps'][0,:,:],dtype='float')))
        stacked_maps_tot2=np.zeros(np.shape(np.array(submaps['submaps'][0,:,:],dtype='float')))
        divsarray=[]
       
        for i in range(len(clum)): #selecting lums of sources 
        	if (l<=clum[i] and S16ILC[i]==1 and incompsepregion[i]==1):# and (zd[i]<0.47354):# and divcut[i]==2 and pscut[i]==1):
                	divmapl=np.array(divmaps['divmaps'][i,:,:],dtype='float')
                        weight=divmapl[modrmap<r].mean()
                        stacked_maps_tot+=np.multiply(np.array(submaps['submaps'][i,:,:],dtype='float'),weight)   
                        divsarray.append(weight)

	for i in range(len(clum)): #selecting lums of sources 
                if (l<=clum[i] and S16ILC[i]==1 and incompsepregion[i]==2):# and (zd[i]<0.47354):# and divcut[i]==2 and pscut[i]==1):
		        divmapl=np.array(divmaps['divmaps'][i,:,:],dtype='float')
                        weight=divmapl[modrmap<r].mean()
                        stacked_maps_tot+=np.multiply(np.array(submaps2['submaps'][i,:,:],dtype='float'),weight)
                        divsarray.append(weight)

	print 'len divarray:'
	print len(divsarray)
#
#	stacked_maps_tot=stacked_maps_tot/np.sum(divsarray)
#
#        print stacked_maps_tot[modrmap<r].mean()-stacked_maps_tot[(modrmap>=r) & (modrmap<rout)].mean()
#
#        #normalize
#        stacked_maps_tot=np.subtract(stacked_maps_tot,stacked_maps_tot[modrmap<r].mean())
#
#        np.savetxt('Bin'+str(counter)+'_S16ILC_20201121.txt',stacked_maps_tot)
	

#counter=5
#for l in range(len(lupper)):
#        counter=counter+1
#
#        stacked_maps_tot=np.zeros(np.shape(np.array(submaps['submaps'][0,:,:],dtype='float')))
#        stacked_maps_tot2=np.zeros(np.shape(np.array(submaps['submaps'][0,:,:],dtype='float')))
#        divsarray=[]
#
#        for i in range(len(clum)): #selecting lums of sources 
#                if (llower[l+1]<=clum[i]<lupper[l] and S16ILC[i]==1 and incompsepregion[i]==1):# and (zd[i]<0.47354):# and divcut[i]==2 and pscut[i]==1):
#                        divmapl=np.array(divmaps['divmaps'][i,:,:],dtype='float')
#                        weight=divmapl[modrmap<r].mean()
#                        stacked_maps_tot+=np.multiply(np.array(submaps['submaps'][i,:,:],dtype='float'),weight)
#                        divsarray.append(weight)
#
#        for i in range(len(clum)): #selecting lums of sources 
#                if (llower[l+1]<=clum[i]<lupper[l] and S16ILC[i]==1 and incompsepregion[i]==2):# and (zd[i]<0.47354):# and divcut[i]==2 and pscut[i]==1):
#                        divmapl=np.array(divmaps['divmaps'][i,:,:],dtype='float')
#                        weight=divmapl[modrmap<r].mean()
#                        stacked_maps_tot+=np.multiply(np.array(submaps2['submaps'][i,:,:],dtype='float'),weight)
#                        divsarray.append(weight)
#
#        print 'len divarray:'
#        print len(divsarray)

        stacked_maps_tot=stacked_maps_tot/np.sum(divsarray)

        print stacked_maps_tot[modrmap<r].mean()-stacked_maps_tot[(modrmap>=r) & (modrmap<rout)].mean()

        #normalize
        #stacked_maps_tot=np.subtract(stacked_maps_tot,stacked_maps_tot[modrmap<r].mean())

        np.savetxt('Bin'+str(counter)+'_S16ILC_20201123.txt',stacked_maps_tot)


