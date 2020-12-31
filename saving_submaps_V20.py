import argparse
import numpy as np
import pandas as pd
import h5py

###################### Adding arguments to run code ##############################
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--catalog", dest="catalog", help="cluster catalog to use", metavar="FILE")
parser.add_argument("-a","--analysis",dest="analysis",help="what analysis to run",metavar="FILE",action='append')
args = parser.parse_args()
###################################################################################

###################### Setting various parameters ################################ 
deg2rad = 0.01745 # deg to radian conversion factor
##################################################################################


######################### Submap locations #######################################
if (args.analysis[0] == 'DR5f150'):
	submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f150_night_map.fits_submaps.h5'
	divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f150_night_map.fits_divmaps.h5'
if (args.analysis[0] == 'DR5f090'):
	submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f090_night_map.fits_submaps.h5'
	divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f090_night_map.fits_divmaps.h5'
if (args.analysis[0] == 'ILC'):
	submapsLocation2='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_cmb_map_v1.2.0_joint.fits_submaps.h5'
	divmapsLocation2='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_cmb_map_v1.2.0_joint.fits_divmaps.h5'
	submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_cmb_map_v1.2.0_joint.fits_submaps.h5'
	divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_cmb_map_v1.2.0_joint.fits_divmaps.h5'

psubmaps = h5py.File(submapsLocation, 'r')
pdivmaps = h5py.File(divmapsLocation,'r')

if (args.analysis[0] == 'ILC'):
	psubmaps2 = h5py.File(submapsLocation2, 'r')
	pdivmaps2 = h5py.File(divmapsLocation2,'r')
####################################################################################


##################### Load your favorite catalog ####################################
df=pd.read_csv(args.catalog, comment = '#')

rad=np.array(df["ra"])
decd=np.array(df["dec"])
lumd=np.array(df["lum"])
zd=np.array(df["z"])
S18coadd=np.array(df["S18coadd"])
incompsepregion=np.array(df["incompsepregion"])
S16ILC=np.array(df["S16ILC"])
####################################################################################
fmax = np.sqrt(2)
radmax = 2.1 #arcmin 
radius = float(radmax)
mult = float(fmax)
radius_out = radius*mult
r=np.deg2rad(radius/60.)
rout=np.deg2rad(radius_out/60.)
modrmap=psubmaps['modr']


larray=[11.6e10,9.8e10,7.9e10,6.1e10,4.3e10]

counter=0

print 'hey'

if (args.analysis[0] == 'DR5f150' or args.analysis[0] == 'DR5f090'):
	for l in range(11):
		counter=counter+1
		stacked_maps_tot=np.zeros(np.shape(np.array(psubmaps['submaps'][0,:,:],dtype='float')))
		divsarray=[]
	
		if l<5:
			for i in range(len(lumd)):
				if (larray[l]<=lumd[i] and S18coadd[i]==1):
					divmapl=np.array(pdivmaps['divmaps'][i,:,:],dtype='float')
	                        	weight=divmapl[modrmap<r].mean()
	                        	stacked_maps_tot+=np.multiply(np.array(psubmaps['submaps'][i,:,:],dtype='float'),weight)
	                        	divsarray.append(weight)
	
		if 5<=l:
			for i in range(len(lumd)):
				if (larray[l+1-5]<=lumd[i]<larray[l-5] and S18coadd[i]==1):
					divmapl=np.array(pdivmaps['divmaps'][i,:,:],dtype='float')
	                                weight=divmapl[modrmap<r].mean()
	                                stacked_maps_tot+=np.multiply(np.array(psubmaps['submaps'][i,:,:],dtype='float'),weight)
	                                divsarray.append(weight)	
		stacked_maps_tot=stacked_maps_tot/np.sum(divsarray)
		
		print 'Nsources:',len(divsarray)
		print stacked_maps_tot[modrmap<r].mean()-stacked_maps_tot[(modrmap>=r) & (modrmap<rout)].mean() #Check value if you want
	
		np.savetxt('Bin'+str(counter)+'_DR5f090_20201228.txt',stacked_maps_tot)
	
if (args.analysis[0] == ('ILC')):
        for l in range(11):
                counter=counter+1
                stacked_maps_tot=np.zeros(np.shape(np.array(psubmaps['submaps'][0,:,:],dtype='float')))
                stacked_maps_tot2=np.zeros(np.shape(np.array(psubmaps2['submaps'][0,:,:],dtype='float')))
		divsarray=[]

                if l<5:
                        for i in range(len(lumd)):
                                if (larray[l]<=lumd[i] and S16ILC[i]==1 and incompsepregion[i]==1):
                                        divmapl=np.array(pdivmaps['divmaps'][i,:,:],dtype='float')
                                        weight=divmapl[modrmap<r].mean()
                                        stacked_maps_tot+=np.multiply(np.array(psubmaps['submaps'][i,:,:],dtype='float'),weight)
                                        divsarray.append(weight)
				if (larray[l]<=lumd[i] and S16ILC[i]==1 and incompsepregion[i]==2):
					divmapl=np.array(pdivmaps2['divmaps'][i,:,:],dtype='float')
                                        weight=divmapl[modrmap<r].mean()
                                        stacked_maps_tot+=np.multiply(np.array(psubmaps2['submaps'][i,:,:],dtype='float'),weight)
                                        divsarray.append(weight)

                if 5<=l:
                        for i in range(len(lumd)):
                                if (larray[l+1]<=lumd[i]<larray[l] and S16ILC[i]==1 and incompsepregion[i]==1):
                                        divmapl=np.array(pdivmaps['divmaps'][i,:,:],dtype='float')
                                        weight=divmapl[modrmap<r].mean()
                                        stacked_maps_tot+=np.multiply(np.array(psubmaps['submaps'][i,:,:],dtype='float'),weight)
                                        divsarray.append(weight)
				if (larray[l+1]<=lumd[i]<larray[l] and S16ILC[i]==1 and incompsepregion[i]==2):                
					divmapl=np.array(pdivmaps2['divmaps'][i,:,:],dtype='float')
                                        weight=divmapl[modrmap<r].mean()
                                        stacked_maps_tot+=np.multiply(np.array(psubmaps2['submaps'][i,:,:],dtype='float'),weight)
                                        divsarray.append(weight)

		stacked_maps_tot=stacked_maps_tot/np.sum(divsarray)
		
		print 'Nsources:',len(divsarray)
                print stacked_maps_tot[modrmap<r].mean()-stacked_maps_tot[(modrmap>=r) & (modrmap<rout)].mean() #Check value if you want

                np.savetxt('Bin'+str(counter)+'_DR4ILC_20201123.txt',stacked_maps_tot)

	

