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
if (args.analysis[0] == 'ILCBN'):
	submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_comptony_map_v1.2.0_joint.fits_submaps.h5'
	divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_comptony_map_v1.2.0_joint.fits_divmaps.h5'
if (args.analysis[0] == 'ILCD56'):
	submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_comptony_map_v1.2.0_joint.fits_submaps.h5'
	divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_comptony_map_v1.2.0_joint.fits_divmaps.h5'



psubmaps = h5py.File(submapsLocation, 'r')
pdivmaps = h5py.File(divmapsLocation,'r')
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

############################### Select sources #####################################
lumselect, raselect, decselect, zselect, inds = [], [], [], [],[]


if (args.analysis[0] == 'DR5f150' or args.analysis[0] == 'DR5f090'):
	for i in range(len(lumd)): #selecting sources
	        if (4.3e10<=lumd[i] and S18coadd[i]==1):
	                inds.append(i)

if (args.analysis[0] == 'ILCD56'):
	for i in range(len(lumd)): #selecting sources
	        if (4.3e10<=lumd[i] and S16ILC[i]==1 and incompsepregion[i]==1):
	                inds.append(i)

if (args.analysis[0] == 'ILCBN'):
	for i in range(len(lumd)): #selecting sources
		if (4.3e10<=lumd[i] and S16ILC[i]==1 and incompsepregion[i]==2):
			inds.append(i)

lumselect=[lumd[i] for i in inds]
raselect=[rad[i] for i in inds]
decselect=[decd[i] for i in inds]
zselect=[zd[i] for i in inds]

print 'Number of sources selected:',len(zselect)
print '<z>:', np.mean(zselect)
print '<l>:',np.mean(lumselect)
print min(lumselect)
print max(lumselect)
##################################################################################

########### Select radii for AP filtering ########################################
fmax = np.sqrt(2) #factor which defines outer AP radius
radmax = 2.1 #arcmin 
radbeam= 0.0#1.05 #arcmin, radius of any beam you'll use for core excising

radius = float(radmax)
radiusbeam = float(radbeam)
mult = float(fmax)
radius_out = radius*mult #defining AP radii 
########################################################################################################

### Getting data from maps and calculating temperatures  ################################################ 
disks, rings, disk_stds, ring_stds, divs = [], [], [], [], []
r=np.deg2rad(radius/60.)
rbeam=np.deg2rad(radiusbeam/60.)
rout=np.deg2rad(radius_out/60.)
modrmap=psubmaps['modr']

print 'Creating submaps...'
for l in range(len(raselect)):
	if (l == int(.1*len(raselect))):
	    print '10% complete'
	if (l == int(.2*len(raselect))):
	    print '20% complete'
	if (l == int(.3*len(raselect))):
	    print '30% complete'
	if (l == int(.4*len(raselect))):
	    print '40% complete'
	if (l == int(.5*len(raselect))):
	    print '50% complete'
	if (l == int(.6*len(raselect))):
	    print '60% complete'
	if (l == int(.7*len(raselect))):
	    print '70% complete'
	if (l == int(.8*len(raselect))):
	    print '80% complete'
	if (l == int(.9*len(raselect))):
	    print '90% complete'
	if (l == int(1.*len(raselect))):
	    print 'done'

	submapl=np.array(psubmaps['submaps'][inds[l],:,:],dtype='float') #make arrays from submaps and dive maps
	divmapl=np.array(pdivmaps['divmaps'][inds[l],:,:],dtype='float')

	disks.append(submapl[(modrmap>=rbeam) & (modrmap<r)].mean()) #record averages within disks and rings
        rings.append(submapl[(modrmap>=r) & (modrmap<rout)].mean())
        disk_stds.append(submapl[(modrmap>=rbeam) & (modrmap<r)].std()) #record standard deviation within disk and ring
        ring_stds.append(submapl[(modrmap>=r) & (modrmap<rout)].std())
	divs.append(divmapl[modrmap<r].mean()) #record average noise value in disk

print 'Length of array'
print len(disks)

np.savetxt('ILCD56_20201223.txt',np.array([raselect,decselect,lumselect,zselect,disks,disk_stds,rings,ring_stds,divs]))

