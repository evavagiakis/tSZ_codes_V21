import argparse
import numpy as np
import pandas as pd
import h5py

###################### Adding arguments to run code ##############################
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--catalog", dest="catalog", help="cluster catalog to use", metavar="FILE")
args = parser.parse_args()
###################################################################################

###################### Setting various parameters ################################ 
deg2rad = 0.01745 # deg to radian conversion factor
##################################################################################


######################### Submap locations #######################################
#submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f150_night_map.fits_submaps.h5'
#divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f150_night_map.fits_divmaps.h5'
#submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f090_night_map.fits_submaps.h5'
#divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f090_night_map.fits_divmaps.h5'

#submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_cmb_map_v1.2.0_joint.fits_submaps.h5'
#divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_cmb_map_v1.2.0_joint.fits_divmaps.h5'
submapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_cmb_map_v1.2.0_joint.fits_submaps.h5'
divmapsLocation='/home/ev66/samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_cmb_map_v1.2.0_joint.fits_divmaps.h5'

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

################### Select sources by richness or luminosity #######################
lumselect, raselect, decselect, zselect,foreg,inds = [], [], [], [], [],[]

for i in range(len(lumd)): #selecting lums of sources
	if (4.3e10<=lumd[i] and S16ILC[i]==1 and incompsepregion[i]==1):
		lumselect.append(lumd[i])
        	raselect.append(rad[i])
        	decselect.append(decd[i])
        	zselect.append(zd[i])
		inds.append(i)

#dec_rad,ra_rad=np.deg2rad(np.array((decselect,raselect)))

print 'Nsource cleaned and cropped:',len(zselect)
print '<z> cleaned and cropped:', np.mean(zselect)
print '<l>',np.mean(lumselect)
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


####################Printing our clean patch bin info, making submaps ###########################################
print 'Clean patch lums', 'N',len(raselect),'<L>',np.mean(lumselect),'<z>',np.mean(zselect)

### Getting data from maps and calculating temperatures  ################################################ 
disks=[]
rings=[]
disk_stds=[]
ring_stds=[]
divs=[]
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
        disk_stds.append(submapl[modrmap<r].std()) #record standard deviation within disk and ring
        ring_stds.append(submapl[(modrmap>=r) & (modrmap<rout)].std())
	divs.append(divmapl[modrmap<r].mean()) #record average noise value in disk

print 'Length of array'
print len(disks)

np.savetxt('DR4ILC_tilec_single_tile_deep56_cmb_map_v1.2.0_joint_AP_v3cat_20201223.txt',np.array([raselect,decselect,lumselect,zselect,disks,disk_stds,rings,ring_stds,divs]))

