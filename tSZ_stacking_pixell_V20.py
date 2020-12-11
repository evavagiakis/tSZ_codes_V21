import argparse
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
import numpy
import pylab
import pandas as pd
from pylab import *
from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel
from scipy import signal
from scipy.integrate import quad
from scipy.integrate import dblquad
from pixell import enmap
from pixell import reproject
import h5py

###################### adding arguments to run code ##############################
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--catalog", dest="catalog", help="cluster catalog to use", metavar="FILE")
parser.add_argument("-m", "--map", dest="maps", help="CMB maps to use", metavar="FILE", action = 'append')
args = parser.parse_args()
##################################################################################


###################### setting various parameters ################################ 
rad2deg = 57.2957795 # 1 radian in degrees
deg2rad = 0.01745 # deg to radian conversion factor
sigt = 0.00000000000000000000000066524587#in cm^2
melectron = 510.9989 #in keV/c^2
om = 0.315 #cosmological params to feed function later
ol = 0.685
h70 = 0.673
P0 = 8.403*h70**(3.0/2.0) #params from Hasselfield et al.
c500 = 1.177
gamma = 0.3081
alpha = 1.051
beta = 5.4905
clight = 299792 #speed of light in km/s
H0 = 67.3
t500 = 5.9/60.0*deg2rad #in radians
tmaxarcmin = 1.3 #radius of Nick's aperture in arcmin
thetamax = (tmaxarcmin/60.)*deg2rad #radius from Nick's aperture in radians 
##################################################################################


##################### defining functions #########################################
def deg2hr(degrees): #degrees to hours
    degrees_per_second = 15./360000.
    divisions = int(degrees/degrees_per_second) # seconds 
    hours = divisions / 360000.
    minutes = (divisions - 360000.*hours) / 6000.
    seconds = ((divisions - 360000.*hours) % 6000.)/100.
    return hours, minutes, seconds

def hr2deg(hours, minutes, seconds=0): #hours to degrees
    return (hours+minutes/60.0+seconds/3600.00)*15

def InvEz(z,H0,om,ol): #cosmology
    return 1.0/(np.sqrt(om*(1.0+z)**3.0 + ol)*H0)

def Pgnfw(s,t,r500,t500): #Pgnfw profile
    x = np.sqrt((s*3.086*10.**24.)**2.0+(r500*t/t500)**2.0)/r500
    ap = 0.22/(1.0+8.*x**3.0)
    Ez = np.sqrt(om*(1.0+z)**3.0 + ol)#*H0
    P500 = 0.00165*h70**2.*m**(2.0/3.0)*Ez**(8.0/3.0) #in keV/cm^3
    Px = P0*(c500*x)**(-gamma)*(1.0+(c500*x)**alpha)**((gamma-beta)/alpha)
    Pgnfw = P500*m**ap*Px
    return Pgnfw    

def coord2pix(ra,dec):
        pix = flipmap.skyToPix(ra,dec)
        pix=[round(pix[0])+0.5,round(pix[1])+0.5]
        return pix
###################################################################################

####################### Open the map you're using ##################################
path ='samba/actpol/map_coadd/20190813/out/act_planck_s08_s16_cmb_f150_night_map.fits'

#s18 new catalog10192020
#submapsLocation='samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f150_night_map.fits_submaps.h5'
#divmapsLocation='samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f150_night_map.fits_divmaps.h5'
submapsLocation='samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f090_night_map.fits_submaps.h5'
divmapsLocation='samba/V20_DR15Catalog_submaps/act_planck_s08_s18_cmb_f090_night_map.fits_divmaps.h5'

#submapsLocation='samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_cmb_map_v1.0.0_rc_joint.fits_submaps.h5'
#divmapsLocation='samba/V20_DR15Catalog_submaps/tilec_single_tile_boss_cmb_map_v1.0.0_rc_joint.fits_divmaps.h5'
#submapsLocation='samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_cmb_map_v1.0.0_rc_joint.fits_submaps.h5'
#divmapsLocation='samba/V20_DR15Catalog_submaps/tilec_single_tile_deep56_cmb_map_v1.0.0_rc_joint.fits_divmaps.h5'

psubmaps = h5py.File(submapsLocation, 'r')
pdivmaps = h5py.File(divmapsLocation,'r')
####################################################################################


##################### Load here your favorite catalog and mass/lum relations  ######
df=pd.read_csv(args.catalog, comment = '#')

rad=np.array(df["ra"])
decd=np.array(df["dec"])
lumd=np.array(df["lum"])
zd=np.array(df["z"])
#pscut=np.array(df["pscut"])#"PScut5arcmin"])
#pscut3=np.array(df["pscut3"])
#pscut2=np.array(df["PScut8arcmin"])
#divcut=np.array(df["divcut"])
#divcut90=np.array(df["divcut90"])
#galcut=np.array(df["galcut"])
stellarmasses=np.array(df["stellarmasses"])
halomasses=np.array(df["halomasses"])
#compsepmaps=np.array(df["compsepmaps"])
#coaddedmap=np.array(df["coaddedmap"])
S18coadd=np.array(df["S18coadd"])
incompsepregion=np.array(df["incompsepregion"])
S16ILC=np.array(df["S16ILC"])
use_PS90_and_S18coadd=np.array(df["use_PS90_and_S18coadd"])
#useforcompsepanalysis=np.array(df["useforcompsepanalysis"])

####################################################################################

################### Select sources by richness or luminosity #######################
lumselect, raselect, decselect, zselect,foreg,inds = [], [], [], [], [],[]

for i in range(len(lumd)): #selecting lums of sources
	if (4.3e10<=lumd[i] and use_PS90_and_S18coadd[i]==1):# and (zd[i]<0.47354):# and divcut[i]==2 and pscut[i]==1):
		lumselect.append(lumd[i])
        	raselect.append(rad[i])
        	decselect.append(decd[i])
        	zselect.append(zd[i])
		inds.append(i)

dec_rad,ra_rad=np.deg2rad(np.array((decselect,raselect)))
#print 'how many sources in d56 area:'
#print len([x for x in useforcompsepanalysis if x==1])
#print 'how many sources in bn area:'
#print len([x for x in useforcompsepanalysis if x==2])
#
#np.savetxt('useforcompsepanalysisflags.txt',useforcompsepanalysis)
#stop

#pixels=enmap.sky2pix(enmapmap.shape,enmapmap.wcs,np.vstack((decselect,raselect))*np.pi/180.,safe=True)
#print 'the length for reg'
#print len(decselect)
#rapi=pixels[1,:]
#decpi=pixels[0,:]
#coords=np.array([rapi,decpi])
#coords=coords.T
#np.savetxt('CompSepBNarea20200528.txt',coords,delimiter=',',newline='\n')
#stop


print 'Nsource cleaned and cropped:',len(zselect)
print '<z> cleaned and cropped:', np.mean(zselect)
print '<l>',np.mean(lumselect)
print min(lumselect)
print max(lumselect)
##################################################################################

########### Select radii and f (mult factor for outer radius) for AP filtering or pix selection ########
fmax = sqrt(2)
radmax = 2.1 #arcmin 
radbeam= 0.0 #arcmin

radius = float(radmax)
radiusbeam = float(radbeam)
mult = float(fmax)
radius_out = radius*mult

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

	submapl=np.array(psubmaps['submaps'][inds[l],:,:],dtype='float')
	divmapl=np.array(pdivmaps['divmaps'][inds[l],:,:],dtype='float')

	disks.append(submapl[(modrmap>=rbeam) & (modrmap<r)].mean())
        rings.append(submapl[(modrmap>=r) & (modrmap<rout)].mean())
        disk_stds.append(submapl[modrmap<r].std())
        ring_stds.append(submapl[(modrmap>=r) & (modrmap<rout)].std())
	divs.append(divmapl[modrmap<r].mean())

print 'Length of array'
print len(disks)

np.savetxt('S18f090_2p1arcmin_V20cat_additionalf090PScutPato_20201209.txt',np.array([raselect,decselect,lumselect,zselect,disks,disk_stds,rings,ring_stds,divs]))

