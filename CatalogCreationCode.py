import argparse
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
import numpy
import pandas as pd
import pylab
from pylab import *
from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel
from scipy import signal
from scipy.integrate import quad
from scipy.integrate import dblquad
from scipy import interpolate
#from flipper import liteMap
from enlib import enmap
from enlib import coordinates
import enlib
from orphics import io
#from orphics import catalogs 
from orphics import maps
import healpy as hp
from astropy_healpix import HEALPix
from astropy import units as u
from astropy.coordinates import SkyCoord

########################################## Mass Relationships ######################################################
#Variables from Kravtsov Appendix A Table 3 [arXiv 1401.7329]
log10M1 = 11.39
log10e = -1.685
alpham = -1.740   # minus sign missing in Table 3
deltam = 4.335
gammam = 0.531

def f(x):
      """Fitting function, Kravtsov+14, eq A4
      """
      result = np.log10(1.+np.exp(x))**gammam
      result *= deltam
      result /= 1. + np.exp(10.**(-x))
      result += -np.log10(10.**(alpham*x) + 1.)
      return result

def fmStar(mVir):
      """Computes stellar mass [M_sun]
      from halo mass Mvir [Msun].
      """
      result = log10e + log10M1
      result += f(np.log10(mVir) - log10M1)
      result -= f(0.)
      result = 10.**result
      return result

#To get function to interpolate 
mVir_sample = np.logspace(np.log10(1.e5), np.log10(1.e20), 501, 10.)
mStar_sample = np.array(map(fmStar, mVir_sample))  # [M_sun]

#Interpolating
fmStarTomVir = interpolate.interp1d(mStar_sample, mVir_sample, kind='cubic', bounds_error=True, fill_value=0.)
fmVirTomStar = interpolate.interp1d(mVir_sample, mStar_sample, kind='cubic', bounds_error=True, fill_value=0.)
###################################################################################################################

###################################### Planck Galactic Plane Masks ################################################
#nside=2048
#hel=HEALPix(nside=nside,order='ring') #change for appropriate order from map header
#
#hmap_40=hp.read_map('product-action?MAP.MAP_ID=COM_CompMap_Compton-SZMap-masks_2048_R2.01.fits',field=0)#HFI_Mask_GalPlane-apo0_2048_R2.00.fits',field=8) #load healpix map
#hmap_50=hp.read_map('product-action?MAP.MAP_ID=COM_CompMap_Compton-SZMap-masks_2048_R2.01.fits',field=1)
###################################################################################################################

#Load catalogs and maps to work with 
original=pd.read_csv('samba/DR15_actplanck_catalog_wbestObjID_PetrANDcModel_20200902_EMV_evavagiakis.csv',comment='#') #DR15 catalog downloaded from casjobs
kcorrected=np.loadtxt('DR15_actplanck_catalog_wbestObjID_PetrANDcModel_20200902_EMV_evavagiakis_Kcorrected_PAG20200913.csv',skiprows=1) #kcorrected catalog produced by Pato
pathdiv='samba/actpol/map_coadd/20190813/out/act_planck_s08_s16_cmb_f150_night_div.fits' #Season 16 div map for div cuts
#pathdiv90='cosmo/kSZ/actpol/map_coadd/20190813/out/act_planck_s08_s16_cmb_f090_night_div.fits' #Season 16 90 GHz div map for div cuts
pathdivs18='samba/actpol/map_coadd/20200228/release/act_planck_s08_s18_cmb_f150_night_ivar.fits' #Season 18 div map for div cuts 
pathPS15mJy='samba/actpol/map_coadd/masks_20190826/mask_s13s16_0.015mJy_10.0arcmin_20190826_ext.fits' #PS mask 15 mJy 10 arcmin
pathPS100mJy='samba/actpol/map_coadd/masks_20190826/mask_s13s16_0.1mJy_5.0arcmin_20190826_ext.fits' #PS mask 100 mJy 5 arcmin
pathbnmask = 'samba/actpol/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/tilec_mask.fits'
pathd56mask = 'samba/actpol/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/tilec_mask.fits'
PSadditional = pd.read_csv('samba/PAG20200416_DR15_point_source_removal_PAG09232020.csv')
V20old= pd.read_csv('V20_DR15_Catalog.csv',comment='#')

#Define arrays for catalog 
ra=np.array(V20old["ra"])
dec=np.array(V20old["dec"])
z=np.array(V20old["z"])
bestObjID=np.array(V20old["bestObjID"])
cModelMag_u=np.array(V20old["cModelMag_u"])
cModelMag_g=np.array(V20old["cModelMag_g"])
cModelMag_r=np.array(V20old["cModelMag_r"])
cModelMag_i=np.array(V20old["cModelMag_i"])
cModelMag_z=np.array(V20old["cModelMag_z"])
cModelMagErr_u=np.array(V20old["cModelMagErr_u"])
cModelMagErr_g=np.array(V20old["cModelMagErr_g"])
cModelMagErr_r=np.array(V20old["cModelMagErr_r"])
cModelMagErr_i=np.array(V20old["cModelMagErr_i"])
cModelMagErr_z=np.array(V20old["cModelMagErr_z"])
petroMag_u=np.array(V20old["petroMag_u"])
petroMag_g=np.array(V20old["petroMag_g"])
petroMag_r=np.array(V20old["petroMag_r"])
petroMag_i=np.array(V20old["petroMag_i"])
petroMag_z=np.array(V20old["petroMag_z"])
petroMagErr_u=np.array(V20old["petroMagErr_u"])
petroMagErr_g=np.array(V20old["petroMagErr_g"])
petroMagErr_r=np.array(V20old["petroMagErr_r"])
petroMagErr_i=np.array(V20old["petroMagErr_i"])
petroMagErr_z=np.array(V20old["petroMagErr_z"])
extinction_u=np.array(V20old["extinction_u"])
extinction_g=np.array(V20old["extinction_g"])
extinction_r=np.array(V20old["extinction_r"])
extinction_i=np.array(V20old["extinction_i"])
extinction_z=np.array(V20old["extinction_z"])

use_PS=np.array(V20old["use_PS"])#(PSadditional["use_PS"])

lum=np.array(V20old["lum"])
#lum=kcorrected[:,3]

#stellarmasses=np.array(V20old["stellarmasses"])
#halomasses=np.array(V20old["halomasses"])
PS15mJy_cut=np.array(V20old["PS15mJy_cut"])
PS100mJy_cut=np.array(V20old["PS100mJy_cut"])
divcut=np.array(V20old["divcut"])
divcuts18=np.array(V20old["divcuts18"])
galcut=np.array(V20old["galcut"])
incompsepregion=np.array(V20old["incompsepregion"])
S18coadd=np.array(V20old["S18coadd"])
S16ILC=np.array(V20old["S16ILC"])

########### Estimate masses
stellarmasses=[3.0*lum[x] for x in range(len(lum))]
halomasses=fmStarTomVir(stellarmasses)
#print halomasses
#stop

#############Perform gal cuts
#nside=2048
#hel=HEALPix(nside=nside,order='ring') #change for appropriate order from map header
#
#hmap_40=hp.read_map('product-action?MAP.MAP_ID=COM_CompMap_Compton-SZMap-masks_2048_R2.01.fits',field=0)#HFI_Mask_GalPlane-apo0_2048_R2.00.fits',field=8) #load healpix map
#hmap_50=hp.read_map('product-action?MAP.MAP_ID=COM_CompMap_Compton-SZMap-masks_2048_R2.01.fits',field=1)
#
#coord_icrs=SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
#coord_gal=coord_icrs.galactic
#
#hpixels=hel.lonlat_to_healpix(coord_gal.l.degree*u.deg,coord_gal.b.degree*u.deg)
#values_40=hmap_40[hel.lonlat_to_healpix(coord_gal.l.degree*u.deg,coord_gal.b.degree*u.deg)]
#values_50=hmap_50[hel.lonlat_to_healpix(coord_gal.l.degree*u.deg,coord_gal.b.degree*u.deg)]
#
#galcut=np.zeros(len(lum))
#galcut40=np.zeros(len(lum))
#galcut50=np.zeros(len(lum))
#
#for i in range(len(lum)):
#        if values_40[i]>0.99:
#		galcut40[i]=1
#        if values_50[i]>0.99:
#                galcut50[i]=1
#for i in range(len(lum)):
#	if galcut40[i]==1 and galcut50[i]==0:
#                galcut[i]=1
#        if galcut40[i]==1 and galcut50[i]==1:
#                galcut[i]=2
#
#print len([x for x in galcut if x==0])
#print len([x for x in galcut if x==1])
#print len([x for x in galcut if x==2])
#######################


##############Perform div cuts #############
#divmap=enmap.read_fits(pathdiv)
#divmaps18=enmap.read_fits(pathdivs18)
#
#subs = 18 #in number of pixels
##map_rep =  np.zeros((len(rad),subs,subs)) # makes empty submap array
#divmap_rep =  np.zeros((len(ra),subs,subs))
#divmaps18_rep = np.zeros((len(ra),subs,subs))
##ps_rep=np.zeros((len(rad),subs,subs))
##ps2_rep=np.zeros((len(rad),subs,subs))
##ps3_rep=np.zeros((len(rad),subs,subs))
##array = coadded_data
#
#pixarray = enmap.sky2pix(divmap.shape,divmap.wcs,np.vstack((dec,ra))*np.pi/180.,safe=True)
#radpixels=pixarray[1,:]
#decdpixels=pixarray[0,:]
#
#print 'Creating submaps...'
#for l in range(len(lum)):
#        if (l == int(.1*len(ra))):
#            print '10% complete'
#        if (l == int(.2*len(ra))):
#            print '20% complete'
#        if (l == int(.3*len(ra))):
#            print '30% complete'
#        if (l == int(.4*len(ra))):
#            print '40% complete'
#        if (l == int(.5*len(ra))):
#            print '50% complete'
#        if (l == int(.6*len(ra))):
#            print '60% complete'
#        if (l == int(.7*len(ra))):
#            print '70% complete'
#        if (l == int(.8*len(ra))):
#            print '80% complete'
#        if (l == int(.9*len(ra))):
#            print '90% complete'
#        if (l == int(1.*len(ra))):
#            print 'done'
##        map_rep[l][:][:]=maps.cutout(enmapmap,arcmin_width=(subs-2.)*0.5,ra=None,dec=None,iy=decdpixels[l],ix=(radpixels[l]),pad=1)
#        divmap_rep[l][:][:]=maps.cutout(divmap[0][:][:],arcmin_width=(subs-2.)*0.5,ra=None,dec=None,iy=decdpixels[l],ix=(radpixels[l]),pad=1)
#        divmaps18_rep[l][:][:]=maps.cutout(divmaps18[0][:][:],arcmin_width=(subs-2.)*0.5,ra=None,dec=None,iy=decdpixels[l],ix=(radpixels[l]),pad=1)
#
#divcut=np.zeros(len(lum))
#divcuts18=np.zeros(len(lum))
#
#for i in range(len(lum)):
#        if 0.0<np.mean(divmap_rep[i][:][:]) and (1./65.0**2.)<=np.mean(divmap_rep[i][:][:]):
#                divcut[i]=1
#        if 0.0<np.mean(divmap_rep[i][:][:]) and (1./45.0**2.)<=np.mean(divmap_rep[i][:][:]):
#                divcut[i]=2
#        if 0.0<np.mean(divmaps18_rep[i][:][:]) and (1./65.0**2.)<=np.mean(divmaps18_rep[i][:][:]):
#                divcuts18[i]=1
#        if 0.0<np.mean(divmaps18_rep[i][:][:]) and (1./45.0**2.)<=np.mean(divmaps18_rep[i][:][:]):
#                divcuts18[i]=2
#
#print 'len divs', len(divcut)
#print '45uK:', len([x for x in divcut if x==2])
#print '65uK:', len([x for x in divcut if x==1])
#print 'none:', len([x for x in divcut if x==0])
#print 'len divs18', len(divcuts18)
#print 's18 45uK:', len([x for x in divcuts18 if x==2])
#print 's18 65uK:', len([x for x in divcuts18 if x==1])
#print 's18 none:', len([x for x in divcuts18 if x==0])
#
#######



###### Perform PS cuts 
#PS15mJy=enmap.read_fits(pathPS15mJy)
#PS100mJy=enmap.read_fits(pathPS100mJy)
#
#pixarray = enmap.sky2pix(PS15mJy.shape,PS15mJy.wcs,np.vstack((dec,ra))*np.pi/180.,safe=True)
#radpixels=pixarray[1,:]
#decdpixels=pixarray[0,:]
#
#subs = 18 #in number of pixels
#PS15mJy_rep=np.zeros((len(ra),subs,subs))
#PS100mJy_rep=np.zeros((len(ra),subs,subs))
#
#print 'Creating submaps...'
#for l in range(len(lum)):
#        if (l == int(.1*len(ra))):
#            print '10% complete'
#        if (l == int(.2*len(ra))):
#            print '20% complete'
#        if (l == int(.3*len(ra))):
#            print '30% complete'
#        if (l == int(.4*len(ra))):
#            print '40% complete'
#        if (l == int(.5*len(ra))):
#            print '50% complete'
#        if (l == int(.6*len(ra))):
#            print '60% complete'
#        if (l == int(.7*len(ra))):
#            print '70% complete'
#        if (l == int(.8*len(ra))):
#            print '80% complete'
#        if (l == int(.9*len(ra))):
#            print '90% complete'
#        if (l == int(1.*len(ra))):
#            print 'done'
#        PS15mJy_rep[l][:][:]=maps.cutout(PS15mJy,arcmin_width=(subs-2.)*0.5,ra=None,dec=None,iy=decdpixels[l],ix=(radpixels[l]),pad=1)
#        PS100mJy_rep[l][:][:]=maps.cutout(PS100mJy,arcmin_width=(subs-2.)*0.5,ra=None,dec=None,iy=decdpixels[l],ix=(radpixels[l]),pad=1)
#
#PS15mJy_cut=np.zeros(len(lum))
#PS100mJy_cut=np.zeros(len(lum))
#
#for i in range(len(lum)):
#        if PS15mJy_rep[i][:][:].all()==True:
#                PS15mJy_cut[i]=1
#        if PS100mJy_rep[i][:][:].all()==True:
#                PS100mJy_cut[i]=1
#
#
#print 'len PS15', len(PS15mJy)
#print 'Pass:', len([x for x in PS15mJy_cut if x==1])
#print 'Fail:', len([x for x in PS15mJy_cut if x==0])
#print 'len PS100', len(PS100mJy)
#print 'Pass:', len([x for x in PS100mJy_cut if x==1])
#print 'Fail:', len([x for x in PS100mJy_cut if x==0])



######## Mask Component Separated Map Areas (D56 and BN)
#
#d56mask=enmap.read_fits(pathd56mask, hdu=None, sel=None, sel_threshold=10e6)
#bnmask=enmap.read_fits(pathbnmask, hdu=None, sel=None, sel_threshold=10e6)
#
#subs = 18 #in number of pixels
#
#pixarrayd56 = enmap.sky2pix(d56mask.shape,d56mask.wcs,np.vstack((dec,ra))*np.pi/180.,safe=True)
#radpixelsd56=pixarrayd56[1,:]
#decdpixelsd56=pixarrayd56[0,:]
#
#pixarraybn = enmap.sky2pix(bnmask.shape,bnmask.wcs,np.vstack((dec,ra))*np.pi/180.,safe=True)
#radpixelsbn=pixarraybn[1,:]
#decdpixelsbn=pixarraybn[0,:]
#
#d56mask_rep=np.zeros((len(ra),subs,subs))
#bnmask_rep=np.zeros((len(ra),subs,subs))
#
#
#print 'Creating submaps...'
#for l in range(len(lum)):
#        if (l == int(.1*len(ra))):
#            print '10% complete'
#        if (l == int(.2*len(ra))):
#            print '20% complete'
#        if (l == int(.3*len(ra))):
#            print '30% complete'
#        if (l == int(.4*len(ra))):
#            print '40% complete'
#        if (l == int(.5*len(ra))):
#            print '50% complete'
#        if (l == int(.6*len(ra))):
#            print '60% complete'
#        if (l == int(.7*len(ra))):
#            print '70% complete'
#        if (l == int(.8*len(ra))):
#            print '80% complete'
#        if (l == int(.9*len(ra))):
#            print '90% complete'
#        if (l == int(1.*len(ra))):
#            print 'done'
#	if (0.0<ra[l]<45.3 or 348.7<ra[l]<360.0) and (-9.1<dec[l]<5.8):
#                d56mask_rep[l][:][:]=maps.cutout(d56mask[:][:],arcmin_width=(subs)*0.5,ra=None,dec=None,iy=decdpixelsd56[l],ix=(radpixelsd56[l]),pad=1)
#        if (146.0<ra[l]<247.0) and (-5.9<dec[l]<23):
#                bnmask_rep[l][:][:]=maps.cutout(bnmask[:][:],arcmin_width=(subs)*0.5,ra=None,dec=None,iy=decdpixelsbn[l],ix=(radpixelsbn[l]),pad=1)
#	
#incompsepregion=np.zeros(len(lum))
#compsepregionvalsd56=np.zeros(len(lum))
#compsepregionvalsbn=np.zeros(len(lum))
#
#
#for i in range(len(lum)):
#        compsepregionvalsd56[i]=np.mean(d56mask_rep[i][:][:])
#        compsepregionvalsbn[i]=np.mean(bnmask_rep[i][:][:])
#        if np.mean(d56mask_rep[i][:][:])>0.9:
#                incompsepregion[i]=1
#        if np.mean(bnmask_rep[i][:][:])>0.9:
#                incompsepregion[i]=2
#
#print 'Total',len(incompsepregion)
#print 'In D56',len([x for x in incompsepregion if x==1])
#print 'In BN',len([x for x in incompsepregion if x==2])


####### Define overall analysis flags ####

#S18coadd=np.zeros(len(lum))
#S16ILC=np.zeros(len(lum))
#
#for i in range(len(z)):
#	if PS15mJy_cut[i]==1 and PS100mJy_cut[i]==1 and divcuts18[i]==2 and galcut[i]==2 and use_PS[i]==True:
#		S18coadd[i]=1
#	if PS15mJy_cut[i]==1 and PS100mJy_cut[i]==1 and divcut[i]==2 and galcut[i]==2 and (incompsepregion[i]==1 or incompsepregion[i]==2) and use_PS[i]==True:
#		S16ILC[i]=1
#
#print 'S18 Coadd=',len([x for x in S18coadd if x==1])
#print 'S16 ILC=',len([x for x in S16ILC if x==1])

np.savetxt('V20_DR15_Catalog_v2.csv',np.array([ra,dec,lum,z,stellarmasses,halomasses,cModelMag_u,cModelMag_g,cModelMag_r,cModelMag_i,cModelMag_z,cModelMagErr_u,cModelMagErr_g,cModelMagErr_r,cModelMagErr_i,cModelMagErr_z,petroMag_u,petroMag_g,petroMag_r,petroMag_i,petroMag_z,petroMagErr_u,petroMagErr_g,petroMagErr_r,petroMagErr_i,petroMagErr_z,extinction_u,extinction_g,extinction_r,extinction_i,extinction_z,bestObjID,PS15mJy_cut,PS100mJy_cut,divcut,divcuts18,galcut,incompsepregion,use_PS,S18coadd,S16ILC]).T,delimiter=',',fmt="%s") 








