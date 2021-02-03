V20 analysis pipeline: 

1) CatalogCreationCode.py is how the DR15 catalog is trimmed and flagged for use
	Takes: SDSS Casjobs query, kcorrected luminosities, inverse variance white noise maps, PS masks, ILC area masks and/or old catalog in same format; filepaths are hard coded in script
	Returns: .csv catalog file
	To run: >>> python CatalogCreationCode.py

2) Submaps from the analyzed maps are extracted using the catalog produced by 1) and codes by Patricio Gallardo. 
	The script to call to generate submaps is: https://github.com/patogallardo/iskay/blob/master/misc/qsub_exportSubmaps.sh
	Which will call the following script for each map: https://github.com/patogallardo/iskay/blob/master/misc/iskay_exportSubmaps.py
	From here, this script will use pixell and iskay (the pairwise code) to generate the submaps

3) tSZ_stacking_V20.py is used to stack these submaps for the full catalog sample and save raselect,decselect,lumselect,zselect,disks,disk_stds,rings,ring_stds,divs for the selected source sample via aperture photometry 
	Takes: .csv catalog via 1), submaps extracted in 2), option of which map you are analyzing. Submap paths are hard coded in script, catalog and map analysis option are taken as arguments 
	Returns: .txt file of raw aperture photometry values for entire catalog 
	To run (e.g. for DR5 f150 analysis): >>> python tSZ_stacking_V20.py -c ~/samba/V20_DR15_Catalog_v3.csv -a DR5f150

4) tSZ_AP_withjackknife_andDivWeighting_V20.py (for the coadded maps) or tSZ_AP_withjackknife_andDivWeighting_component_sep_mask_V20.py (for the ILC maps) takes the output from 3) and computes dT or y and JK error bars per bin, saving them into a .txt file 
	Takes: Output .txt file from 3) (two in the case of the ILC maps -- one for each area, BN and D56). Paths are hard coded in script
	Returns: Prints out and/or saves .txt file of aperture photometry results for each luminosity bin
	To run: >>> python tSZ_AP_withjackknife_andDivWeighting_V20.py

5) DustAndBeamCorrections.py takes the output from 4) and corrects it via Herschel dust corrections, beam corrections, then calculates tau and tau uncertainties, and compares these to theory tau 
	Takes: Output .txt files from 4) for the 3 maps analyzed, dust correction .txt files, both hard coded in script
	Returns: Prints out results (can comment out/uncomment which sections you want printed in script)



