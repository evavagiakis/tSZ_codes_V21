V20 analysis pipeline: 

1) CatalogCreationCode.py is how the DR15 catalog is trimmed and flagged for use
2) Submaps from the analyzed maps are extracted using the catalog produced by 1)
3) tSZ_stacking_V20.py is used to stack these submaps for the full catalog sample and save raselect,decselect,lumselect,zselect,disks,disk_stds,rings,ring_stds,divs for the selected source sample via aperture photometry 
4) tSZ_AP_withjackknife_andDivWeighting_V20.py (for the coadded maps) or tSZ_AP_withjackknife_andDivWeighting_component_sep_mask_V20.py (for the ILC maps) takes the output from 3) and computes dT or y and JK error bars per bin, saving them into a .txt file 
5) DustAndBeamCorrections.py takes the output from 4) and corrects it via Herschel dust corrections, beam corrections, then calculates tau and tau uncertainties, and compares these to theory tau 

