# EPTI
EPTI (Echo Planar Time-resolved Imaging) raw-data processing, image reconstruction and data preprocessing.
Reconstruction version date: 09/2023; 

We will continue to optimize the reconstruction and the sequence. Updated reconstruction code, and available sequences will be posted on our website https://martinos.org/~fw089/.

Fuyixue Wang <fwang18@mgh.harvard.edu>; Zijing Dong <zdong3@mgh.harvard.edu>;  July 2023,  MGH

-------------------------------------------------------------------------------------------------------------
The MATLAB scripts and functions for EPTI raw-data processing, image reconstruction and data preprocessing.
It supports:
- 2D EPTI sequence
- gradient echo(GE) and/or asymmetric spin-echo(ASE)/spin-echo(SE) EPTI;
- single-shot and multi-shot EPTI;
- EPTI with or without SMS acquisition;
- 3T and 7T EPTI.

The main MATLAB scripts included in the folders 'main_scripts_VE11C_3T' and 'main_scripts_VE12U_7T' were prepared with preset parameters for the example EPTI protocols released in the EPTI C2P sequence package. (You can request the EPTI sequence for Siemens scanners through MGH C2P Program http://nmr.mgh.harvard.edu/c2p or Siemens teamplay C2P Exchange).

The scripts produce reconstructed EPTI images from Twix rawdata, and saves NIFTI output and some intermediate MATLAB files (see details about the output below).

All MATLAB functions needed are included in the 'funcs_EPTI' folder.

The main MATLAB scripts will call some BART functions, so BART installation is required, please download bart-0.8.00, copy the EPTI customized files included in 'bart-0.8.00_EPTI/' to the BART folder, and overwrite the files including the 'Makefile' and source files in 'src/'.

All sequence and reconstruction parameters were optimized for in-vivo human brain data. Phantom data may need different set of parameters.


Please cite the following work on EPTI for this, including: 

1. Wang F, Dong Z, Reese TG, et al. Echo planar time-resolved imaging (EPTI). Magn Reson Med. 2019 Jun;81(6):3599-3615. doi: 10.1002/mrm.27673.

2. Dong Z, Wang F, Reese TG, et al. Echo planar time-resolved imaging with subspace reconstruction and optimized spatiotemporal encoding. Magn Reson Med. 2020 Nov;84(5):2442-2455. doi: 10.1002/mrm.28295.

3. Dong Z, Wald LL, Polimeni JR, Wang F. Single-shot echo planar time-resolved imaging for multi-echo functional MRI and distortion-free diffusion imaging. Magn Reson Med. 2024 Oct 20. doi: 10.1002/mrm.30327.

Other relevant citations that may warrant consideration (EPTI's application in fMRI, dMRI and qMRI), see EPTI website: https://martinos.org/~fw089/


-----------------------------------------------------------------------------------------------------
The main output includes (note that no motion correction was performed across dynamics in the output):
- EPTIdata/2_Recon_nii/*filename/*filename_im_sos.nii: final output time-series images (all echoes combined via sos), [nx,ny,nz,ndynamics]; (check this data first for image quality, GESE acquisition has separated sos images)
- EPTIdata/2_Recon_nii/*filename/*filename_img_Tave_echoes.nii: time-resolved echoes (average all dynamics), [nx,ny,nz,nechoes];
- EPTIdata/2_Recon_nii/*filename/echoes/*filename_im_echoN.nii: time-series images of the Nth echo, [nx,ny,nz,ndynamics];
- EPTIdata/2_Recon_nii/*filename/*filename__TEs.txt: TE values of all the echoes.

Additional output if it is a gradient echo acquisition:
- EPTIdata/2_Recon_nii/*filename/*filename_im_comb.nii: final output time-series images after T2*-weighted optimal echo-combination, [nx,ny,nz,ndynamics]; (optimized CNR for fMRI)
- EPTIdata/2_Recon_nii/*filename/*filename__T2s.nii: fitted T2* map, [nx,ny,nz];

Additional output if it is a gradient-echo & spin-echo (GESE) acquisition:
- EPTIdata/2_Recon_nii/*filename/*filename__T2s.nii: fitted T2* map; [nx, ny, nz]
- EPTIdata/2_Recon_nii/*filename/*filename__T2.nii: fitted T2 map; [nx, ny, nz]

Acquisition parameters
- EPTIdata/0_Data_acq/meas_prot_*filename.mat: acquisition parameters for imaging scan;
- EPTIdata/0_Data_acq/meas_prot_Calib_*filename.mat: acquisition parameters for calibration scan;

Other intermediate data for debugging
- EPTIdata/1_Recon_mat/*filename/Recon_EPTI_*filename_Dyn_*N.mat: reconstructed complex image for each dynamic saved in MATLAB format;
- EPTIdata/0_Data_acq/Basis_*filename.mat: subspace basis used for reconstruction;
- EPTIdata/0_Data_acq/PreRecon_SensB0_*filename.mat: calibrated coil sensitivity, B0 and phase information;

