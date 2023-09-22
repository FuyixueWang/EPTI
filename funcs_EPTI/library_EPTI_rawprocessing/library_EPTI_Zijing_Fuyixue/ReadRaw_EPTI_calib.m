function [kdata_GE,parameters,linear_fit_coeff,K_nav,meas] = ReadRaw_EPTI_calib(filename,output_ghost_coeff,SMS_factor,PhasCorSegSwap)
% EPTI rawdata processing for calibraiton data
% Fuyixue Wang, MGH 2021

% optimized for C2P, adding support to estimate, apply, and output 
% eddy-current related even-odd echo phase
% Fuyixue Wang, Zijing Dong, MGH, 2023

%%
if nargin < 2
    output_ghost_coeff = 0;
end
if nargin < 3
    SMS_factor = [];
end
if nargin < 4
    PhasCorSegSwap = 0;
end
nRepToRead=1; BeginRep=1; SMS_data=0; ascendingSlice_acq=0; pf_echo=0; 
[meas] = read_meas_dat_memmap_EPTI_GESE(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo);
Nreps = meas.prot.lRepetitions+1;
parameters=meas.prot;
%% pre estimation of ghost correction coefficients (SMS slice group average)
RepsToRead=1:Nreps;
if output_ghost_coeff == 1
    slice=1:1:parameters.sSliceArray_lSize;
    slice_sep = parameters.sSliceArray_lSize/SMS_factor;
    [linear_fit_coeff] = EPTI_SMS_Estimate_Mean_Linear_Coeff(filename,1,RepsToRead,SMS_data,slice,pf_echo,SMS_factor,PhasCorSegSwap);
    [kdata_GE,K_nav] = EPTI_SMS_Preprocess_Imgscan_Allslice_Calib_GC(filename,0,RepsToRead,SMS_data,[],pf_echo,linear_fit_coeff,slice_sep);
else
    linear_fit_coeff=[];
    [kdata_GE,K_nav]  = EPTI_SMS_Preprocess_Imgscan_Allslice_Calib(filename,0,RepsToRead,SMS_data,[],pf_echo);
end

kdata_GE=single(kdata_GE);
K_nav=single(K_nav);

parameters.nechoGE=size(kdata_GE,1);

meas.data =[];
meas.data_phascor1d = []; 
end
