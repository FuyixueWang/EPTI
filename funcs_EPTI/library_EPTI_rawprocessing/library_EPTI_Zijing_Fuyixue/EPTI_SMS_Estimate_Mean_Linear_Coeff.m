function  [linear_fit_coeff ] = EPTI_SMS_Estimate_Mean_Linear_Coeff(filename,SelectiveReps,RepsToRead,SMS_data,OverlapGroupToUnfold,pf_echo,SMS_factor,PhasCorSegSwap)
% Estimate linear coefficient for eddy-current related phase correction
% between even and odd echoes from calibration data, and can later apply
% to EPTI data
% Fuyixue Wang & Zijing Dong, 2022
%% parameters
% data to reconstruct
subgroupUnfolding = 1;                          % set to 1 if want to only do unfolding of the group specified by 'OverlapGroupToUnfold'...
ascendingSlice_acq = 0;                         % set to 1 for non-interleave data

%% load data 
if SelectiveReps == 0
   nRepsToRead=size(RepsToRead,2);
   [meas] = read_meas_dat_memmap_EPTI_GESE(filename,nRepsToRead,1,SMS_data,ascendingSlice_acq,pf_echo,PhasCorSegSwap);
else      
    for count = 1:length(RepsToRead)
        nRepToRead = 1;
        BeginRep = RepsToRead(count);
        [meas_Current] = read_meas_dat_memmap_EPTI_GESE_selective(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo,PhasCorSegSwap);
        if count ==1
            meas = meas_Current;
        else
            meas.data = [];
            meas.data_phascor1d = cat(7,meas.data_phascor1d, meas_Current.data_phascor1d);
        end
    end
end

%% separate SMS data from reference 
meas_Indiv = meas;  
NslicesEX = 1;
meas.data=[];
meas.data_phascor1d=[];
meas.patrefscan=[];
meas.patrefscan_phascor=[];
meas.smsrefscan=[];
meas.smsrefscan_phascor=[];
clear meas_Collapsed_all;
N_slice = meas.prot.sSliceArray_lSize;
slice_sep  =  N_slice/SMS_factor;
%% figuring out parameters from data
if subgroupUnfolding ~=1
    OverlapGroupToUnfold = 1:size(meas_Indiv.data,10);
end
%% Grappa Recon Loop for each Overlap group
% for OverlapGroupCount = 1:length(OverlapGroupToUnfold)
for OverlapGroupCount = 1:slice_sep
%     disp('*********************************************************************')
    %% Ghost correct and grid, initial 
%     disp('*********************************************************************')
%     disp('GhostCorrect and Grid ')
    tic
    slice_index=OverlapGroupToUnfold(OverlapGroupCount);
    meas = meas_Indiv;
    meas.data_phascor1d=mean(meas.data_phascor1d(:,:,:,:,:,:,:,:,:,slice_index:slice_sep:end),10);
%     meas.data_phascor1d=meas.data_phascor1d(:,:,:,:,:,:,:,:,:,slice_index);

    dims = size(meas.data_phascor1d);
    if length(dims) < 10
        dims(end+1:10) = 1; 
    end
    linear_fit_coeff_tmp = mrir_artifact_ghost_compute(meas.data_phascor1d);
    linear_fit_coeff_tmp = reshape( median( reshape(linear_fit_coeff_tmp, [2 1 dims(3:7) 1 1 dims(10)]), 7) , 2,[]) ; % for saving to use for tailor ghost correction   
    linear_fit_coeff(:,:,slice_index) = linear_fit_coeff_tmp;
end  %end OverlapGroupCount      
end




