function  [ K,K_nav,linear_fit_coeff] = EPTI_SMS_Preprocess_Imgscan_Allslice_Calib(filename,SelectiveReps,RepsToRead,SMS_data,OverlapGroupToUnfold,pf_echo  )
% EPTI data raw data processing, all slices
% Gradient Echo calibration data with different number of echoes processed
% Fuyixue Wang, Jan 2018 
%% parameters
% nRepToRead                                % number of repetition to read
% BeginRep                                  % the first rep position
% OverlapGroupToUnfold                      % Group to unfold
% SelectiveReps = 1;                        % if want to read non cosecutive reps then turn this flag on
% RepsToRead = [1:1:nRepToRead];            % non consecutive reps that you want to read
%%
subgroupUnfolding = 0;                          % set to 1 if want to only do unfolding of the group specified by 'OverlapGroupToUnfold'...
ascendingSlice_acq = 0;                         % set to 1 for non-interleave data
%% load data 
if SelectiveReps == 0
   nRepsToRead=size(RepsToRead,2);
   [meas] = read_meas_dat_memmap_EPTI_GESE(filename,nRepsToRead,1,SMS_data,ascendingSlice_acq,pf_echo);
else      
    for count = 1:length(RepsToRead)
        nRepToRead = 1;
        BeginRep = RepsToRead(count);
        [meas_Current] = read_meas_dat_memmap_EPTI_GESE_selective(filename,nRepToRead,BeginRep,SMS_data,ascendingSlice_acq,pf_echo);
        if count ==1
            meas = meas_Current;
        else
            meas.data = cat(7,meas.data, meas_Current.data);
            meas.data_phascor1d = cat(7,meas.data_phascor1d, meas_Current.data_phascor1d);
        end
    end
end

%% separate SMS data from reference 
meas_Indiv_all = meas;  
meas_Indiv = meas_Indiv_all;
if subgroupUnfolding ~=1
    OverlapGroupToUnfold = 1:size(meas_Indiv_all.data,10);
end
%% Grappa Recon Loop for each Overlap group
for OverlapGroupCount = 1:length(OverlapGroupToUnfold)
    slice_index=OverlapGroupToUnfold(OverlapGroupCount);
    meas_Indiv.data=meas_Indiv_all.data(:,:,:,:,:,:,:,:,:,slice_index);
    meas_Indiv.data_phascor1d=meas_Indiv_all.data_phascor1d(:,:,:,:,:,:,:,:,:,slice_index);
    [meas_Indiv, linear_fit_coeff] = GhostCorrectAndGrid_mkm_PSFv2(meas_Indiv,[],1);  % no meas_Collapsed
    meas_Indiv.data = sum(meas_Indiv.data,8);
    meas_Indiv.data_phascor1d = sum(meas_Indiv.data_phascor1d,8);

    K_all(:,:,:,:,:,slice_index) = squeeze(meas_Indiv.data);
    K_nav(:,:,:,:,:,slice_index) = squeeze(meas_Indiv.data_phascor1d);
end  %end OverlapGroupCount      

K = K_all;
k_size=size(K); k_size(1)=k_size(1)/2;
K=fftc((crop(ifftc(K,1),k_size)),1);
K=permute(K,[2 1 4 3 5 6]);

k_size=size(K_nav); k_size(1)=k_size(1)/2;
K_nav=fftc((crop(ifftc(K_nav,1),k_size)),1);
K_nav=permute(K_nav,[2 1 4 3 5 6]);
end
%%



