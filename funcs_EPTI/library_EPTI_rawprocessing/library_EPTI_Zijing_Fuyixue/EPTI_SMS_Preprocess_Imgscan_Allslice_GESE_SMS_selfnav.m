function  [ K,K_nav,linear_fit_coeff ] = EPTI_SMS_Preprocess_Imgscan_Allslice_GESE_SMS_selfnav(filename,RepsToRead,OverlapGroupToUnfold,linear_fit_coeff_all,seq_type)
% EPTI data raw data processing
% GE-SE data with different number of echoes processed
% raw data processing for GESE EPTI sequence 
% Fuyixue Wang - Jan 2018

% make it work for ACE-EPTI sequence 
% Zijing Dong - Jan 2021

% single-shot EPTI sequence with SMS
% Fuyixue Wang- Nov 2021

% single-shot EPTI sequence with self phase cor navigator
% Fuyixue Wang- Dec 2022

% generalized for single-shot and multi-shot EPTI sequence processing for
% 3T and 7T data for C2P
% Fuyixue Wang & Zijing Dong, updated 2023
%% parameters
% nRepToRead                                  % number of repetition to read
% BeginRep                                    % the first rep position
% OverlapGroupToUnfold                        % Group to unfold
% RepsToRead = [1:1:nRepToRead];              % non consecutive reps that you want to read
%%
% data to reconstruct
subgroupUnfolding = 0;                          % set to 1 if want to only do unfolding of the group specified by 'OverlapGroupToUnfold'...
ascendingSlice_acq = 0;                         % set to 1 for non-interleave data
%% load data 

for count = 1:length(RepsToRead)
    nRepToRead = 1;
    BeginRep = RepsToRead(count);
    [meas_Current] = read_meas_dat_memmap_EPTI_SMS_Generalized(filename,nRepToRead,BeginRep,ascendingSlice_acq);
    if count == 1
        meas = meas_Current;
    else
        meas.data = cat(7,meas.data, meas_Current.data);
        if isfield(meas,'data_phascor1d') 
            meas.data_phascor1d = cat(7,meas.data_phascor1d, meas_Current.data_phascor1d);
        end
    end
end
SMS_data = meas.prot.SMSdata;

GC_apply_flag = 1-isempty(linear_fit_coeff_all);
if (meas.prot.navigatorOn == 0) && (seq_type==0)
    meas.data_phascor1d = meas.data(:,1:3,:,:,:,:,:,:,:,:,:,:); % first 3 lines are self phasecor navigator
end

%% separate SMS data from reference 
if (SMS_data)
    meas_Collapsed_all.data = meas.data;
    meas_Collapsed_all.data_phascor1d = meas.data_phascor1d;
    meas_Collapsed_all.prot = meas.prot; 
    meas_Collapsed_all.evp = meas.evp; 

    NslicesEX = meas.prot.sWipMemBlock_adFree(1); %number of overlapping slices
    if isempty(meas_Collapsed_all.evp.NFirstLin)
        meas_Collapsed_all.evp.NFirstLin = 1;
    end
    meas_Collapsed = meas_Collapsed_all;
else
    meas_Indiv_all = meas;  
    NslicesEX = 1;
end

if (mod(NslicesEX,2) == 0)
  meas_Collapsed.data = cat(10,meas_Collapsed.data(:,:,:,:,:,:,:,:,:,end),meas_Collapsed.data(:,:,:,:,:,:,:,:,:,1:end-1));
  if seq_type == 0 || meas.prot.navigatorOn == 1
      meas_Collapsed.data_phascor1d = cat(10,meas_Collapsed.data_phascor1d(:,:,:,:,:,:,:,:,:,end),meas_Collapsed.data_phascor1d(:,:,:,:,:,:,:,:,:,1:end-1));
  end
end

meas.data=[];
meas.data_phascor1d=[];
meas.patrefscan=[];
meas.patrefscan_phascor=[];
meas.smsrefscan=[];
meas.smsrefscan_phascor=[];
clear meas_Collapsed_all;
%% figuring out parameters from data
if (SMS_data) % SMS parameters
    if subgroupUnfolding ~=1
        OverlapGroupToUnfold = 1:size(meas_Collapsed.data,10);
    end
    NslicesEX = meas.prot.sWipMemBlock_adFree(1); %number of overlapping slices
    if meas.prot.sWipMemBlock_adFree(3) == 1
        PhaseShiftBtwSimulSlices = 0;
    else
        PhaseShiftBtwSimulSlices = 2*pi/meas.prot.sWipMemBlock_adFree(3);
        % phase shift between undersampled ky lines (i.e. FOV/4 with IPAT2 => pi,FOV/6 with IPAT3 => pi)
    end

    if sum(meas_Collapsed.data(:,1,1,1,1,1,1,1,1,1,1,1),1) ~= 0
        DataCollapsedSegOne_OddLine = 1;
    else
        DataCollapsedSegOne_OddLine = 0;
    end
    nEchos = mrir_ice_dimensions(meas_Collapsed.data, 'eco');
    SliceSep=meas.prot.SliceSep;
    meas_Collapsed_all = meas_Collapsed;

else % no SMS 
    if subgroupUnfolding ~=1
        OverlapGroupToUnfold = 1:size(meas_Indiv_all.data,10);
    end
    nEchos = mrir_ice_dimensions(meas_Indiv_all.data, 'eco');
    meas_Indiv = meas_Indiv_all;
end

%% Grappa Recon Loop for each Overlap group
for OverlapGroupCount = 1:length(OverlapGroupToUnfold)
    %% Ghost correct and grid, initial 
%     disp('GhostCorrect and Grid ')
    slice_index=OverlapGroupToUnfold(OverlapGroupCount);
    if SMS_data
        meas_Collapsed.data=meas_Collapsed_all.data(:,:,:,:,:,:,:,:,:,slice_index);
        if (seq_type==0) || (meas.prot.navigatorOn == 1)
            meas_Collapsed.data_phascor1d=meas_Collapsed_all.data_phascor1d(:,:,:,:,:,:,:,:,:,slice_index);   
        end
    else
        meas_Indiv.data=meas_Indiv_all.data(:,:,:,:,:,:,:,:,:,slice_index);
        if (seq_type==0) || (meas.prot.navigatorOn == 1)
            meas_Indiv.data_phascor1d=meas_Indiv_all.data_phascor1d(:,:,:,:,:,:,:,:,:,slice_index);
        end
    end
    if GC_apply_flag == 1
        linear_fit_coeff = linear_fit_coeff_all(:,:,slice_index);
%         if SMS_data
%             linear_fit_coeff = repmat(linear_fit_coeff,[1 1 size(meas_Collapsed.data,5)]); % multi echo
%         else
%             linear_fit_coeff = repmat(linear_fit_coeff,[1 1 size(meas_Indiv.data,5)]); % multi echo
%         end

        if SMS_data
            [meas_Collapsed] = GhostCorrectAndGrid_mkm_PSFv2(meas_Collapsed,linear_fit_coeff,0);
        else 
            [meas_Indiv] = GhostCorrectAndGrid_mkm_PSFv2(meas_Indiv,linear_fit_coeff,0);  % no meas_Collapsed
        end
    else
        if SMS_data
            [meas_Collapsed,linear_fit_coeff(:,:,slice_index)] = GhostCorrectAndGrid_mkm_PSFv2(meas_Collapsed,[],1);
        else 
            [meas_Indiv,linear_fit_coeff(:,:,slice_index)] = GhostCorrectAndGrid_mkm_PSFv2(meas_Indiv,[],1);  % no meas_Collapsed
        end
    end
    if (SMS_data) 
        meas_Collapsed.data = sum(meas_Collapsed.data,8); 
        if meas.prot.navigatorOn == 1 
            meas_Collapsed.data_phascor1d = sum(meas_Collapsed.data_phascor1d,8);
            K_nav(:,:,:,:,:,slice_index) = squeeze(meas_Collapsed.data_phascor1d);
        end
    else
        meas_Indiv.data = sum(meas_Indiv.data,8); 
        if meas.prot.navigatorOn == 1 
            meas_Indiv.data_phascor1d = sum(meas_Indiv.data_phascor1d,8);
            K_nav(:,:,:,:,:,slice_index) = squeeze(meas_Indiv.data_phascor1d);
        end
    end
    %% Caipirinha
    if (SMS_data)
        if (meas.prot.navigatorOn == 0) && (seq_type==0)
            start_point = 3;
        else
            start_point = 1;
        end
%         start_point = 1;
        K_CollapsedFull = CaipirinhaDeblur_SingleSlice_EPTI(meas_Collapsed.data, meas_Collapsed.prot, meas_Collapsed.evp , PhaseShiftBtwSimulSlices, SliceSep, slice_index,start_point);
    else
        K_CollapsedFull = meas_Indiv.data;
    end
    K_all(:,:,:,:,:,slice_index) = squeeze(K_CollapsedFull);
end  %end OverlapGroupCount      

K=K_all;
k_size=size(K); k_size(1)=k_size(1)/2;
K=fftc((crop(ifftc(K,1),k_size)),1);
if meas.prot.pf_echo == 0
    K=permute(K,[2 1 4 3 6 5]);
else
    K=permute(K,[2 1 5 3 6 4]);
end
if meas.prot.navigatorOn == 0 
    K_nav = K(1:3,:,:,:,:,1);
else
    k_size=size(K_nav); k_size(1)=k_size(1)/2;
    K_nav=fftc((crop(ifftc(K_nav,1),k_size)),1);
    K_nav=permute(K_nav,[2 1 4 3 6 5]);
end

end
