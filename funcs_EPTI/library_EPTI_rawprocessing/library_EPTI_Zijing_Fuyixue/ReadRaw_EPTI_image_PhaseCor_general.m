function [kdata,k_nav,parameters,linear_fit_coeff] = ReadRaw_EPTI_image_PhaseCor_general(filename,parameters,dyn,process_type,seq_type,gccmtx_aligned,linear_fit_coeff,dt,k_nav_first,TEs_all)
% EPTI rawdata processing for imaging data
% including reading rawdata, raw data processing, and put into EPTI 
% encoding pattern
% Fuyixue Wang, 2021, MGH

% Generalzed code for both single-shot and multi-shot EPTI data
% support applying odd-even echo eddy-current phase correction
% support to both GE and SE EPTI data with SMS
% Fuyixue Wang, Zijing Dong, MGH, 2023

% added phase correction and shot-to-shot correction optimized to reduce 
% physiological noise/artifacts for higher temporal stability
% Fuyixue Wang, Zijing Dong, MGH, Aug, 2023
%%
if nargin < 6
    gccmtx_aligned=[];
end
if nargin < 7
    linear_fit_coeff=[];
end

RepsToRead=1:parameters.Nseg;
RepsToRead=RepsToRead+(dyn-1)*parameters.Nseg;
[kdata_all,k_nav,linear_fit_coeff]  = EPTI_SMS_Preprocess_Imgscan_Allslice_GESE_SMS_selfnav(filename,RepsToRead,[],linear_fit_coeff,seq_type);

switch process_type
    case 0
        if parameters.pf_echo > 0
            kdata_all = kdata_all(1:floor(parameters.pf_echo*size(kdata_all,1)),:,:,:,:,1);
        else
            kdata_all = kdata_all(:,:,:,:,:,1);
        end
    case 1  % 0: SE only
        kdata_all=kdata_all(1:end,:,:,:,:,2);
end

if seq_type == 1
    B0_drift_Hz_all = DriftB0_esti_nav_ms_coil_aveSlice2group(k_nav_first,k_nav,dt); % field drift/global B0 change estimation
    B0_drift_Hz = B0_drift_Hz_all(:,:,:,1);
    correct_phase = B0_drift_Hz.*permute(TEs_all(:),[2 3 4 1])*2*pi;
    kdata_all(:,:,:,:,1:2:end) = kdata_all(:,:,:,:,1:2:end).*permute(repmat(exp(-1i*correct_phase),[1,1,ceil(size(kdata_all,5)/2),1,size(kdata_all,2)]),[4 5 1 2 3]);
    B0_drift_Hz = B0_drift_Hz_all(:,:,:,2);
    correct_phase = B0_drift_Hz.*permute(TEs_all(:),[2 3 4 1])*2*pi;
    kdata_all(:,:,:,:,2:2:end) = kdata_all(:,:,:,:,2:2:end).*permute(repmat(exp(-1i*correct_phase),[1,1,floor(size(kdata_all,5)/2),1,size(kdata_all,2)]),[4 5 1 2 3]);
else
    B0_drift_Hz = DriftB0_esti_nav_ms(k_nav_first,k_nav,dt); % field drift/global B0 change estimation
    correct_phase = B0_drift_Hz.*permute(TEs_all(:),[2 3 4 1])*2*pi;
    kdata_all = kdata_all.*permute(repmat(exp(-1i*correct_phase),[1,1,1,1,size(kdata_all,2)]),[4 5 3 2 1]);
end

if (isempty(gccmtx_aligned)==0)
    s = size(kdata_all);
    s(4) = size(gccmtx_aligned,2);
    kdata_cc=zeros(s);
    for t=1:size(kdata_all,1)
        for echo = 1:size(kdata_all,6)
            CCDATA_aligned = CC(permute(kdata_all(t,:,:,:,:,echo),[2 3 5 4 1 6]),gccmtx_aligned, 1);
            kdata_cc(t,:,:,:,:,echo)=permute(CCDATA_aligned,[1 2 4 3]);
        end
    end
    kdata_all=kdata_cc;
    clear kdata_cc;
end

switch process_type
    case 0  % 0: GE only
        if seq_type == 0 % zigzag/single-shot version
            kdata = putRawtoPattern_SMS_zigzag_PF_selfnav(kdata_all,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift,parameters.PF_shift);
        else             % multi-shot version
            kdata = putRawtoPattern_EPTI_SMS(kdata_all,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift);
            kdata = kdata(:,:,end:-1:1,:,:,:);
        end
        
    case 1  % 0: SE only
        if seq_type == 0 % zigzag/single-shot version
            kdata = putRawtoPattern_SMS_zigzag_PF_selfnav(kdata_all,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift,parameters.PF_shift);
        else             % multi-shot version
            kdata = putRawtoPattern_EPTI_SMS(kdata_all,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift);
            kdata = kdata(:,:,end:-1:1,:,:,:);
        end
%     case 2 %  GE-SE joint
%         kdata_SE=kdata_all(1:end,:,:,:,:,2);
%         kdata_GE=kdata_all(1:floor(parameters.pf_echo*size(kdata_all,1)),:,:,:,:,1);
% 
%         if seq_type == 0 % zigzag/single-shot version
%             kdata_GE = putRawtoPattern_SMS_zigzag_PF_selfnav(kdata_GE,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift,parameters.PF_shift);
%             kdata_SE = putRawtoPattern_SMS_zigzag_PF_selfnav(kdata_SE,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift,parameters.PF_shift);
%         else             % multi-shot version
%             kdata_GE = putRawtoPattern_EPTI_SMS(kdata_GE,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift);
%             kdata_SE = putRawtoPattern_EPTI_SMS(kdata_SE,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift);
%         end
%         kdata = cat(1,kdata_GE,kdata_SE);
end

kdata = single(kdata);
k_nav = single(k_nav);

