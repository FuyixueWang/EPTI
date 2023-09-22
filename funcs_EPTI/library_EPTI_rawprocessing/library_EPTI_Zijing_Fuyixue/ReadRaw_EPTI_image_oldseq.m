function [kdata,k_nav,parameters,linear_fit_coeff] = ReadRaw_EPTI_image_oldseq(filename,parameters,dyn,process_type,seq_type,gccmtx_aligned,linear_fit_coeff)
% EPTI rawdata processing for imaging data
% including reading rawdata, raw data processing, and put into EPTI 
% encoding pattern
% Fuyixue Wang, 2021, MGH

% Generalzed code for both single-shot and multi-shot EPTI data
% support applying odd-even echo eddy-current phase correction
% support to both GE and SE EPTI data with SMS
% Fuyixue Wang, Zijing Dong, MGH, 2023
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
        if parameters.pf_echo > 0
            kdata = kdata_all(1:floor(parameters.pf_echo*size(kdata_all,1)),:,:,:,:,1);
        else
            kdata = kdata_all(:,:,:,:,:,1);
        end
        if seq_type == 0 % zigzag/single-shot version
            kdata = putRawtoPattern_SMS_zigzag_PF_selfnav_oldseq(kdata,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift,parameters.PF_shift);
        else             % multi-shot version
            kdata = putRawtoPattern_EPTI_SMS(kdata,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift);
            kdata = kdata(:,:,end:-1:1,:,:,:);
        end
        
    case 1  % 0: SE only
        kdata=kdata_all(1:end,:,:,:,:,2);
        if seq_type == 0 % zigzag/single-shot version
            kdata = putRawtoPattern_SMS_zigzag_PF_selfnav_oldseq(kdata,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift,parameters.PF_shift);
        else             % multi-shot version
            kdata = putRawtoPattern_EPTI_SMS(kdata,parameters.Nseg,parameters.Rseg,parameters.Rpe,parameters.SMS_shift);
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

