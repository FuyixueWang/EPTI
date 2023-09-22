function [meas, linear_fit_coeff, linear_fit_coeff_acs] = GhostCorrectAndGrid_mkm_PSFv2(meas,linear_fit_coeff,processNav)
% ghost correct and grid the data and navigator 
% note: output data is collapsed in segment direction but navigator is not

% Zijing Dong, 2018, make it work for PSF data (based on MaryKate's previous version)

% Fuyixue Wang, 2019, make it work for EPTI data

% updated to work for single-shot and multi-shot EPTI data
% Fuyixue Wang, 2023

if nargin <3
    processNav = 0;
end
if nargin <2
    linear_fit_coeff = [];
end

if isempty(linear_fit_coeff)
    precompute_coeff = 0;
    linear_fit_coeff = mrir_artifact_ghost_compute(meas.data_phascor1d);
    %linear_fit_coeff = mrir_artifact_ghost_compute_v2(meas.data_phascor1d);
else
    precompute_coeff = 1;
end

mask_FWD = squeeze(meas.data(end/2,:,1,1,end,1,1,1,1,1) ~= 0);
mask_REV = squeeze(meas.data(end/2,:,1,1,end,1,1,2,1,1) ~= 0);

data_hybrid = mrir_iDFT_freqencode(meas.data); 
meas.data = [];
data_hybrid = mrir_artifact_ghost_correct(data_hybrid, linear_fit_coeff);
data_hybrid = sum(data_hybrid,8); 
data_hybrid = mrir_regrid_trapezoid(data_hybrid, meas.prot);

meas.data = mrir_fDFT_freqencode(data_hybrid); 
clear data_hybrid

s = size(meas.data); s(8) = 2; s(s==0)=1;
measTemp = zeros(s);
measTemp(:,mask_FWD,:,:,:,:,:,1,:,:) = meas.data(:,mask_FWD,:,:,:,:,:,:,:,:);
measTemp(:,mask_REV,:,:,:,:,:,2,:,:) = meas.data(:,mask_REV,:,:,:,:,:,:,:,:);
meas.data = measTemp; clear measTemp;

if processNav == 1
    nav_hybrid = mrir_iDFT_freqencode(meas.data_phascor1d);
    meas.data_phascor1d = [];
    nav_hybrid = mrir_artifact_ghost_correct( nav_hybrid, linear_fit_coeff);
    nav_hybrid = mrir_regrid_trapezoid(nav_hybrid, meas.prot);
    meas.data_phascor1d = mrir_fDFT_freqencode(nav_hybrid);
    clear nav_hybrid
end

if (~precompute_coeff)
    %linear_fit_coeff3 =
    %mean(reshape(linear_fit_coeff,2,[],size(meas.data_phascor1d,7)),3); %this is wrong!!
    dims = size(meas.data_phascor1d);
    if length(dims) < 10
        dims(end+1:10) = 1; 
    end
    linear_fit_coeff = reshape( median( reshape(linear_fit_coeff, [2 1 dims(3:7) 1 1 dims(10)]), 7) , 2,[]) ; % for saving to use for tailor ghost correction   
end

%%  add section to process offline data (DPG for SMS ref scan) 

if isfield(meas, 'offline') 
    if (~isempty(meas.offline))
    
        mask_FWD = squeeze(meas.offline(end/2,:,1,1,end,1,1,1,1,1) ~= 0);
        mask_REV = squeeze(meas.offline(end/2,:,1,1,end,1,1,2,1,1) ~= 0);
       % meas.offline = mrir_artifact_ghost_correct_CorrelationMethod_v4(meas.offline, linear_fit_coeff);
        data_hybrid = mrir_iDFT_freqencode(meas.offline); 
        meas.offline = [];
        data_hybrid = mrir_artifact_ghost_correct(data_hybrid, linear_fit_coeff);
        data_hybrid = sum(data_hybrid,8); 
        data_hybrid = mrir_regrid_trapezoid(data_hybrid, meas.prot);
        meas.offline = mrir_fDFT_freqencode(data_hybrid); 
        clear data_hybrid

        s = size(meas.offline); s(8) = 2; measTemp = zeros(s);
        measTemp(:,mask_FWD,:,:,:,:,:,1,:,:) = meas.offline(:,mask_FWD,:,:,:,:,:,:,:,:);
        measTemp(:,mask_REV,:,:,:,:,:,2,:,:) = meas.offline(:,mask_REV,:,:,:,:,:,:,:,:);
        meas.offline = measTemp; clear measTemp;
    end
end
%% for accelerated data
if isfield(meas, 'patrefscan')
    if (~isempty(meas.patrefscan))
        if (~precompute_coeff)
            linear_fit_coeff_acs = mrir_artifact_ghost_compute(meas.patrefscan_phascor);
        else
            linear_fit_coeff_acs = linear_fit_coeff;
        end
          
        mask_FWD = squeeze(meas.patrefscan(end/2,:,1,1,end,1,1,1,1,1) ~= 0);
        mask_REV = squeeze(meas.patrefscan(end/2,:,1,1,end,1,1,2,1,1) ~= 0);
        acs_hybrid = mrir_iDFT_freqencode(meas.patrefscan);
        meas.patrefscan = [];
        acs_hybrid = mrir_artifact_ghost_correct(acs_hybrid, linear_fit_coeff_acs);
        acs_hybrid = sum(acs_hybrid,8);
        acs_hybrid = mrir_regrid_trapezoid(acs_hybrid, meas.prot);
        meas.patrefscan = mrir_fDFT_freqencode(acs_hybrid);
        clear acs_hybrid

        s = size(meas.patrefscan); s(8) = 2; 
        s(s==0)=1; measTemp = zeros(s);
        measTemp(:,mask_FWD,:,:,:,:,:,1,:,:) = meas.patrefscan(:,mask_FWD,:,:,:,:,:,:,:,:);
        measTemp(:,mask_REV,:,:,:,:,:,2,:,:) = meas.patrefscan(:,mask_REV,:,:,:,:,:,:,:,:);
        meas.patrefscan = measTemp; clear measTemp;

        if isfield(meas, 'patrefscan_dpg') 
            mask_FWD = squeeze(meas.patrefscan_dpg(end/2,:,1,1,end,1,1,1,1,1) ~= 0);
            mask_REV = squeeze(meas.patrefscan_dpg(end/2,:,1,1,end,1,1,2,1,1) ~= 0);
            %meas.patrefscan_dpg = mrir_artifact_ghost_correct_CorrelationMethod_v4(meas.patrefscan_dpg, linear_fit_coeff_acs);
            acs_hybrid = mrir_iDFT_freqencode(meas.patrefscan_dpg);
            meas.patrefscan_dpg = [];
            acs_hybrid = mrir_artifact_ghost_correct(acs_hybrid, linear_fit_coeff_acs);
            acs_hybrid = sum(acs_hybrid,8);
            acs_hybrid = mrir_regrid_trapezoid(acs_hybrid, meas.prot);
            meas.patrefscan_dpg = mrir_fDFT_freqencode(acs_hybrid);
            clear acs_hybrid

            s = size(meas.patrefscan_dpg); s(8) = 2; measTemp = zeros(s);
            measTemp(:,mask_FWD,:,:,:,:,:,1,:,:) = meas.patrefscan_dpg(:,mask_FWD,:,:,:,:,:,:,:,:);
            measTemp(:,mask_REV,:,:,:,:,:,2,:,:) = meas.patrefscan_dpg(:,mask_REV,:,:,:,:,:,:,:,:);
            meas.patrefscan_dpg = measTemp; clear measTemp;
        end
        
        if processNav == 1
            acs_nav_hybrid = mrir_iDFT_freqencode(meas.patrefscan_phascor);
            meas.patrefscan_phascor = [];
            acs_nav_hybrid = mrir_artifact_ghost_correct(acs_nav_hybrid, linear_fit_coeff_acs);
            acs_nav_hybrid = mrir_regrid_trapezoid(acs_nav_hybrid, meas.prot);
            meas.patrefscan_phascor = mrir_fDFT_freqencode(acs_nav_hybrid);
            clear acs_nav_hybrid
        end
    end
end














