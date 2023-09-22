% EPTI raw-data processing, image reconstruction and preprocessing pipeline, which supports:
% - 2D-EPTI sequence
% - gradient echo(GE) and/or asymmetric spin-echo(ASE)/spin-echo(SE) EPTI;
% - single-shot and multi-shot EPTI;
% - EPTI with or without SMS acquisition;
% - 3T and 7T EPTI data.

% This script produces reconstructed images from raw Twix data, and saves
% NIFTI output and some intermediate MATLAB files.
% See details about the output in README.txt.

% The parameters were pre-set for the example protocol released in the EPTI C2P sequence package.
% (You can request the EPTI sequence for Siemens scanners through MGH C2P Program http://nmr.mgh.harvard.edu/c2p or Siemens teamplay C2P Exchange).

% All MATLAB functions needed were included in the 'funcs_EPTI' folder.

% The script will call some BART functions, so BART installation is required. 
% Please download bart-0.8.00, then copy the EPTI customized files included
% in 'bart-0.8.00_EPTI/' to the BART folder, and overwrite the files including the 'Makefile' and source files in 'src/'.

% Before running this script, put in 
% - 'bart_path'; the path for BART;
% - 'directory': current working directory, data will be saved here;
% - 'filename_calib': the path for calibration rawdata
% - 'filename_EPTI': the path for imaging rawdata;
% - 'save_filename': output filename;
% Check these parameters to see if you are using the correct script that matches your data:
% - 'seq_type': single-shot or multi-shot;
% - 'process_type': gradient echo or spin-echo;
% - 'GPU_flag': GPU is recommended. if not available, set to 0;
% Adjust other parameters accordingly if using different EPTI protocols other than the provided example protocols.

% All sequence and reconstruction parameters were optimized for in-vivo
% human brain data. Phantom data may need different set of parameters.

% Please cite the following work on EPTI for this, including: 
% 1. "Echo planar time-resolved imaging (EPTI)", Fuyixue Wang et al., MRM 2019 Jun;81(6):3599-3615;
% 2. "Echo planar time-resolved imaging with subspace reconstruction and optimized spatiotemporal encoding", Zijing Dong et al., MRM 2020 Nov;84(5):2442-2455;
% 3. "Improving fMRI acquisition using single-shot EPTI with distortion-free high-SNR high-CNR multi-echo imaging", Fuyixue Wang et al., ISMRM 2022 p3330;
% Other relevant citations that may warrant consideration (EPTI's application in fMRI, dMRI and qMRI), see EPTI website: https://martinos.org/~fw089/

% Reconstruction version date: 09/2023; 
% We will continue to optimize the reconstruction and the sequence.
% Updated reconstruction code, and available sequences will be posted on our website https://martinos.org/~fw089/.

% Fuyixue Wang <fwang18@mgh.harvard.edu>; Zijing Dong <zdong3@mgh.harvard.edu>; July 2023, MGH
%%
clear;
addpath(genpath('../funcs_EPTI'));
%% bart path
bart_path = '';                        % path of BART toolbox
setenv('TOOLBOX_PATH', bart_path);
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));
%% write data path and filename
directory = '';           % please put the working directory, data will be saved here
filename_calib = fullfile(directory,'meas_MID00176_FID17047_EPTI_p2_calibration_2iso.dat');       % path for rawdata file of the calibration scan
filename_EPTI = fullfile(directory,'meas_MID00178_FID17049_EPTI_p2_singleshot_SE_2iso.dat');      % path for rawdata file of the imaging scan
save_filename = 'EPTI_p2_singleshot_SE_2iso';                                                     % filename for output data

[sub_dir,subdir1_acq,subdir2_recon,subdir3_nii] = EPTI_makeDir(directory,save_filename);
%% set parameters
% set basic parameters
seq_type = 0;           % '0': single-shot EPTI version;  '1': multi-shot EPTI version
process_type = 1;       % '0': process GE EPTI data;      '1': process SE EPTI data
Dyn_Dummy = 2;          % the number of early dynamcis to skip as dummy scans which have not reached steady state
GPU_flag = 1;           % '1': GPU; '0': CPU. GPU is recommended for faster reconstruction
show_image_flag = 0;    % show reconstructed images for each slice

% reconstruction parameters

recon_algorithm.regularization = 1;         % regularization for subspace reconstruction. '0': linear, no regularization; '1': locally low-rank
if recon_algorithm.regularization == 1
    recon_algorithm.rho = 1e-2;             % rho for ADMM
    recon_algorithm.lambda = 1e-3;          % lambda for regularization
    recon_algorithm.block_size = 8;         % block size for LLR
end
recon_algorithm.PF_recontype = 1;           % If there is partial fourier in the data, what type of partial fourier reconstruction to use. '0': zero-filling; '1': POCS
recon_algorithm.apodization = 1;            % perform apodization to reduce Gibbs artifact. '1': yes; '0': no
recon_algorithm.apodize_factor = 0.1;       % apodization factor
espirit_num_acs = 24;                       % ESPIRiT parameters for 3D coil sensitivity estimation, calbiration size, should be larger than # of slices
espirit_c = 0.1;                            % ESPIRiT parameters for 3D coil sensitivity estimation, threshold for eigenvalue

im_scaling = 100;                           % global scaling factor for image intensity        

% phase correction parameters
output_ghost_coeff = 1;                     % estimate odd-even echo phase difference from calibraiton

% if (process_type == 0) && (seq_type == 1)
%     Time_varing_phasecor = 1;               % apply time-varing odd-even phase difference correction to early echoes, can improve image quality for 7T GE acquisition
% else
%     Time_varing_phasecor = 0;               % apply time-constant odd-even phase difference (the time-varing phase from calibration does not apply to SE)
% end
Time_varing_phasecor = 0;                   % Turn off for 3T data
StartEcho_forPhasecor = 3;                  % start from which echo in calibration data to calculate odd-even echo phase difference for all echoes (skip early echoes which may have strong variation at 7T)
%% preparation and preprocessing
BART_savepath = [directory,'tmp/',save_filename,'_BART_tmp/'];
mkdir(BART_savepath);

% data name
if process_type == 0
    save_filename = [save_filename,'_GE'];
else
    save_filename = [save_filename,'_SE'];
end

% read rawdata info and calibration data
disp('---------------- preparation, reading EPTI raw data info ----------------');
meas = ReadRaw_EPTI_image_info_v2(filename_EPTI);    % read the image scan parameters

% reconstruction parameters
if (((meas.prot.lBaseResolution > 100) && (meas.prot.Nseg<=5))) || ((process_type == 0) && (meas.prot.nechoGE<=70))
   K = 3;   % number of subsapce bases, increase if more contrast/data are acquired, smaller K has higher SNR
else
   K = 4;  
end

switch  recon_algorithm.regularization
    case 0
        apx = ['_linear_K',num2str(K)];
    case 1
        apx = ['_lambda',num2str(recon_algorithm.lambda),'_rho',num2str(recon_algorithm.rho),'_LLR_K',num2str(K)];
end

% reconstruction basis
switch process_type  
    case 0
        [Phi,TEs_all,S] = EPTI_gen_Subspace_Basis_GE(meas.prot,K);      % generate subspace bases: GE
        scaling_for_Phi = (S(1:K).^(1/3))';
        Phi = Phi.*scaling_for_Phi;
    case 1
        [Phi,TEs_all,S] = EPTI_gen_Subspace_Basis_SE(meas.prot,K);      % generate subspace bases: SE
end
save(fullfile(subdir1_acq,['Basis_',save_filename,'.mat']),'Phi','S');

if (seq_type == 0) && (meas.prot.PF_shift>0)  % k-t partial fourier phase shift
    npe = meas.prot.Rseg*meas.prot.Nseg;
    mid=ceil((npe)/2);
    shift = ceil((npe+meas.prot.PF_shift)/2)-mid;
    linearphase_PF = permute(exp(-1i*2*pi* shift/npe*((0:1:npe-1)-mid)'),[2 1]);
else
    linearphase_PF = 1;
end

dt = meas.prot.iEffectiveEpiEchoSpacing*1e-6;
echo_cut1 = 0;
echo_cut2 = 0;
TEs_all_save = TEs_all(echo_cut1+1:end-echo_cut2);
meas.prot.echo_cut1 = echo_cut1;
meas.prot.echo_cut2 = echo_cut2;
if process_type == 1
    meas.prot.TEs_all_save = TEs_all_save + meas.prot.alTE(2)*1e-6;
else
    meas.prot.TEs_all_save = TEs_all_save;
end

shift_y = getEPTIimageshift_NCOoff(meas.prot)/( meas.prot.dPhaseFOV / (meas.prot.Nseg*meas.prot.Rseg+meas.prot.PF_shift) );

meas.prot.Dyn_Dummy = Dyn_Dummy;
save(fullfile(subdir1_acq,['meas_prot_',save_filename,'.mat']),'meas');
%%
disp('---------------- reading calibration data ----------------');
[kdata_calib,parameters_calib,linear_fit_coeff] = ReadRaw_EPTI_calib(filename_calib,output_ghost_coeff,meas.prot.NslicesEX);
kdata_calib = squeeze(kdata_calib);
save(fullfile(subdir1_acq,['meas_prot_Calib_',save_filename,'.mat']),'parameters_calib');

disp('---------------- Calculating sensitivity and B0 ----------------');
[sens_map_all,P_dB_all,Phase0_odd_all,Phase0_even_all,gccmtx_aligned,im_GRE_comb] = EPTI_Preprocessing_B0_sens_BARTv5(kdata_calib,meas.prot,parameters_calib,espirit_num_acs,espirit_c);
save(fullfile(subdir1_acq,['PreRecon_SensB0_',save_filename,'.mat']),'sens_map_all','P_dB_all','Phase0_odd_all','Phase0_even_all','gccmtx_aligned','im_GRE_comb','-v7.3');
%% image reconstruction (dynamic loop)
disp('---------------- EPTI Image Reconstruction start ----------------');
dyn_to_recon = Dyn_Dummy+1:meas.prot.N_dyn;
slice_to_recon = 1:meas.prot.N_slice;

nx = meas.prot.lBaseResolution;
ny = meas.prot.Nseg*meas.prot.Rseg + meas.prot.PF_shift; 
nsl = meas.prot.N_slice*meas.prot.NslicesEX; 
nDyn = length(dyn_to_recon);
necho = length(meas.prot.TEs_all_save);

[kdata,k_nav_first] = ReadRaw_EPTI_image(filename_EPTI,meas.prot,max(Dyn_Dummy+1,1),process_type,seq_type,gccmtx_aligned,linear_fit_coeff); % read the first dynamic of image scan
scaling_all = calculate_imagescaling_EPTI_useCalib(im_GRE_comb);
scaling_all = scaling_all/mean(scaling_all(:));
for dif_dyn=dyn_to_recon
    disp(['-------- reconstructing dynamic #',num2str(dif_dyn),'-------- ']);

    [kdata,k_nav] = ReadRaw_EPTI_image_PhaseCor_general(filename_EPTI,meas.prot,dif_dyn,process_type,seq_type,gccmtx_aligned,linear_fit_coeff,dt,k_nav_first,TEs_all);
    slice_sep = meas.prot.sSliceArray_lSize/meas.prot.sWipMemBlock_adFree(1);
    scaling_dyn = calculate_imagescaling_EPTIv2(kdata).*scaling_all;
    fprintf('slice: ');
    im_epti = zeros(nx,ny,nsl,necho);
    for dif_slice=slice_to_recon            % slice loop
        fprintf('%d ',dif_slice);
        if meas.prot.NslicesEX > 1
            slice_group = dif_slice:slice_sep:slice_sep*meas.prot.NslicesEX;
            slice_group = slice_group(end:-1:1);
            kdata_tmp = kdata(:,:,:,end:-1:1,:,dif_slice);
        else
            slice_group = dif_slice;
            kdata_tmp = permute(kdata(:,:,:,:,dif_slice),[1 2 3 5 4]);
        end
        scaling = max(scaling_dyn(slice_group));
        sens_map = sens_map_all(:,:,slice_group,:);
        P_dB = P_dB_all(:,:,slice_group);
        
        if seq_type == 1
            Phase0_odd = Phase0_odd_all(:,:,:,slice_group);
            Phase0_even = Phase0_even_all(:,:,:,slice_group);
            Phase0_odd=exp(1i*Phase0_odd).*linearphase_PF;
            Phase0_even=exp(1i*Phase0_even).*linearphase_PF;
            Phase0_use_t = zeros(size(Phase0_odd,1),size(Phase0_odd,2),size(Phase0_odd,4),size(kdata_tmp,1));
            Phase0_use_t(:,:,:,1:2:end) = permute(repmat(exp(1i*angle(mean(Phase0_odd(:,:,min(StartEcho_forPhasecor,size(Phase0_odd,3)):end,:),3))),[1 1 length(1:2:size(Phase0_use_t,4)) 1]),[1 2 4 3]);
            Phase0_use_t(:,:,:,2:2:end) = permute(repmat(exp(1i*angle(mean(Phase0_even(:,:,min(StartEcho_forPhasecor,size(Phase0_even,3)):end,:),3))),[1 1 length(2:2:size(Phase0_use_t,4)) 1]),[1 2 4 3]);
            if Time_varing_phasecor == 1
                Phase0_use_t(:,:,:,1:2:size(Phase0_odd,3)*2) = permute(Phase0_odd,[1 2 4 3]);
                Phase0_use_t(:,:,:,2:2:size(Phase0_even,3)*2) = permute(Phase0_even,[1 2 4 3]);
            end
        else
            Phase0_use_t = ones(size(kdata_tmp,2),size(kdata_tmp,3),size(kdata_tmp,4),size(kdata_tmp,1)).*linearphase_PF;  % skip high order phase correction
        end

        [im_recon,a] = SubspaceRecon_Process_v4(kdata_tmp,sens_map,P_dB,...
            Phase0_use_t,TEs_all,Phi,meas.prot,scaling,BART_savepath,recon_algorithm,GPU_flag); 
        im_recon = im_recon(:,:,:,echo_cut1+1:end-echo_cut2);
        if show_image_flag == 1
            target_echo = 5:size(im_recon,4)-5;
            im_recon_show = abs(permute(im_recon,[2 1 3 4]));  
            im_recon_show = sos(im_recon_show(:,:,:,target_echo),4);
            figure; imshow3(im_recon_show(end:-1:1,:,:),[],[1,meas.prot.NslicesEX]);
        end
        im_epti(:,:,dif_slice+(0:meas.prot.NslicesEX-1)*meas.prot.N_slice,:)=im_recon(:,:,end:-1:1,:);
    end
    im_epti = single(im_epti);
    save(fullfile(subdir2_recon,['Recon_EPTI_',save_filename,apx,'_Dyn_',num2str(dif_dyn),'.mat']),'im_epti');
    fprintf('\n');
end
disp('---------------- image reconstruction finished ----------------')
%% EPTI data preprocessing
clear im_epti;
disp('---------------- image preprocessing start ----------------')
im_epti_final = single(zeros(nx,ny,nsl,necho,nDyn));
for dif_dyn=1:nDyn
    disp(['dyn:',num2str(dyn_to_recon(dif_dyn))]);
    load(fullfile(subdir2_recon,['Recon_EPTI_',save_filename,apx,'_Dyn_',num2str(dyn_to_recon(dif_dyn)),'.mat']));
    im_epti_final(:,:,:,:,dif_dyn) = im_epti;
end
im_epti_final = im_epti_final*(im_scaling/mean(abs(im_epti_final(:))));
im_epti_final = circshift(im_epti_final,round(shift_y),2);

im_epti_sos = squeeze(sos(im_epti_final,4));
saveniftidata_EPTI(im_epti_sos,meas,fullfile(subdir3_nii,[save_filename,'_im_sos.nii']));

im_epti_Tave = mean(abs(im_epti_final),5);
saveniftidata_EPTI(im_epti_Tave,meas,fullfile(subdir3_nii,[save_filename,'_im_Tave_echoes.nii']));
writematrix(meas.prot.TEs_all_save,fullfile(subdir3_nii,[save_filename,'__TEs.txt']));

for dif = 1:necho
    disp(['---------------------saving echo',num2str(dif)]);
    saveniftidata_EPTI(squeeze(abs(im_epti_final(:,:,:,dif,:))),meas,fullfile(subdir3_nii,'echoes',[save_filename,'_im_echo',num2str(dif),'.nii']));
end

if(strcmp(save_filename(end-1:end),'GE'))
    [img_comb,T2s_map] = optimalCombineEchoes(im_epti_final,meas.prot.TEs_all_save);
    saveniftidata_EPTI(T2s_map,meas,fullfile(subdir3_nii,[save_filename,'__T2s.nii']));
    saveniftidata_EPTI(img_comb,meas,fullfile(subdir3_nii,[save_filename,'_im_comb.nii'])); 
end
%%
disp('done!');

disp(['Note: the first ',num2str(Dyn_Dummy),' dynamics (',num2str(Dyn_Dummy*meas.prot.alTR*meas.prot.Nseg*1e-6),' seconds) were removed as dummy scans, please adjust your fMRI analysis accordingly']);