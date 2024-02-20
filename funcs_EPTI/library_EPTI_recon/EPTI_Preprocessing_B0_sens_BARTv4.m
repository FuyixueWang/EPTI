function [sens_map_all,P_dB_all,Phase0_odd_all,Phase0_even_all,gccmtx_aligned,im_GRE_combo] = EPTI_Preprocessing_B0_sens_BARTv4(directory,filename,kdata_calib0,parameters,parameters_calib,save_data_flag,coilcompression_flag) 
% Calculate the sensitivity map and B0 map for EPTI reconstruction
% estimating using EPTI k-t calibraiton data
% Support accelerated EPTI calibration scan
% This version uses BART for sensitivity estimation

% Please cite the following work on EPTI for this, including: 
% 1. "Echo planar time-resolved imaging (EPTI)", Fuyixue Wang et al., MRM 2019 Jun;81(6):3599-3615;
% 2. "Echo planar time-resolved imaging with subspace reconstruction and optimized spatiotemporal encoding", Zijing Dong et al., MRM 2020 Nov;84(5):2442-2455;
% 3. "Single-shot Echo Planar Time-resolved Imaging for multi-echo functional MRI and distortion-free diffusion imaging", Dong Z, Wald LL, Polimeni JR, Wang F. bioRxiv. 2024 Jan 26:2024.01.24.577002.
% Other relevant citations that may warrant consideration (EPTI's application in fMRI, dMRI and qMRI), see EPTI website: https://martinos.org/~fw089/

% Fuyixue Wang, Zijing Dong, 2022, MGH

% v2: output the odd-even echo phase for all the echoes of the calibration
% data, this can be useful for estimation time-varing eddy-current changes
% Fuyixue Wang, Zijing Dong, 2023, MGH

if nargin<7
    coilcompression_flag = 1;
end
N_slice = parameters_calib.sSliceArray_lSize;
nc = parameters.iMaxNoOfRxChannels;
nx = parameters.lBaseResolution;
npe = parameters.Rseg.*parameters.Nseg;
dt_calib = parameters_calib.iEffectiveEpiEchoSpacing*1e-6; % echo spacing
t0 = parameters_calib.alTE(1)*1e-6 - (parameters_calib.nechoGE/2)*dt_calib;
TEs_calib=(dt_calib:dt_calib:parameters_calib.nechoGE*dt_calib)+t0;

%% v4
ACS_npe = parameters_calib.ACS_npe; %%% Fuyixue testing
PAT_factor = parameters_calib.sWipMemBlock_alFree(35);
Nseg = parameters_calib.sWipMemBlock_alFree(34);

n_pe = ACS_npe + ceil((Nseg-ACS_npe)/PAT_factor);
index=[];
for i = 1:n_pe
    if i<(ceil((Nseg-ACS_npe)/PAT_factor)/2)+1
        index(i) = -(Nseg-1)/2 + PAT_factor*(i-1);
    elseif i >= ((ceil((Nseg-ACS_npe)/PAT_factor)/2) + ACS_npe + 1)
        index(i) = -(Nseg-1)/2 + PAT_factor*(i-ACS_npe+round(ACS_npe/PAT_factor)-1);
    else
        index(i) = -(ACS_npe)/2 + (i-floor(ceil((Nseg-ACS_npe)/PAT_factor)/2)-1);
    end
%     index(i) = index(i)+1;
end
index = index + (Nseg-1)/2 +1;

% % v2
% ACS_npe = 12;
% PAT_factor = parameters_calib.sWipMemBlock_alFree(35);
% Nseg = parameters_calib.sWipMemBlock_alFree(34);
% 
% n_pe = ACS_npe + floor((Nseg-ACS_npe)/PAT_factor);
% for i = 1:n_pe
%     if i<(floor((Nseg-ACS_npe)/PAT_factor)/2)+1
%         index(i) = -(Nseg-1)/2 + PAT_factor*(i-1);
%     elseif i >= ((floor((Nseg-ACS_npe)/PAT_factor)/2) + ACS_npe + 1)
%         index(i) = -(Nseg-1)/2 + PAT_factor*(i-ACS_npe+round(ACS_npe/PAT_factor)-1);
%     else
%         index(i) = -(ACS_npe)/2 + (i-floor(floor((Nseg-ACS_npe)/PAT_factor)/2)-1-1);
%     end
%     index(i) = index(i)+1;
% end
% index = index + round((Nseg+1)/2);

% v3
% ACS_npe = 12;
% PAT_factor = parameters_calib.sWipMemBlock_alFree(35);
% Nseg = parameters_calib.sWipMemBlock_alFree(34);
% 
% n_pe = ACS_npe + floor((Nseg-ACS_npe)/PAT_factor);
% index=[];
% for i = 1:n_pe
%     if i<(floor((Nseg-ACS_npe)/PAT_factor)/2)+1
%         index(i) = -(Nseg-1)/2 + PAT_factor*(i-1);
%     elseif i >= ((floor((Nseg-ACS_npe)/PAT_factor)/2) + ACS_npe + 1)
%         index(i) = -(Nseg-1)/2 + PAT_factor*(i-ACS_npe+round(ACS_npe/PAT_factor)-1);
%     else
%         index(i) = -(ACS_npe)/2 + (i-floor(floor((Nseg-ACS_npe)/PAT_factor)/2)-1);
%     end
% %     index(i) = index(i)+1;
% end
% index = index + round((Nseg+1)/2);

kdata_calib = zeros(size(kdata_calib0,1),size(kdata_calib0,2),Nseg,size(kdata_calib0,4),size(kdata_calib0,5));
kdata_calib(:,:,index,:,:) = kdata_calib0(:,:,1:n_pe,:,:);
clear kdata_calib0;
%%
kdata_calib = permute(kdata_calib,[1 2 3 5 4]);

if coilcompression_flag == 1
    dim = 1;
    ncalib_p = ACS_npe; %
    ncalib_z = 24; %
    if nc <= 40
        ncc = round(nc/2);
    else
        ncc = round(nc/3);
    end
    slwin = 7; % sliding window length for GCC
    calib = fftc(squeeze(kdata_calib(1,:,:,:,:)),3);
    ncalib_z = min(ncalib_z,size(calib,3));
    ncalib_p = min(ncalib_p,size(calib,2));
    
    calib = crop(calib,[nx,ncalib_p,ncalib_z,nc]);
    gccmtx = calcGCCMtx(calib,dim,slwin);
    gccmtx_aligned = alignCCMtx(gccmtx(:,1:ncc,:));
    
    % disp('Compressing calibration data ...');
    kdata_calib_cc = zeros(size(kdata_calib,1), size(kdata_calib,2), size(kdata_calib,3), size(kdata_calib,4), ncc);
    for dif_Necho=1:size(kdata_calib,1)
        for slice = 1:N_slice
            CCDATA_aligned = CC(squeeze(kdata_calib(dif_Necho,:,:,slice,:)),gccmtx_aligned, dim);
            kdata_calib_cc(dif_Necho,:,:,slice,:)=CCDATA_aligned;
        end
    end
else
    kdata_calib_cc = kdata_calib;
    gccmtx_aligned = [];
    ncc = nc;
end
%%
lambda_grappa = 0.1;
kdata_calib_cc = permute(kdata_calib_cc,[3 2 5 1 4]);
CalibSize = [ACS_npe,48];
k_calib_rec = zeros(size(kdata_calib_cc));
pe_fill = min(index(:)):max(index(:));

for slice = 1:size(kdata_calib_cc,5)
%     disp([num2str(slice),' / ',num2str(size(kdata_calib_cc,5))]);
    kcalib = crop(squeeze(kdata_calib_cc(:,:,:,1,slice)),[CalibSize,ncc]);
    [ kernel ] = calibration_GRAPPA(kcalib,PAT_factor,[2 5],lambda_grappa);    
    for echo  = 1:size(kdata_calib_cc,4)
        k_calib_rec(pe_fill,:,:,echo,slice) = GRAPPAREC( kdata_calib_cc(pe_fill,:,:,echo,slice),kernel,PAT_factor,[2 5] );
    end
end
%% Sens map and B0
kdata_calib_cc = permute(k_calib_rec,[2 1 5 3 4]);
% im_GRE = ifft2c(zpad(kdata_calib_cc,[nx,npe,N_slice,ncc,size(kdata_calib_cc,5)]));

im_calib_cc = ifft2c(kdata_calib_cc);
im_GRE = zeros([nx,npe,N_slice,ncc,size(kdata_calib_cc,5)]);
for dif_sz3=1:N_slice
    for dif_sz4=1:ncc
        for dif_sz5=1:size(kdata_calib_cc,5)
            im_GRE(:,:,dif_sz3,dif_sz4,dif_sz5) = imtranslate(imresize( im_calib_cc(:,:,dif_sz3,dif_sz4,dif_sz5), [nx,npe]),[0.5,0]);
        end
    end
end

num_acs = 20;
c = 0.01;
pad_edgeslice = 5;
im_echo1 = zeros(size(im_GRE,1),size(im_GRE,2),size(im_GRE,3)+pad_edgeslice*2,size(im_GRE,4));
im_echo1(:,:,pad_edgeslice+1:end-pad_edgeslice,:) = im_GRE(:,:,:,:,1); 
% im_echo1(:,:,1:pad_edgeslice,:) = repmat(im_echo1(:,:,pad_edgeslice+1,:),[1 1 pad_edgeslice,1]);
% im_echo1(:,:,end-pad_edgeslice+1:end,:) = repmat(im_echo1(:,:,end-pad_edgeslice,:),[1 1 pad_edgeslice,1]);

[calib,~] = bart('ecalib -d 0 -r ', num2str(num_acs), ' -c ', num2str(c), ' ', fft3c(im_echo1));
sens_map=calib(:,:,pad_edgeslice+1:end-pad_edgeslice,:,1);

% sens_map=zeros([nx,npe,N_slice,ncc]);
% volume_svd = svd_compress3d( im_echo1, 1 );
% for slice=1:N_slice
%     fprintf('%d/%d ',slice,N_slice);
%     [calib,~] = bart('ecalib -d 0 -r ', num2str(num_acs), ' -c ', num2str(c), ' ', fft2c(cat(4, 1e-6 * volume_svd(:,:,slice), im_echo1(:,:,slice,:))));
%     sens = bart('slice 4 0', calib);
%     sens=sens(:,:,:,2:end);
%     sens_map(:,:,slice,:) = sens;   
% end
% figure,imshow3(squeeze(abs(sens_map(:,:,2,:))),[],[4 4]);
im_GRE_combo = squeeze( sum(conj(sens_map).* im_GRE, 4) ./ (eps +  sum(abs(sens_map).^2, 4)) );
%%
PHS = angle(im_GRE_combo);
im_mean=abs(im_GRE_combo(:,:,:,1));
MSK_extended=(im_mean)>0.3*mean(im_mean(:));

warning off;
%%
echo_select = 1:size(PHS,4);
PHS = PHS(:,:,:,echo_select);
TEs_calib_use = TEs_calib(echo_select);
P_dB_all=zeros(nx,npe,N_slice);
Necho = size(im_GRE_combo,4);

Phase0_odd_all=zeros(nx,npe,length(1:2:Necho),N_slice);
Phase0_even_all=zeros(nx,npe,length(2:2:Necho),N_slice);

for slice=1:N_slice
    P_dB_all(:,:,slice) = dB_fitting_JumpCorrect(squeeze(PHS(:,:,slice,:)),TEs_calib_use(:),logical(MSK_extended(:,:,slice)),1);
    num=0;
    tmp=[];
    for t = 1:2:Necho
        num=num+1;
        tmp(:,:,num)=squeeze(im_GRE_combo(:,:,slice,t)).*exp(-1i*2*pi*P_dB_all(:,:,slice)*TEs_calib(t));
    end
    Phase0_odd_all(:,:,:,slice)=angle(tmp).*MSK_extended(:,:,slice);
    num=0;     
    tmp=[];
    for t = 2:2:Necho
        num=num+1;
        tmp(:,:,num)=squeeze(im_GRE_combo(:,:,slice,t)).*exp(-1i*2*pi*P_dB_all(:,:,slice)*TEs_calib(t));
    end
    Phase0_even_all(:,:,:,slice)=angle(tmp).*MSK_extended(:,:,slice);
end
%% save
sens_map_all = single(sens_map);
if save_data_flag
    save_path = [directory,'data/1_Recon_subspace/'];
    save([save_path,'PreRecon_SensB0Dephase_SMS_',filename,'.mat'],'sens_map_all','P_dB_all','Phase0_odd_all','Phase0_even_all',...
        'gccmtx_aligned','-v7.3');
end
end
