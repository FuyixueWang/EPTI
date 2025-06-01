function [sens_map_all,P_dB_all,Phase0_odd_all,Phase0_even_all,gccmtx_aligned,im_GRE_combo] = EPTI_Preprocessing_B0_sens_BARTv5(kdata_calib0,parameters,parameters_calib,num_acs,c,coilcompression_flag) 
% Calculate the sensitivity map and B0 map for EPTI reconstruction
% estimating using EPTI k-t calibraiton data
% Support accelerated EPTI calibration scan
% This version uses BART for sensitivity estimation

% Please cite the following work on EPTI for this, including: 
% 1. Dong Z, Wald LL, Polimeni JR, Wang F. Single-shot echo planar time-resolved imaging for multi-echo functional MRI and distortion-free diffusion imaging. Magn Reson Med. 2025; 93(3): 993-1013. doi: 10.1002/mrm.30327
% 2. Wang F, Dong Z, Reese TG, et al. Echo planar time-resolved imaging (EPTI). Magn Reson Med. 2019 Jun;81(6):3599-3615. doi: 10.1002/mrm.27673.
% 3. Dong Z, Wang F, Reese TG, et al. Echo planar time-resolved imaging with subspace reconstruction and optimized spatiotemporal encoding. Magn Reson Med. 2020 Nov;84(5):2442-2455. doi: 10.1002/mrm.28295.
% Other relevant citations that may warrant consideration (EPTI's application in fMRI, dMRI and qMRI), see EPTI website: https://martinos.org/~fw089/

% v2: output the odd-even echo phase for all the echoes of the calibration
% data, this can be useful for estimation time-varing eddy-current changes

% v3: chagne from 2D Espirit to 3D Espirit to improve the consistency of 
% sensitivity map across slices

% v4-5: Optimized GRAPPA reconstruction, optimized zero padding to
% interpolation for B0 calculation to avoid ringings

% Fuyixue Wang, Zijing Dong, 2023, MGH


if nargin<6
    coilcompression_flag = 1;
end
N_slice = parameters_calib.sSliceArray_lSize;
nc = parameters.iMaxNoOfRxChannels;
nx = parameters.lBaseResolution;
npe = parameters.Rseg.*parameters.Nseg;
dt_calib = parameters_calib.iEffectiveEpiEchoSpacing*1e-6; % echo spacing
t0 = parameters_calib.alTE(1)*1e-6 - (parameters_calib.nechoGE/2)*dt_calib;
TEs_calib=(dt_calib:dt_calib:parameters_calib.nechoGE*dt_calib)+t0;

ACS_npe = 18; 
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

%% GRAPPA conv
KernalSize(1) = 3;
KernalSize(2) = 3;
lam = 0.01;
kdata_calib_cc = permute(kdata_calib_cc,[2 3 5 1 4]);

if PAT_factor>1
    k_calib_rec = zeros(size(kdata_calib_cc));
    CalibSize = [nx,ACS_npe];
    
    for slice = 1:size(kdata_calib_cc,5)
        acs = crop(squeeze(kdata_calib_cc(:,:,:,1,slice)),[CalibSize,ncc]);
        G = mrir_array_GRAPPA_2d_kernel_improved_dev(acs, PAT_factor, 1, KernalSize(1), KernalSize(2) , 1, 1, 1,0,lam,nx);    
        k = mrir_array_GRAPPA_conv_kernel_v2(G);
        for echo = 1:size(kdata_calib_cc,4)
            k_fill = zeros(size(kdata_calib_cc,1),size(kdata_calib_cc,2),size(kdata_calib_cc,3));
            k_fill(:,1:PAT_factor:end,:) = kdata_calib_cc(:,1:PAT_factor:end,:,echo,slice);
            k_conv_cha = complex(zeros(nx, Nseg, ncc, ncc));
            for trg = 1:ncc
              for src = 1:ncc
                  k_conv_cha(:,:,src,trg) = conv2( k_fill(:,:,src), k(:,:,1,src,trg), 'same');
              end
            end
            k_fill = squeeze(sum(k_conv_cha, 3));
            k_fill(:,index,:) = kdata_calib_cc(:,index,:,echo,slice);
            k_calib_rec(:,:,:,echo,slice) = k_fill;
        end
    end
else
    k_calib_rec = kdata_calib_cc;
end
%% Sens map and B0
kdata_calib_cc = permute(k_calib_rec,[1 2 5 3 4]);
im_GRE = ifft2c(zpad(kdata_calib_cc,[nx,npe,N_slice,ncc,size(kdata_calib_cc,5)]));
% im_GRE = zeros([nx,npe,N_slice,ncc,size(kdata_calib_cc,5)]);
% im_calib_cc = ifft2c(kdata_calib_cc);
% for dif_sz3=1:N_slice
%     for dif_sz4=1:ncc
%         for dif_sz5=1:size(kdata_calib_cc,5)
%         if mod(size(kdata_calib_cc,2),2) == 1
%             im_GRE(:,:,dif_sz3,dif_sz4,dif_sz5) = imtranslate(imresize( im_calib_cc(:,:,dif_sz3,dif_sz4,dif_sz5), [nx,npe]),[0.5,0]);
%         else
%             im_GRE(:,:,dif_sz3,dif_sz4,dif_sz5) = imtranslate(imresize( im_calib_cc(:,:,dif_sz3,dif_sz4,dif_sz5), [nx,npe]),[-0.5,0]);
%         end
%         end
%     end
% end

% num_acs = 24;
% c = 0.1;

pad_edgeslice = 5;
if (size(im_GRE,3)+pad_edgeslice*2)<num_acs
    pad_edgeslice = pad_edgeslice+ceil((num_acs-size(im_GRE,3)+pad_edgeslice*2)/2);
end
pad_edgeslice_end = pad_edgeslice;
if mod(size(im_GRE,3),2)==1   % odd matrix size may cause issue for BART
    pad_edgeslice_end = pad_edgeslice_end+1;
end

im_echo1 = zeros(size(im_GRE,1),size(im_GRE,2),size(im_GRE,3)+pad_edgeslice+pad_edgeslice_end,size(im_GRE,4));
im_echo1(:,:,pad_edgeslice+1:end-pad_edgeslice_end,:) = im_GRE(:,:,:,:,1); 
% im_echo1(:,:,1:pad_edgeslice,:) = repmat(im_echo1(:,:,pad_edgeslice+1,:),[1 1 pad_edgeslice,1]);
% im_echo1(:,:,end-pad_edgeslice_end+1:end,:) = repmat(im_echo1(:,:,end-pad_edgeslice_end,:),[1 1 pad_edgeslice_end,1]);

[calib,~] = bart('ecalib -d 0 -r ', num2str(num_acs), ' -c ', num2str(c), ' -S -m1 ', fft3c(im_echo1));
sens_map=calib(:,:,pad_edgeslice+1:end-pad_edgeslice_end,:,1);

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

end
