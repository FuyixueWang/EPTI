function [im_recon,a,P_dB] = SubspaceRecon_Process_B0update(kdata,sens_map,P_dB,Phase0_all,TEs,Phi,Phi_B0,parameters,scaling,bart_savepath,recon_algorithm,GPU_flag)
% EPTI subspace reconstruction using BART
% Zijing Dong, MGH, 2019

% Generalized for EPTI sequences
% Adding multiple options
% Fuyixue Wang, MGH, 2022

% self B0 updating for field drift correction
% Zjing Dong, MGH, 2022

% Support time varing background phase input
% Fuyixue Wang, MGH, 2023

% Please cite the following work on EPTI for this, including: 
% 1. "Echo planar time-resolved imaging (EPTI)", Fuyixue Wang et al., MRM 2019 Jun;81(6):3599-3615;
% 2. "Echo planar time-resolved imaging with subspace reconstruction and optimized spatiotemporal encoding", Zijing Dong et al., MRM 2020 Nov;84(5):2442-2455;
% 3. "Single-shot Echo Planar Time-resolved Imaging for multi-echo functional MRI and distortion-free diffusion imaging", Dong Z, Wald LL, Polimeni JR, Wang F. bioRxiv. 2024 Jan 26:2024.01.24.577002.
% Other relevant citations that may warrant consideration (EPTI's application in fMRI, dMRI and qMRI), see EPTI website: https://martinos.org/~fw089/

% Fuyixue Wang, Zijing Dong, 2023, MGH
%% set parameters
if GPU_flag==1
    gpu = '-g';
else
    gpu = '';
end
PF_shift = parameters.PF_shift;
ns = parameters.NslicesEX;
K=size(Phi,2);
K_B0=size(Phi_B0,2);

%% Generate Trajectory
[nt,nx,npe,SMS_shift,~] = size(kdata);
Phi=permute(Phi,[3 4 5 1 6 2 7]); 
Phi_B0=permute(Phi_B0,[3 4 5 1 6 2 7]); 
sens_map = permute(sens_map,[1 2 3 5 6 7 8 4]);
%% Phase evluation
Phase_T=zeros(nx,npe,ns,nt);
for t=1:1:nt
    Phase_T(:,:,:,t)=exp(1i*2*pi*P_dB*(TEs(t))).*Phase0_all(:,:,:,t);    % B0 phase of SE
end

kdata = permute(kdata,[2 3 8 1 6 7 4 5]);
mask_sample = kdata~=0;
Matrix_encode = Generate_SMS_encode_Matrix(nx,npe,ns,SMS_shift);

AHSigFN=[bart_savepath,'AHSig'];
SMS_FN = [bart_savepath,'SMS_encode'];
Phi_FN=[bart_savepath,'PhiUse'];
sample_mask_FN=[bart_savepath,'sample_mask'];
sens_FN=[bart_savepath,'sens_map'];
ISizeByCh_FN = [bart_savepath,'ISize'];
Ops={'name perCh','_Phi fmac 1 32','_SMS fmac 3 4','_FT fftc 3',...
    '_Mask fmac 2 0','a FT','a SMS','a Phi',...
    'nextlinop','setdefop','ident','normal',...
    '_OpenCh fmac 0 0','byslice perCh 7','a OpenCh'};
ScriptByCh=[bart_savepath, 'EPTI_CompsToSig_ByCh.txt'];
WriteLinopToFile(ScriptByCh,Ops);
writecfl(sens_FN,sens_map);
writecfl(SMS_FN,Matrix_encode);

%% B0 update
if recon_algorithm.B0update ~= 0

Echo_use_B0 = 1:size(Phi_B0,4);
ISize=FillOnesTo16_2([nx,npe,ns,1,1,K_B0]);
Phi_use = Phi_B0(:,:,:,Echo_use_B0,:,:).*Phase_T(:,:,:,Echo_use_B0);
writecfl(Phi_FN,Phi_use);
writecfl(sample_mask_FN,mask_sample(:,:,:,Echo_use_B0,:,:,:,1));

AHSig=sum(sum( sum(sum(ifft2c(kdata(:,:,:,Echo_use_B0,:,:,:,:)).*conj(sens_map),8) .*conj(Matrix_encode),7).* ...
    conj(Phi_use(:,:,:,Echo_use_B0,:,:)),4),5);
writecfl(AHSigFN,AHSig);
ISizeByCh=[ISize.' ISize.'];
writecfl(ISizeByCh_FN,ISizeByCh);

% Rec=bart(['picsS -g -S -d 0 -w ',num2str(scaling),' ' ScriptByCh],ISizeByCh_FN,AHSigFN,sens_FN,Phi_FN,sample_mask_FN,SMS_FN);
Rec=bart(['picsS ',gpu,' -S -d 0 -i 100 -u ',num2str(recon_algorithm.rho),' -C 8 -N  -m -w ',num2str(scaling),' -b ',num2str(recon_algorithm.block_size),' -R L:7:7:',num2str(recon_algorithm.lambda),...
    ' ' ScriptByCh],ISizeByCh_FN,AHSigFN,sens_FN,Phi_FN,sample_mask_FN,SMS_FN);

im_recon=squeeze(sum(Rec.*Phi_B0,6));   
if ns >1
    im_recon = permute(im_recon,[1 2 4 3]);
end
% wind=tukeywin2([npe,nx],0,16);
wind=tukeywin2([npe,nx],16,32);
im_recon = ifft2c(wind.*fft2c(im_recon));
PHS=angle(im_recon);
MSK_extended = P_dB~=0;
% P_dB0 = P_dB;
% P_dB_rec = [];
dB_global = [];
Echo_use_fit = 6:Echo_use_B0(end)-6;
for slice_calib = 1:ns 
    MSK = MSK_extended(:,:,slice_calib);
    tmp = dB_fitting_JumpCorrect(PHS(:,:,Echo_use_fit,slice_calib),TEs(Echo_use_fit(:)),logical(MSK),1);
    tmp(isnan(tmp))=0;
    tmp2 = tmp(MSK(:));
    dB_global(slice_calib) = mean(tmp2(abs(tmp2(:))<30));
    if isnan(dB_global(slice_calib)) 
        dB_global(slice_calib)=0;
    end
    if recon_algorithm.B0update == 2
        P_dB(:,:,slice_calib) = tmp+P_dB(:,:,slice_calib);
    elseif recon_algorithm.B0update == 1
        P_dB(:,:,slice_calib) = P_dB(:,:,slice_calib)+dB_global(slice_calib).*MSK;
    end
end
%         figure; imshow3((permute(P_dB_rec(:,end:-1:1,:),[2 1 3])),[-150 150],[1 ns]); colormap('jet');
%         figure; imshow3((permute(P_dB(:,end:-1:1,:),[2 1 3])),[-150 150],[1 ns]); colormap('jet');
%         figure; imshow3((permute(P_dB_rec(:,end:-1:1,:)-P_dB0(:,end:-1:1,:),[2 1 3])),[-50 50],[1 ns]); colormap('jet');
Phase_T=zeros(nx,npe,ns,nt);
for t=1:1:nt
    Phase_T(:,:,:,t)=exp(1i*2*pi*P_dB*(TEs(t))).*Phase0_all(:,:,:,t);    % B0 phase of SE
end
end
%%  image recon
ISize=FillOnesTo16_2([nx,npe,ns,1,1,K]);
Phi_use = Phi.*Phase_T;
writecfl(Phi_FN,Phi_use);
writecfl(sample_mask_FN,mask_sample(:,:,:,:,:,:,:,1));
AHSig=sum(sum( sum(sum(ifft2c(kdata).*conj(sens_map),8) .*conj(Matrix_encode),7).*conj(Phi_use),4),5);
writecfl(AHSigFN,AHSig);

ISizeByCh=[ISize.' ISize.'];
writecfl(ISizeByCh_FN,ISizeByCh);
if recon_algorithm.regularization == 1
    Rec=bart(['picsS ',gpu,' -S -d 0 -i 100 -u ',num2str(recon_algorithm.rho),' -C 8 -N -m -w ',num2str(scaling),' -b ',num2str(recon_algorithm.block_size),' -R L:7:7:',num2str(recon_algorithm.lambda),...
        ' ' ScriptByCh],ISizeByCh_FN,AHSigFN,sens_FN,Phi_FN,sample_mask_FN,SMS_FN);
else
    Rec=bart(['picsS ',gpu,' -S -d 0 -w ',num2str(scaling),' ' ScriptByCh],ISizeByCh_FN,AHSigFN,sens_FN,Phi_FN,sample_mask_FN,SMS_FN);
end
a = squeeze(Rec);
%% PF recon
im_recon=squeeze(sum(Rec.*Phi.*Phase_T,6));
if ns==1
    im_recon = permute(im_recon,[1 2 4 3]);
end
if PF_shift>0
    PF_factor = npe/(PF_shift+npe);
    k_recon = fft2c(im_recon);
    im_recon=[];
    k_recon_pf(:,PF_shift+1:npe+PF_shift,:,:) = k_recon;
    if recon_algorithm.PF_recontype == 0
        im_recon = ifft2c(k_recon_pf);
    else
        for sub_slice = 1:ns
            niter = 4;
            [ dp_rec ] = POCS_PF_EPTI( k_recon_pf(:,:,sub_slice,:),PF_factor,niter );
            im_recon(:,:,sub_slice,:) = ifft2c(dp_rec);
        end
    end
end

if (recon_algorithm.apodization == 1)
    im_recon = ifft2c(mrir_filter_raw_apodize_1d(mrir_filter_raw_apodize_1d(fft2c(im_recon), 1, recon_algorithm.apodize_factor),2,recon_algorithm.apodize_factor)) ;
end

im_recon = single(im_recon);
a = single(a);


