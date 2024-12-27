function [im_recon,a] = SubspaceRecon_MATLAB(kdata,sens_map,P_dB,Phase0_odd,Phase0_even,gccmtx_aligned,TEs,Phi,parameters,scaling,recon_algorithm)
% EPTI subspace reconstruction using MATLAB 

% Please cite the following work on EPTI for this, including: 
% 1. Wang F, Dong Z, Reese TG, et al. Echo planar time-resolved imaging (EPTI). Magn Reson Med. 2019 Jun;81(6):3599-3615. doi: 10.1002/mrm.27673.
% 2. Dong Z, Wang F, Reese TG, et al. Echo planar time-resolved imaging with subspace reconstruction and optimized spatiotemporal encoding. Magn Reson Med. 2020 Nov;84(5):2442-2455. doi: 10.1002/mrm.28295.
% 3. Dong Z, Wald LL, Polimeni JR, Wang F. Single-shot echo planar time-resolved imaging for multi-echo functional MRI and distortion-free diffusion imaging. Magn Reson Med. 2024 Oct 20. doi: 10.1002/mrm.30327.
% Other relevant citations that may warrant consideration (EPTI's application in fMRI, dMRI and qMRI), see EPTI website: https://martinos.org/~fw089/

% Fuyixue Wang, Zijing Dong, 2022, MGH

% modifying for generalized EPTI sequences (unfinished)
%% set parameters
PF_shift = parameters.PF_shift;
ns = parameters.NslicesEX;
K=size(Phi,2);
coilcompresson_on = 1;
dim = 1; % dimension of readout for coil compression
%% Generate Trajectory
[nt,nx,npe,SMS_shift,~] = size(kdata);
ncc=size(sens_map,4);
if coilcompresson_on == 1
    kdata_cc=zeros(nt,nx,npe,SMS_shift,ncc);
    for t=1:nt
        CCDATA_aligned = CC(permute(kdata(t,:,:,:,:),[2 3 4 5 1]),gccmtx_aligned, dim);
        kdata_cc(t,:,:,:,:)=CCDATA_aligned;
    end
    kdata=kdata_cc;
    clear kdata_cc;
end
%% Phase evluation
Phase_T=zeros(nx,npe,ns,nt);
for t=1:2:nt
    Phase_T(:,:,:,t)=exp(1i*2*pi*P_dB*(TEs(t))).*Phase0_odd;    % B0 phase of SE
end
for t=2:2:nt
    Phase_T(:,:,:,t)=exp(1i*2*pi*P_dB*(TEs(t))).*Phase0_even;    % B0 phase of SE
end
%%  image recon
mask_sample = kdata~=0;
Matrix_encode = Generate_SMS_encode_Matrix(nx,npe,ns,SMS_shift);
if recon_algorithm.regularization == 1
    iter_ops.max_iter = 30;
    llr_ops.lambda = .005;
    iter_ops.rho = 0.1;
    llr_ops.block_dim = [8, 4];
    lsqr_ops.max_iter = 10;
    lsqr_ops.tol = 1e-4;
    [im_recon,a]=EPTI_Image_Recon_Subspace_GESE_LLR_SMS(kdata,mask_sample,...
        sens_map,Phase_T,Phi,Matrix_encode,iter_ops,llr_ops,lsqr_ops,scaling);
else
    a0=zeros([nx,npe,K]);
    lambda = 0.01;   % tikhonov regularization
    N_iteration = 60;
    [im_recon,a]=EPTI_Image_Recon_Subspace_GESE_Tik_hybrid(kdata,mask_sample,sens_map,Phase_T,Phi,a0,N_iteration,lambda);
end
%% PF recon
% im_recon=squeeze(sum(Rec.*Phi.*Phase_T,6));
if ns==1
    im_recon = permute(im_recon,[1 2 4 3]);
end
if PF_shift>0
    PF_factor = npe/(PF_shift+npe);
    k_recon = fft2c(im_recon);
    im_recon=[];
    k_recon_pf(:,PF_shift+1:npe+PF_shift,:,:) = k_recon;
    if recon_algorithm.PF_recon == 0
        im_recon = ifft2c(k_recon_pf);
    else
        for sub_slice = 1:ns
            niter = 4;
            [ dp_rec ] = POCS_PF_EPTI( k_recon_pf(:,:,sub_slice,:),PF_factor,niter );
            im_recon(:,:,sub_slice,:) = ifft2c(dp_rec);
        end
    end
end

im_recon = single(im_recon);
a = single(a);


