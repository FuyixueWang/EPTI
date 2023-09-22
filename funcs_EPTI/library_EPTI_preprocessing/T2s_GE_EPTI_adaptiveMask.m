function [S0, T2s, fitting_map, mask_adapt] = T2s_GE_EPTI_adaptiveMask(imgs,TEs,mask,threshold)

% Fit T2* from multi-echo GE-EPTI data
% Fuyixue Wang, MGH, 2021
% Inputs
% imgs:          input images at each echo time 
% TEs:           echo times 
% mask:          mask of voxels to fit 

[nx,ny,np]=size(imgs);

A1 = [ones(np,1),-TEs];   % ln(s0) 
% A1 = [ones(np,1),-TEs(1:np)];   % ln(s0) 
A = A1;

S0=zeros(nx,ny);
T2s=zeros(nx,ny);
mask_adapt = zeros(nx,ny);
fitting_map=zeros(nx,ny,np);
% rmse=zeros(nx,ny);
for ii=1:nx
    for jj=1:ny
        if mask(ii,jj) == 1
           % fit voxel
           signal=squeeze(imgs(ii,jj,:));
           maskk=signal>(threshold*signal(1));
           mask_adapt(ii,jj) = sum(maskk(:));
           signal=signal(maskk);
           AA=A(repmat(maskk(:),[1 2]));
           AA=reshape(AA,[size(signal,1),2]);
           param = AA\log(signal);
           %  result  
           S0(ii,jj) = exp(param(1));
           T2s(ii,jj) = 1/param(2);
           fitting_map(ii,jj,:)=S0(ii,jj)*exp(-TEs*param(2));
        end
    end % end jj 
end % end ii
S0(isnan(S0)) = 0;
T2s(isnan(T2s)) = 0; 

