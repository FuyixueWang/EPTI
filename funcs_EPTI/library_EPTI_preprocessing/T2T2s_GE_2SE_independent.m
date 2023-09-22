function [S0,S02,T2,T2s,fitting_map_GE,fitting_map_SE,mask_adapt_T2s] = T2T2s_GE_2SE_independent(imgs,TEs_GE,TEs_SE,te_SE,mask,threshold)

%  fit T2 and T2* 
% Inputs
% imgs:          input images at each echo time 
% TEs:            echo times 
% te2:            time of spin-echo 
% mask:          mask of voxels to fit 
% np_SE1:        time point index of the spin echo from the SE series
% Outputs
% So T2 T2s rmse : three parameter fit and rmse 
% Fuyixue Wang, MGH

[nx,ny,~]=size(imgs);
np_GE=length(TEs_GE);
np_SE1 = find(TEs_SE>te_SE,1);
np_SE2=length(TEs_SE)-np_SE1;

[S0, T2s, fitting_map_GE, mask_adapt_T2s] = T2s_GE_EPTI_adaptiveMask(imgs(:,:,1:np_GE),TEs_GE,mask,threshold);

im_SE=imgs(:,:,np_GE+1:end);
R2s=1./T2s;
R2s(isnan(R2s))=0;
A1 = [ones(np_SE1,1), te_SE-2*TEs_SE(1:np_SE1)];   % ln(s0)  R2
A2 = [ones(np_SE2,1), -ones(np_SE2,1)*te_SE];  % ln (s02)  R2
A = [A1;A2];

S02=zeros(nx,ny);
T2=zeros(nx,ny);
fitting_map_SE=zeros(nx,ny,np_SE1+np_SE2);

for ii=1:nx
    for jj=1:ny
        if mask(ii,jj) == 1
           % fit voxel
           signal=squeeze(im_SE(ii,jj,:));
           maskk=signal>threshold;
           signal=signal(maskk);
           T2s_component=exp([(TEs_SE(1:np_SE1)-te_SE);te_SE-(TEs_SE(np_SE1+1:end))]*R2s(ii,jj));
           signal=signal./T2s_component;
           AA=A(repmat(maskk(:),[1 2]));
           AA=reshape(AA,[size(signal,1),2]);
           param = AA\log(signal);
%          param = lsqr(AA,log(signal),[],[],[],[],50*ones(4,1));
           %  result  
            S02(ii,jj) = exp(param(1));
            T2(ii,jj) = 1/param(2);
            fitting_map_SE(ii,jj,1:np_SE1)=T2s_component(1:1:np_SE1).*S02(ii,jj).*exp((te_SE-2*TEs_SE(1:np_SE1))*param(2));
            fitting_map_SE(ii,jj,np_SE1+1:end)=T2s_component(np_SE1+1:end).*S02(ii,jj).*exp(-te_SE*param(2));
        end
    end % end jj 
end % end ii

S02(isnan(S02)) = 0;
S0(isnan(S0)) = 0;
T2(isnan(T2)) = 0;
T2s(isnan(T2s)) = 0; 


