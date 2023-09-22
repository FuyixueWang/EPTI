function [S0, T2, T2s, fitting_map] = T2T2s_SE_EPTI(imgs,TEs,te2,mask,threshold)

% Fit T2 and T2* from SE EPTI readout
% Zijing Dong, MGH, 2021

% Inputs
% imgs:          input images at each echo time 
% TEs:            echo times 
% te2:            time of spin-echo 
% mask:          mask of voxels to fit 
% displayON;    enables display 
% Outputs
% So T2 T2s rmse : three parameter fit and rmse 
[nx,ny,np]=size(imgs);

ind = find(TEs>te2,1);
A1 = [ones(ind-1,1), TEs(1:ind-1)-te2, -TEs(1:ind-1)];   % ln(s0) R' R2
A2 = [ones(np-ind+1,1),te2-TEs(ind:np), -TEs(ind:np)];
A = [A1; A2];

S0=zeros(nx,ny);
T2=zeros(nx,ny);
T2s=zeros(nx,ny);
fitting_map=zeros(nx,ny,np);
% rmse=zeros(nx,ny);
for ii=1:nx
    for jj=1:ny
        if mask(ii,jj) == 1
           % fit voxel
           signal=squeeze(imgs(ii,jj,:));
           maskk=signal>threshold;
           signal=signal(maskk);
           AA=A(repmat(maskk(:),[1 3]));
           AA=reshape(AA,[size(signal,1),3]);
           param = AA\log(signal);
           %  result  
            S0(ii,jj) = exp(param(1));
            T2(ii,jj) = 1/param(3);
            T2s(ii,jj) = 1/(param(2)+param(3));
            fitting_map(ii,jj,:)=S0(ii,jj)*[exp(-te2*param(2)).*exp(-TEs(1:ind-1)*(param(3)-param(2)));exp(te2*param(2)).*exp(-TEs(ind:np)*(param(3)+param(2)))];
        end
    end % end jj 
end % end ii
S0(isnan(S0)) = 0;
T2(isnan(T2)) = 0;
T2s(isnan(T2s)) = 0; 

