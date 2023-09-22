function [ dB,fit_error,phase ] = dB_fitting_JumpCorrect( phase,TEs,mask,unwrap_flag)
% Fitting of the B0 inhomogeneity field (dB) using multi-echo phase data
% Fuyixue Wang, 2018

% updated, clean up
% Zijing Dong, 2019

if(unwrap_flag==1)
    phase=unwrap(phase,[],3);
end
phase=phase./(2*pi);
[nx,ny,nt]=size(phase);
A = [ones(nt,1),TEs(1:nt)];   
dB=zeros(nx,ny);
fit_error=zeros(nx,ny);
for ii=1:nx
    for jj=1:ny
        if mask(ii,jj) == 1
           % fit voxel
           signal=double(squeeze(phase(ii,jj,:)));
%            signal_diff=signal(2:end)-signal(1:end-1);
%            Jump_index=find(abs(signal_diff)>0.3);
%            if (~isempty(Jump_index))&&(length(Jump_index)<3)               
%            end
           param = A\signal;
           phase_fitting=A*param;
           error = abs(phase_fitting-signal);
%            mask_temp=(error > thershold*mean(error(:)));
%            jump_num=sum(mask_temp(:));
%            if (jump_num>0)&&(jump_num<2) 
%                mask_temp=logical(1-mask_temp);
%                param = A(mask_temp,:)\signal(mask_temp);
%                phase_fitting=A(mask_temp,:)*param;
%                error = abs(phase_fitting-signal(mask_temp));
%            end
           dB(ii,jj) = param(2);
           fit_error(ii,jj)=sqrt(sum(error.^2)/size(error,1));   
        end
    end % end jj 
end % end ii
dB(isnan(dB)) = 0;
% dB=dB-mean(dB(mask(:)));
dB = dB.*mask;
phase=phase.*(2*pi);
end

