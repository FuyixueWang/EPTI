function [ kernel ] = calibration_GRAPPA( kdata,R,ksize,lambda)
% Calculate GRAPPA kernel for GRAPPA reconstruction of the accelerated EPTI 
% calibration scan. Calibration used to be fully-sampled k-t space, now 
% accelerated along PE to reduce the scan time
% Fuyixue Wang, 2016, applied to EPTI calibration, 2022

if nargin < 4
	lambda = 0;
end
%% calib
  nc=size(kdata,3);
  ksize(1)=(ksize(1)-1)*R+1;
  mask=zeros(ksize(1),ksize(2),nc);
  mask([1 ksize(1)],:,:)=1;
  mask=mask==1;
%    figure; imshow(mask(:,:,1),[]);

  mask2=zeros(ksize(1),ksize(2),nc);
  mask2(2:end-1,round(ksize(2)/2),:)=1;
  mask2=mask2==1;
%    figure; imshow(mask2(:,:,1),[]);
%   imshow(sos(abs(kdata)),[]);
  calib=[];
  for i=1:nc
      calib(:,:,i)=transpose(im2col(squeeze(kdata(:,:,i)),[ksize(1) ksize(2)]));  % size: (blocks,kernel size,nc)
  end
  calib=calib(1:end,:);        % size: (blocks,kernel size*nc)
%   imshow(abs(calib),[]);
  cali_size=size(calib);
  Y=reshape(calib(repmat(mask2(:)',[size(calib,1) 1])),[cali_size(1) (R-1)*nc]);
  A=reshape(calib(repmat(mask(:)',[size(calib,1) 1])),[cali_size(1) cali_size(2)/ksize(1)*2]);
  AtA=A'*A;
  lambda = norm(AtA,'fro')/size(AtA,1)*lambda;
  if lambda==0
        kernel=AtA\A'*Y;
  else
      kernel=(AtA+ eye(size(AtA))*lambda)\A'*Y;
  end
end

