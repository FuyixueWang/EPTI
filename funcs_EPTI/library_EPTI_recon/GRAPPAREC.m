function [ recon_k ] = GRAPPAREC( kdata,kernel,R,ksize )
% GRAPPA reconstruction
% undersampling must along the first dimension
% Fuyixue Wang, 2016

% Apply to EPTI calibration scan accelerated along PE to reduce the scan time
% Fuyixue Wang, 2022

[nx ny nc]=size(kdata);
ksize(1)=(ksize(1)-1)*R+1;
mask=zeros(ksize(1),ksize(2),nc);
mask([1 ksize(1)],:,:)=1;
mask=mask==1;
recon_k=kdata;
maskk=kdata~=0;
% imshow(maskk(:,:,1),[])
for i=1:R:nx-R
  A=[];
  for j=1:nc
     A(:,:,j)=transpose(im2col(squeeze(kdata(i:R+i,:,j)),[ksize(1) ksize(2)]));
  end
  A=A(1:end,:);
  cali_size=size(A);
  A=reshape(A(repmat(mask(:)',[size(A,1) 1])),[cali_size(1) cali_size(2)/ksize(1)*2]);  
  Y=(A*kernel);
  recon_k(i+1:i+R-1,1+(ksize(2)-1)/2:ny-((ksize(2)-1)/2),:)=permute(reshape(Y,[ny-(ksize(2)-1) R-1  nc]),[2 1 3]);
end
recon_k=recon_k.*(1-maskk)+kdata;
end