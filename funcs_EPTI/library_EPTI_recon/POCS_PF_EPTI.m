function [ dp_rec ] = POCS_PF_EPTI( data,PF_factor,niter )
% POCS partial Fourier reconstruction 
% Zijing Dong, 2018

% POCS partial Fourier reconstruction for k-t EPTI data
% Fuyixue Wang, 2021, MGH

[nx,ny,nt,nc]=size(data);

ny_r=round(ny.*PF_factor);

dp_rec=zeros(nx,ny,nt,nc);

wind=zpad(hamming(2*ny_r-ny),[ny 1]);
win=repmat(permute(wind,[2 1]),[nx,1,nt,nc]);

d_low=data.*win;

ms=ifft2c(d_low);
clear d_low;
clear win;
clear win2;
mask=(data~=0);

for i=1:niter
    dp_comb=dp_rec;
    dp_comb(mask)=data(mask);
    im_cmb=abs(ifft2c(dp_comb));
    im_rec=im_cmb.*(ms./abs(ms));
    dp_rec=fft2c(im_rec);
end

end

