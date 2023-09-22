function [im_mean imgtemp]=cpxaver_mxd(img,winsize)
if ~exist('winsize','var')
    winsize = 1/3;
end
    [nPE,nFE,n1,n2,n3]=size(img);
    n=n1*n2*n3;
    img=reshape(img,[nPE,nFE,n]);
    window=zeros(nPE,nFE);
%     window(nPE/2+1-nPE/4:nPE/2+nPE/4,nPE/2+1-nPE/4:nPE/2+nPE/4)=triang(nPE/2)*triang(nPE/2)';
%     window(nPE/2+1-nPE/8:nPE/2+nPE/8,nPE/2+1-nPE/8:nPE/2+nPE/8)=triang(nPE/4)*triang(nPE/4)';      % More aliasing
    warning off;
%     window(nPE/2+1-nPE/3:nPE/2+nPE/3,nPE/2+1-nPE/3:nPE/2+nPE/3)=triang(nPE*2/3)*triang(nPE*2/3)';  % More noise
%     window(nPE/2+1-nPE*1/2*winsize:nPE/2+nPE*1/2*winsize,nPE/2+1-nPE*1/2*winsize:nPE/2+nPE*1/2*winsize)=triang(nPE*winsize)*triang(nPE*winsize)';
    window=crop(generate_hamm(round(winsize*nPE),nPE),nPE,nFE);

%     opts_phi.TVtype=2;
% opts_phi.aTV=0.1;
imgtemp = zeros(nPE,nFE,n);
    for i=1:n
        imgtemp(:,:,i)=ifft2c(fft2c(img(:,:,i)).*window);

% imgtemp(:,:,i)=TVDenoise(img(:,:,i),opts_phi);
%         imgtemp=imfilter(img(:,:,i),fspecial('gauss'));

        img(:,:,i)=img(:,:,i).*conj(imgtemp(:,:,i)./abs(imgtemp(:,:,i)));
        
%                     pha=medfilt2(angle(img(:,:,i)),[9,9]);
%         imgtemp = abs(img(:,:,i)).*exp(1i*pha);  
%         img(:,:,i)=img(:,:,i).*conj(imgtemp./abs(imgtemp));
        
    end
    masknan=isnan(img);
    img(masknan)=0;
    im_mean=mean(img,3);
end