% method 1
% navigator phase conjugation
function [ ksp_nav_new ] = Phase_inverse( ksp_nav )
[Nfe_nav,Npe_nav,Nc,Nex,Nshot]=size(ksp_nav);
im_nav = ifft2c(ksp_nav);

im_nav_comb = zeros(Nfe_nav,Npe_nav);
ksp_nav_new = zeros(size(ksp_nav));
for nshot = 1:Nshot
    for nex = 1:Nex
    im_nav_comb = sos(im_nav(:,:,:,nex,nshot),3).*exp(-1i*sop(im_nav(:,:,:,nex,nshot),3)); % channel combine
    s_nav0=im_nav(:,:,:,nex,nshot)./repmat((sos(im_nav(:,:,:,nex,nshot),3).*exp(1i*sop(im_nav(:,:,:,nex,nshot),3))),[1 1 Nc]); 
    ksp_nav_new(:,:,:,nex,nshot) = fft2c(s_nav0.*repmat((im_nav_comb),[1 1 Nc]));% image conjugation
    end
end

% method 2
% ksp_nav_new = zeros(size(ksp_nav));
% for nshot = 1:Nshot
%     for nex = 1:Nex
%         ksp_nav_new(:,:,:,nex,nshot) = fft2c(conj(im_nav(:,:,:,nex,nshot)./im_nav0).*im_nav0);
%     end
% end