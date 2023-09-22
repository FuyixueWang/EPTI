function im_shifted = image2dshift(im,shift_x,shift_y)
% shift the image in the first and second dimension by shift_x and shift_y
    nsize = size(im);
    if length(nsize) < 3 
        nsize(3)=1;
    end
    nx = nsize(1);
    ny = nsize(2);
    tmp_y=exp(-1i*2*pi* shift_y/ny*( (1:1:ny)-(ny+1)/2 )');
    tmp_y=permute(repmat(tmp_y,[1,nx,nsize(3:end)]),[2,1,3:length(nsize)]);
    tmp_x=exp(-1i*2*pi* shift_x/nx*( (1:1:nx)-(nx+1)/2 )');
    tmp_x=repmat(tmp_x,[1,ny,nsize(3:end)]);
    im_shifted = ifft2c(tmp_x.*tmp_y.*fft2c(im));

end

