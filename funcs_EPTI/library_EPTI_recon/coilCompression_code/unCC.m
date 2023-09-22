function [res] = unCC(cDATA,dim,mtx)
% [res] = unCC(cDATA,dim,mtx)
%
% Geometric Decomposition Coil Compression Based on: Zhang et. al MRM 2013;69(2):571-82.
% The function uses compression matrices computed by calcECCMtx or
% calcGCCMtx to uncompress compressed virtual coils and rotate back to the
% original physical coils. This is to mainly check the error of compression
% or for processing tht requires the original non-compressed channels
%
%  Inputs:
%           cDATA - a 4D matrix representing [Kx Ky Kz cCoils] data
%                   or a 3D matrix representing [Kx Ky cCOils]
%           dim  - dimension to perform compression on
%           mtx  - Aligned compression matrices from calicECCMtx,
%           calcGCCMtx
%
%  Outputs:
%           res - returned data that is rotated to the compressed space.
%                 The first ncc coils are the compressed & aligned virtual coils.
%
%
% See:
%       calcGCCMtx,  calcECCMtx, alignCCMtx, CC
%
%
% Example Use compressing from 8->4 channels along 2nd dimension:
%
%       load brain_8ch
%       dim = 2;
%       calib = crop(DATA,[24,200,8]);
%       gccmtx =  calcGCCMtx(calib,dim,1);
%       gccmtx = gccmtx(:,1:4,:);
%       gccmtx = alignCCMtx(gccmtx);
%       CCDATA = CC(DATA,dim,gccmtx);
%       UCCDATA = unCC(CCDATA,dim,gccmtx);
%       figure, imshow(sos(ifft2c(DATA-UCCDATA)).^(1/2),[]);
%       title('compression error from 8->4 channels')
%
% (c) Michael Lustig 2013

if nargin < 4
    pr = 0;
end

% check if 2D or 3D data
if length(size(cDATA))==3
    ndims = 2;
    cDATA = permute(cDATA,[1,2,4,3]);
else
    ndims = 3;
end

Nc = size(mtx,1);


% permute data if compression is done on 2 dimension
% done for simple reusible code.
if dim==2
    cDATA = permute(cDATA,[2,1,3,4]);
end

if dim==3
    cDATA = permute(cDATA,[3,1,2,4]);
end

if dim>3 | dim <1
    disp('Error, compression dimension is 1 2 or 3')
end

% perform IFFT
[Nx,Ny,Nz,ncc] = size(cDATA);
im = ifftc(cDATA,1);
res = zeros(Nx,Ny,Nz,Nc);

% rotate data by compression matrices
for n=1:Nx
    tmpc = reshape(squeeze(im(n,:,:,:)),[Ny*Nz,ncc]);
    res(n,:,:,:) = permute(reshape(tmpc*mtx(:,:,n)',[Ny,Nz,Nc]),[4,1,2,3]);
end

% perform fft
res = fftc(res,1);

% permute if necessary
if dim==2
    res = permute(res,[2,1,3,4]);
end

if dim==3
    res = permute(res,[3,1,2,4]);
end

% squeeze back to two dimensions
if ndims ==2
    res = squeeze(res);
end

