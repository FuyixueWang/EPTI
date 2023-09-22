function mtx = calcECCMtx(DATA,dim,ncc, ks)
% mtx = ECC(DATA,[dim,ncc, ks)
%
% ESPIRiT based coils compression.
% Based on D. Bahri et. al, "ESPIRiT-Based Coil Compression for Cartesian
% Sampling" ISMRM 2013 pp:2657
% The function computes and returns compression matrices for each position in space
% along a chosen dimension. To compress, only use the Nth firs entries of the 2st
% dimension of mtx. Use with alignCCMtx to align and CC to perform compression
% ECC is similar to GCC. It is slower, but is less biased by noise
% than GCC.
%
%
%
%  Inputs:
%           DATA - a 4D matrix representing [Kx Ky Kz Coils] data
%                   or a 3D matrix representing [Kx Ky COils]
%           dim  - dimension to perform compression on
%           ncc  - approxmate number of target of virtual coils. This is
%           useful to set the right thresholding parameters. If missing
%           then going to estimate based on Singular Values. 
%           ks   - ESPIRiT kernel size (default 8)
%
%  Outputs:
%
%           mtx -  compression matrix.
%
% See:
%       calcGCCMtx, CC, alignCCMtx
%
% (c) Michael Lustig 2013

if nargin < 3
    ncc = [];
end


if nargin  < 4
    ks = [8];
end


% check if 2D or 3D data
if length(size(DATA))==3
    ndims = 2;
    DATA = permute(DATA,[1,2,4,3]);
else
    ndims = 3;
end

% permute data if compression is done on 2 dimension
% done for simple reusible code.
if dim==2
    DATA = permute(DATA,[2,1,3,4]);
end

if dim==3
    DATA = permute(DATA,[3,1,2,4]);
end

if dim>3 | dim <1
    disp('Error, compression dimension is 1 2 or 3')
end

% perform IFFT
[Nx,Ny,Nz,Nc] = size(DATA);

DATA = reshape(DATA,[Nx,Ny*Nz,Nc]);

% Construct an ESPIRiT kernel size 8x1 -- this is enough
% for most coils arrays.
[k,S] = dat2Kernel(DATA,[ks,1]);

if isempty(ncc)
    idx = min(find(S<S(1)*0.02));
else 
    idx = min(ceil(length(S)*ncc*1.2/Nc),length(S));
end

    
% crop kernels
k = k(:,:,:,1:idx);


% Compute compression matrixes by eigen decomposition in image space.
[mtx,W] = kernelEig(k,[Nx,1]);

% compression matrices must be conjugated -- compatability with ESPIRiT
% SENSE code.
mtx = conj(permute(mtx,[3,4,1,2]));
% set most important coils to be first
mtx = mtx(:,end:-1:1,:);
