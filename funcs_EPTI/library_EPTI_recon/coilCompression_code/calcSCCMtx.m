function mtx = calcSCCMtx(calibDATA)
% mtx = calcGCC(calibDATA )
% Coil Compression using a single compression matrix. Based on Huang et.
% al, MRM 2008;26:133-141
% The function computes and returns compression matrix 
%
%
%  Inputs:
%           calibDATA - a 4D matrix representing [Kx Ky Kz Coils] calibration data
%                   or a 3D matrix representing [Kx Ky COils]
%           
%
%  Outputs:
%           mtx - the compression matrix.
%
% See:
%       calcECCMtx, ECC, CC, alignCCMtx, calcGCCMtx
%
% (c) Michael Lustig 2013


if nargin < 3
    ws = 1;
end


% check if k-space is 2D or 3D 
if length(size(calibDATA))==3
    ndims = 2;
    calibDATA = permute(calibDATA,[1,2,4,3]);
else
    ndims = 3;
end


[Nx,Ny,Nz,Nc] = size(calibDATA);
calibDATA = reshape(calibDATA,[Nx*Ny*Nz,Nc]);

[U,S,V] = svd(calibDATA,'econ');

if size(calibDATA,1) > size(calibDATA,2)
    mtx = V;
else
    mtx = U;
end
