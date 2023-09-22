function res = sop(x ,dim)
% res = sop(x ,dim)
%
% function computes the phase of sum of complex data along dimension dim.
% If dim is not specified, it computes it along the last dimension.
%
% (c) Feiyu Chen 2012

if nargin < 2
    dim = size(size(x),2);
    disp('Must check the dim of coil')
end

diminput = size(size(x),2);

coil=size(x,dim);
y=squeeze(sum(abs(x),dim));
cpx=y.*0;
for ic=1:coil
    if dim == diminput
        if dim == 3
            xcoil = squeeze(x(:,:,ic));
        end
        if dim == 4
            xcoil = squeeze(x(:,:,:,ic));
        end
        if dim == 5
            xcoil = squeeze(x(:,:,:,:,ic));
        end
    end
    if dim == diminput - 1
        if dim == 3
            xcoil = squeeze(x(:,:,ic,:));
        end
        if dim == 4
            xcoil = squeeze(x(:,:,:,ic,:));
        end
        if dim == 5
            xcoil = squeeze(x(:,:,:,:,ic,:));
        end
    end
    cpx = xcoil.*abs(xcoil)./y + cpx; 
%     xcoil=squeeze(x(:,:,ic));
%     ang=angle(xcoil).*(abs(xcoil)./y)+ang;
end
res=atan2(imag(cpx),real(cpx));                                
% res = angle(sum(x,dim));
