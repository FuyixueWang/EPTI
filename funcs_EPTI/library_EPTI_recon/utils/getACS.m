function [ ACS ] = getACS( kData,ACSsize )
%get ACSdata from full kspace data for kernel
[nfe,npe,nch]=size(kData);
ACS{1}=round(nfe/2)-round(ACSsize(1)/2)+1:round(nfe/2)-round(ACSsize(1)/2)+ACSsize(1);
ACS{2}=round(npe/2)-round(ACSsize(2)/2)+1:round(npe/2)-round(ACSsize(2)/2)+ACSsize(2);
ACS{3}=1:nch;

end

