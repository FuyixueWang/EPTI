function [kdata] = putRawtoPattern_SMS_zigzag_PF_selfnav_oldseq(kdata,Nseg,Rseg,Rpe,MB,PFshift)
% put raw k-space EPTI data into EPTI encoding pattern
% change the Ngroup to be the number of echoes instead of echo group
% to allow arbitrary echo # for subspace reconstruction
% Fuyixue Wang 2018

% support EPTI with SMS
% Fuyixue Wang 2019

% put into pattern of zigzag single-shot EPTI data
% Fuyixue Wang 2021, MGH

% works for partial Fourier acquisition
% Zijing Dong 2022, MGH
%% Generate Trajectory
[nt,nx,~,nc,nslice]=size(kdata);
npe = Nseg*Rseg;
ky_acq = zeros(Nseg,nt);
Necho_each = round(Rseg/Rpe/2);
Necho_PFshift = floor(PFshift/Rpe/2);
if Nseg>1
    Start_point = (Rseg:Rseg:Nseg*Rseg) - Rseg/2;
else
    Start_point = (Rseg+PFshift)/2 - PFshift;
end
change_point = (1:2:50)*Necho_each-Necho_PFshift;
for shot_num = 1:Nseg
    PE_currect = Start_point(shot_num);
    Sign_PEdirection = -1;
    for i = 1:nt
      ky_acq(shot_num,i) = PE_currect;
      if sum(i==change_point(:))>0
          PE_currect = PE_currect+Sign_PEdirection*Rpe/2;
          Sign_PEdirection = -1*Sign_PEdirection;
      else
          PE_currect = PE_currect+Sign_PEdirection*Rpe;
      end
    end
end
%%  put data into pattern
kdata_EPTI=zeros(nt,nx,npe,nc,nslice);

for shot=1:Nseg
    for t=1:nt
        kdata_EPTI(t,:,ky_acq(shot,t),:,:)=kdata(t,:,shot,:,:);
    end
end
if MB>1
    kdata=zeros(t,nx,npe,MB,nc,nslice);
    for i = 1:MB
        kdata(i:MB:end,:,:,i,:,:)=kdata_EPTI(i:MB:end,:,:,:,:);
    end
else
    kdata=kdata_EPTI;
end

end


