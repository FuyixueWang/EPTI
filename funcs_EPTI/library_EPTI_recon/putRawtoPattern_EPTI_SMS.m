function [kdata] = putRawtoPattern_EPTI_SMS(kdata,Nseg,Rseg,Rpe,MB)
% put raw k-space EPTI data into EPTI encoding pattern
% change the Ngroup to be the number of echoes instead of echo group
% to allow arbitrary echo # for subspace reconstruction
% Fuyixue Wang 2018

% support EPTI with SMS
% Fuyixue Wang 2019
%% Generate Trajectory
[N_total,nx,~,nc,nslice]=size(kdata);
dk_pe = 1;
N_blipdown = Rseg/Rpe;
dk_blipdown = dk_pe*Rpe;

dk_pre = dk_pe*Rseg;

ky_acq = zeros(Nseg,N_total);
ky_pre=zeros(1,Nseg);
blip = zeros(Nseg,N_total+1);

for diff_s = 1:Nseg
    if(mod(diff_s,2)==1) 
        ky_pre(diff_s) = diff_s*dk_pre;
        for diff_n = 1:N_total
            if( mod(diff_n,2*N_blipdown) == N_blipdown )
                blip(diff_s,diff_n+1)=(N_blipdown-1)*dk_blipdown-dk_pe*Rpe/2;
            elseif( mod(diff_n,2*N_blipdown) == 0)
                blip(diff_s,diff_n+1)=(N_blipdown-1)*dk_blipdown+dk_pe*Rpe/2;
            else
                blip(diff_s,diff_n+1) = -dk_blipdown;
            end
            ky_acq(diff_s,diff_n) = ky_pre(diff_s)+squeeze(sum(blip(diff_s,1:diff_n),2));
        end
    else
        ky_pre(diff_s) = (diff_s)*dk_pre-dk_pe*Rpe/2;
        for diff_n = 1:N_total
            if( mod(diff_n,2*N_blipdown) == N_blipdown )
                blip(diff_s,diff_n+1)=(N_blipdown-1)*dk_blipdown+dk_pe*Rpe/2;
            elseif( mod(diff_n,2*N_blipdown) == 0)
                blip(diff_s,diff_n+1)=(N_blipdown-1)*dk_blipdown-dk_pe*Rpe/2;
            else
                blip(diff_s,diff_n+1) = -dk_blipdown;
            end
            ky_acq(diff_s,diff_n) = ky_pre(diff_s)+squeeze(sum(blip(diff_s,1:diff_n),2));
        end
    end
end
%%  put data into pattern
npe=Rseg*Nseg;
kdata_EPTI=zeros(N_total,nx,npe,nc,nslice);

for shot=1:Nseg
    for nt=1:N_total
        kdata_EPTI(nt,:,ky_acq(shot,nt),:,:)=kdata(nt,:,shot,:,:);
    end
end
kdata_EPTI=kdata_EPTI(:,:,end:-1:1,:,:);

if MB>1
    kdata=zeros(N_total,nx,npe,MB,nc,nslice);
    for i = 1:MB
        kdata(i:MB:end,:,:,i,:,:)=kdata_EPTI(i:MB:end,:,:,:,:);
    end
else
    kdata=kdata_EPTI;
end

end


