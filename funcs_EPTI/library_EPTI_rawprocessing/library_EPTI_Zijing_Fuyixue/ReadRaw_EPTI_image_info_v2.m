function meas = ReadRaw_EPTI_image_info_v2(filename)
% EPTI rawdata processing for imaging data
% quickly read rawdata information (skip reading the whole data)
% Fuyixue Wang, MGH, 2020

% Fuyixue Wang, Zijing Dong, optimized for C2P, MGH, 2023
%%
nRepToRead=1;
BeginRep=1;
ascendingSlice_acq=0;
[meas] = read_meas_dat_memmap_EPTI_SMS_Generalized(filename,nRepToRead,BeginRep,ascendingSlice_acq);

meas.prot.Nseg = meas.prot.sWipMemBlock_alFree(34);
meas.prot.Rseg=meas.prot.sWipMemBlock_alFree(33);
meas.prot.Rpe=meas.prot.sWipMemBlock_alFree(35);
meas.prot.NslicesEX=meas.prot.sWipMemBlock_adFree(1);
meas.prot.PF_shift = meas.prot.sWipMemBlock_alFree(7);
meas.prot.N_slice = size(meas.data,10);

if meas.prot.SMSdata==1
    meas.prot.SMS_shift = meas.prot.sWipMemBlock_adFree(3);
else
    meas.prot.SMS_shift = 1;
end
if meas.prot.pf_echo > 0
    meas.prot.nechoSE=size(meas.data,2);
    meas.prot.nechoGE=floor(meas.prot.pf_echo*meas.prot.nechoSE);
else
    meas.prot.nechoGE=size(meas.data,2);
end
meas.prot.N_dyn=(meas.prot.lRepetitions+1)/meas.prot.Nseg;

meas.data =[];
meas.data_phascor1d = []; 
