function [parameters,meas] = ReadRaw_EPTI_image_info(filename)
% EPTI rawdata processing for imaging data
% quickly read rawdata information (skip reading the whole data)
% Fuyixue Wang, MGH, 2020

% Fuyixue Wang, Zijing Dong, optimized for C2P, MGH, 2023
%%
nRepToRead=1;
BeginRep=1;
ascendingSlice_acq=0;
[meas] = read_meas_dat_memmap_EPTI_SMS_Generalized(filename,nRepToRead,BeginRep,ascendingSlice_acq);
parameters=meas.prot;
parameters.Nseg = meas.prot.sWipMemBlock_alFree(34);
parameters.Rseg=meas.prot.sWipMemBlock_alFree(33);
parameters.Rpe=meas.prot.sWipMemBlock_alFree(35);
parameters.NslicesEX=meas.prot.sWipMemBlock_adFree(1);
parameters.PF_shift = meas.prot.sWipMemBlock_alFree(7);
parameters.N_slice = size(meas.data,10);

if parameters.SMSdata==1
    parameters.SMS_shift = meas.prot.sWipMemBlock_adFree(3);
else
    parameters.SMS_shift = 1;
end
if meas.prot.pf_echo > 0
    parameters.nechoSE=size(meas.data,2);
    parameters.nechoGE=floor(meas.prot.pf_echo*parameters.nechoSE);
else
    parameters.nechoGE=size(meas.data,2);
end
parameters.N_dyn=(meas.prot.lRepetitions+1)/parameters.Nseg;

meas.data =[];
meas.data_phascor1d = []; 
meas.prot = parameters;
