function saveniftidata_EPTI(data,meas,fndata)
% save NIFTI data with correct header info read from rawdata
% data size [nx,ny,nsl,ndyns]

% by Fuyixue Wang based on Jon's mrir
% v2, output NIFTI without mri convert from MGZ

    disp(sprintf('==> [%s]: saving file "%s"...', mfilename, fndata));
    meas.prot.tPatientPosition = 'HFS';

    M0_vox2ras = mrir_measdat_vox2ras_EPTI(meas.prot, meas.evp);

    res = [meas.prot.dReadoutFOV/meas.prot.lBaseResolution, meas.prot.dPhaseFOV/(meas.prot.Nseg*meas.prot.Rseg + meas.prot.PF_shift), meas.prot.sSliceArray(1).dThickness];

    save_nii(make_nii_EPTI(double(data), [res,meas.prot.alTR*meas.prot.Nseg/1000], [], 16,M0_vox2ras), fndata);

end

