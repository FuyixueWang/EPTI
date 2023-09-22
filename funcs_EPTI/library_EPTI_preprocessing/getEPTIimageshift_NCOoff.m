function shift_y = getEPTIimageshift_NCOoff(prot)
% calculate the image shift due to NCO turned off in EPTI sequence
% Fuyixue Wang, 2023

prot.tPatientPosition = 'HFS';

Q1_LPS = [prot.sSliceArray(1).sPosition_dSag; ...
          prot.sSliceArray(1).sPosition_dCor; ...
          prot.sSliceArray(1).sPosition_dTra];

Q2_LPS = [prot.sSliceArray(2).sPosition_dSag; ...
          prot.sSliceArray(2).sPosition_dCor; ...
          prot.sSliceArray(2).sPosition_dTra];
Q_delta_LPS = Q2_LPS - Q1_LPS;

Q_center_LPS = Q1_LPS + Q_delta_LPS*(prot.sSliceArray_lSize-1)/2;

R_LPS2RAS = [-1 0 0; 0 -1 0; 0 0 +1];
Q_center_RAS = R_LPS2RAS * Q_center_LPS;

[COL_dircos_LPS, LIN_dircos_LPS, IMG_dircos_LPS] = mrir_measdat_vox2ras__dircos(prot);
COL_dircos_RAS = R_LPS2RAS * COL_dircos_LPS;
LIN_dircos_RAS = R_LPS2RAS * LIN_dircos_LPS;
IMG_dircos_RAS = R_LPS2RAS * IMG_dircos_LPS;

% temporary hack for debugging
COL_dircos_RAS = -COL_dircos_RAS; % FW uncomment
LIN_dircos_RAS = -LIN_dircos_RAS; % FW uncomment
%IMG_dircos_RAS = -IMG_dircos_RAS;

R = [COL_dircos_RAS, LIN_dircos_RAS, IMG_dircos_RAS];

Q_center_vox = R'*Q_center_RAS;

shift_y = -Q_center_vox(2);

end

