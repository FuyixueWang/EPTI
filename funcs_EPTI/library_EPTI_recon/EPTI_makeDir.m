function [sub_dir,subdir1_acq,subdir2_recon,subdir3_nii] = EPTI_makeDir(directory,filename)
% make directories to save the generated data files
sub_dir = fullfile(directory,'EPTIdata');
subdir1_acq = fullfile(sub_dir,'0_Data_acq',filename);
subdir2_recon = fullfile(sub_dir,'1_Recon_mat',filename);
subdir3_nii = fullfile(sub_dir,'2_Recon_nii',filename);
subdir3_nii_echoes = fullfile(sub_dir,'2_Recon_nii',filename,'echoes');

mkdir(sub_dir);
mkdir(subdir1_acq);
mkdir(subdir2_recon);
mkdir(subdir3_nii);
mkdir(subdir3_nii_echoes);

mkdir([directory,'tmp']);

end