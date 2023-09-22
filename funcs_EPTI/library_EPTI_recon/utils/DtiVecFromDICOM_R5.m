% DtiVecFromDICOM_R5.m
% Extract DTI gradient vectors from DICOM images
% by Xiaodong Ma, 03/25/2015
function [DTiVec,DiffusionBValue] = DtiVecFromDICOM_R5
%%
[fnames, fpathin] = uigetfile(...
                    {'All files (*.*) (Ascending order!!)';},...
                    'Open', 'MultiSelect', 'on'); 
if ~iscell(fnames) 
    fnames={fnames};
end

Ndir = input('How many diffusion directions (b0 excluded)?     ');
Check1 = input('Do iso DWI images exist?  1 yes   2 no     ');
Check2 = input('Do echo2 images exist?  1 yes   2 no     ');

addpath(fpathin);
% for i = 1:numel(fnames)
if Check1==1
    N_1slice = (Ndir+2)*2;
else
    N_1slice = (Ndir+1)*2;
end

idx_BValue = 1;
DiffusionBValue(idx_BValue)=0;

for i = 1:N_1slice
    % calculate echo number
    if Check2==1
        necho = mod(i,2);
        if necho == 1
            % calculate direction number
            ndir = (i+1)/2-1;
            if ndir>0 && ndir<Ndir+1
                im_info = dicominfo(fnames{i});
                DTiVec(:,ndir) =im_info.DiffusionGradientOrientation;
                if im_info.DiffusionBValue~=DiffusionBValue(idx_BValue);
                    idx_BValue = idx_BValue+1;
                    DiffusionBValue(idx_BValue) = im_info.DiffusionBValue;
                end
            end
        end
    else
        % calculate direction number
        ndir = i-1;
        if ndir>0 && ndir<Ndir+1
            im_info = dicominfo(fnames{i});
            DTiVec(:,ndir) =im_info.DiffusionGradientOrientation;
        end
    end
end

% print to .txt file
fid = fopen('Dti_table.txt','wt');
str_b0 = ['0: 0,0,0'];
fprintf(fid,'%s\n',str_b0);
for i=1:size(DTiVec,2)
    fprintf(fid,'%s',[num2str(i),': ']);
    fprintf(fid,'%f',DTiVec(1,i));
    fprintf(fid,'%s',',');
    fprintf(fid,'%f',DTiVec(2,i));
    fprintf(fid,'%s',',');
    fprintf(fid,'%f\n',DTiVec(3,i));
end
fclose(fid);
rmpath(fpathin);
