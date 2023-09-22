function [ img ] = ReadVec( filename )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[hdr, filetype, fileprefix, machine] = load_nii_hdr([filename,'.hdr']);
[img] = load_nii_img(hdr,filetype,fileprefix,machine);
end