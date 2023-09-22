function [ FA ] = Readimg( Name )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
info=analyze75info([Name,'.hdr']);
FA=double(analyze75read([Name,'.img']));
FA=FA(end:-1:1,end:-1:1,:);
figure; imshow(FA,[]);

end

