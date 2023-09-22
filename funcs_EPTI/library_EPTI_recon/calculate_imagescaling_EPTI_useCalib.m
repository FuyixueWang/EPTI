function scaling = calculate_imagescaling_EPTI_useCalib(im)
% Calculation of the image scaling factor for EPTI image reconstruction
% (used for matlab or BART reconstruction)
% Zijing Dong, MGH, 2020

tmp = sos(im, 4);
for slice = 1:size(tmp,3)
    tmpnorm = tmp(:,:,slice);
    tmpnorm2 = sort(tmpnorm(:), 'ascend');
%     p100 = tmpnorm2(end);
%     p90 = tmpnorm2(round(.9 * length(tmpnorm2)));
%     p50 = tmpnorm2(round(.5 * length(tmpnorm2)));
%     if (p100 - p90) < 2 * (p90 - p50)
%         scaling(slice) = p90;
%     else
%     scaling(slice) = p100;
%     end
    p98 = tmpnorm2(round(.98 * length(tmpnorm2)));
    scaling(slice) = p98;

end

end