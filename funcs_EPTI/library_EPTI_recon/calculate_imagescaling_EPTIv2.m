function scaling = calculate_imagescaling_EPTIv2(kdata)
% Calculation of the image scaling factor for EPTI image reconstruction
% (used for matlab or BART reconstruction)
% Zijing Dong, MGH, 2020

if length(size(kdata)) > 5
    kdata = permute(kdata,[2 3 1 4 5 6]);
else
    kdata = permute(kdata,[2 3 1 4 6 5]);
end

tmp = dimnorm((ifft2c(squeeze(sum(kdata(:,:,:,:,:,:),4)))), 3);
tmpnorm = squeeze(dimnorm(tmp, 4));
tmpnorm = mean(tmpnorm,3);
tmpnorm2 = sort(tmpnorm(:), 'ascend');
p100 = tmpnorm2(end);
p90 = tmpnorm2(round(.9 * length(tmpnorm2)));
p50 = tmpnorm2(round(.5 * length(tmpnorm2)));

if (p100 - p90) < 2 * (p90 - p50)
    scaling = p90;
else
    scaling = p100;
end
% fprintf('\nScaling: %f\n\n', scaling);

end