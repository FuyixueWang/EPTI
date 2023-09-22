function [tSNR_img,img_std,img_mean] = tSNR(img,detrend_order)
% calculating tSNR with detrend
disp('calculating tSNR ...');
[nx,ny,nsl,nt] = size(img);
img_detrend = zeros(nx,ny,nsl,nt);
for sl = 1:nsl
    for x = 1:nx
        for y = 1:ny
            img_detrend(x,y,sl,:) = detrend(squeeze(img(x,y,sl,:)),detrend_order);
        end
    end
end
img_mean = mean(img,4);
img_std = std(img_detrend,[],4);
tSNR_img = img_mean./img_std;
tSNR_img(isnan(tSNR_img)) = 0;
end

