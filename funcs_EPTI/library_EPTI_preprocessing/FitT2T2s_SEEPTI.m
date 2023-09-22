function [T2_map,T2s_map,S0_map,fitting_map_SE] = FitT2T2s_SEEPTI(im_epti,TEs,te_se)

% fit T2 T2* from SE-EPTI data
% Zijing Dong, MGH, 2021

    im_Tave = mean(abs(im_epti),5);
%     im_Tave = abs(mean((im_epti),5));

    [nx,ny,nsl,necho0,nDyn] = size(im_epti); 
    echo_select = 1:necho0;
    im_Tave = im_Tave(:,:,:,echo_select);
    TEs = TEs(echo_select);

    T2_map = zeros(nx,ny,nsl);
    T2s_map = zeros(nx,ny,nsl);
    S0_map = zeros(nx,ny,nsl);
    fitting_map_SE = zeros(nx,ny,nsl,length(echo_select));

    im_sos_Tave = mean(squeeze(sos(im_epti,4)),4);         % generate fitting mask
    mask_final = im_sos_Tave>0.8*mean(im_sos_Tave(:));
    mask_final = imfill(mask_final,'holes');
    
    for dif_sl = 1:nsl
        [S0, T2, T2s,fitting_map] = T2T2s_SE_EPTI(squeeze(im_Tave(:,:,dif_sl,:)),TEs(:),te_se,mask_final(:,:,dif_sl),0);
        T2(T2<0) = 500;
        T2s(T2s<0) = 500;
        T2(T2>500) = 500;
        T2s(T2s>500) = 500;
        T2_map(:,:,dif_sl) = T2;
        T2s_map(:,:,dif_sl) = T2s;
        S0_map(:,:,dif_sl) = S0;
        fitting_map_SE(:,:,dif_sl,:) = fitting_map;
    end

end

