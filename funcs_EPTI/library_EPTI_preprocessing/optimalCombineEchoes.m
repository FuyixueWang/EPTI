function [img_comb,T2s_map,S0_map,fitting_map,adapt_mask] = optimalCombineEchoes(im_epti,TEs_GE0)

% T2*-optimized combined
% Fuyixue Wang, MGH, 2021

    im_Tave = mean(abs(im_epti),5);

    [nx,ny,nsl,necho0,nDyn] = size(im_epti); 
    echo_select = 4:necho0;
    im_Tave = im_Tave(:,:,:,echo_select);
    TEs_GE = TEs_GE0(echo_select);

    T2s_map = zeros(nx,ny,nsl);
    S0_map = zeros(nx,ny,nsl);
    adapt_mask = zeros(nx,ny,nsl);
    fitting_map = zeros(nx,ny,nsl,length(echo_select));

    im_sos_Tave = mean(squeeze(sos(im_epti,4)),4);         % generate fitting mask
    mask_final = im_sos_Tave>0.2*mean(im_sos_Tave(:));
    mask_final = imfill(mask_final,'holes');
    
    for dif_sl = 1:nsl
        [S0,T2s,fitting_map_GE,mask_adapt] = T2s_GE_EPTI_adaptiveMask(squeeze(im_Tave(:,:,dif_sl,:)),TEs_GE(:)*1000,mask_final(:,:,dif_sl),0.1);
        T2s(T2s<0) = TEs_GE(end)*1000*4;
        T2s(T2s>TEs_GE(end)*1000*4) = TEs_GE(end)*1000*4;
        T2s_map(:,:,dif_sl) = T2s;
        S0_map(:,:,dif_sl) = S0;
        adapt_mask(:,:,dif_sl) = mask_adapt;
        fitting_map(:,:,dif_sl,:) = fitting_map_GE;
    end
    
    img_weight = zeros(nx,ny,nsl,necho0);
    for dif_echo = 1:necho0
        img_weight(:,:,:,dif_echo) = TEs_GE0(dif_echo)*1000*exp(-TEs_GE0(dif_echo)*1000./T2s_map); 
    end
    img_weight = img_weight./repmat(sum(img_weight,4),[1,1,1,necho0]);
    img_weight(isnan(img_weight)) = 0;
    img_comb = squeeze(sum(abs(im_epti).*repmat(img_weight,[1,1,1,1,nDyn]),4));
end

