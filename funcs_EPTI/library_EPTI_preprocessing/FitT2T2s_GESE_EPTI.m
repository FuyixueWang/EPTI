function [T2_map,T2s_map,S0_map,fitting_map_GE,fitting_map_SE,mask_adapt_T2s] = FitT2T2s_GESE_EPTI(im_epti,TEs_GE,TEs_SE,te_SE,flag_joint)

% fit T2 T2* from GESE-EPTI data
% Fuyixue Wang, MGH, 2021

    if nargin < 5
        flag_joint =0;
    end
    necho_GE = length(TEs_GE);
    necho_SE = length(TEs_SE);
    echo_select_GE = 4:necho_GE;
    echo_select_SE = 1:necho_SE;

    [nx,ny,nsl,~,~] = size(im_epti);
    im_Tave = mean(abs(im_epti),5);
    im_GE = im_Tave(:,:,:,1:necho_GE);
    im_SE = im_Tave(:,:,:,necho_GE+1:end);
    im_GE = im_GE(:,:,:,echo_select_GE);
    im_SE = im_SE(:,:,:,echo_select_SE);
    im_Tave = cat(4,im_GE,im_SE);
    TEs_GE = TEs_GE(echo_select_GE);
    TEs_SE = TEs_SE(echo_select_SE);

    T2_map = zeros(nx,ny,nsl);
    T2s_map = zeros(nx,ny,nsl);
    S0_map = zeros(nx,ny,nsl);
    S02_map = zeros(nx,ny,nsl);
    fitting_map_SE = zeros(nx,ny,nsl,length(echo_select_SE));
    fitting_map_GE = zeros(nx,ny,nsl,length(echo_select_GE));
    mask_adapt_T2s = zeros(nx,ny,nsl);

    im_sos_Tave = sos(im_Tave,4);         % generate fitting mask
    mask_final = im_sos_Tave>0.8*mean(im_sos_Tave(:));
    mask_final = imfill(mask_final,'holes');
    
    for dif_sl = 1:nsl
        if(flag_joint==1)
            [S0,S02,T2,T2s,fitting_map_GE_tmp,fitting_map_SE_tmp] = T2T2s_GE_2SE(squeeze(im_Tave(:,:,dif_sl,:)),TEs_GE(:),TEs_SE(:),te_SE,mask_final(:,:,dif_sl),0);
        else
            [S0,S02,T2,T2s,fitting_map_GE_tmp,fitting_map_SE_tmp,mask_adapt_T2s_tmp] = T2T2s_GE_2SE_independent(squeeze(im_Tave(:,:,dif_sl,:)),TEs_GE(:),TEs_SE(:),te_SE,mask_final(:,:,dif_sl),0);
            mask_adapt_T2s(:,:,dif_sl) = mask_adapt_T2s_tmp;
        end
        T2(T2<0) = 500;
        T2s(T2s<0) = 500;
        T2(T2>500) = 500;
        T2s(T2s>500) = 500;
        T2_map(:,:,dif_sl) = T2;
        T2s_map(:,:,dif_sl) = T2s;
        S0_map(:,:,dif_sl) = S0;
        S02_map(:,:,dif_sl) = S02;
        fitting_map_SE(:,:,dif_sl,:) = fitting_map_SE_tmp;
        fitting_map_GE(:,:,dif_sl,:) = fitting_map_GE_tmp;
    end

end

