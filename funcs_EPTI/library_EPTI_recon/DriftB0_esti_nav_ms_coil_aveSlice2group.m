function b0_avg = DriftB0_esti_nav_ms_coil_aveSlice2group(nav1,nav2,dt)
% estimate the B0 offset using 1D phase cor navigators
% version with shot-to-shot phase for multi-shot acquisition
% Fuyixue Wang 
% 2023, MGH

nav1 = ifftc(nav1,2);
nav2 = ifftc(nav2,2);
nav1 = mean(nav1,3);
temp = nav2.*conj( nav1 ) ./( abs(nav2).*abs(nav1) + eps);
temp = temp(3,:,:,:,:,:).*conj(temp(1,:,:,:,:,:)).* abs(nav1(1,:,:,:,:,:));
b0_avg(:,:,:,1) = angle(squeeze(mean(mean(temp(:,:,:,:,1:2:end),2),5)))/2/dt/2/pi;
b0_avg(:,:,:,2) = angle(squeeze(mean(mean(temp(:,:,:,:,2:2:end),2),5)))/2/dt/2/pi;
