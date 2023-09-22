function b0_avg = DriftB0_esti_nav(nav1,nav2,dt)
% Estimate the global B0 offset / field drift using 1D phasecor navigators

nav1 = ifftc(nav1,2);
nav2 = ifftc(nav2,2);

temp = nav2.*conj( nav1 ) ./( abs(nav2).*abs(nav1) + eps);
temp = temp(3,:,:,:,:,:).*conj(temp(1,:,:,:,:,:)).* abs(nav1(1,:,:,:,:,:));
b0_avg = angle(squeeze(mean(mean(mean(temp,2),4),3)))/2/dt/2/pi;


