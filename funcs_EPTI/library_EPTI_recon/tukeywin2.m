function f=tukeywin2(N,k1,w)
% generate 2D tukey window

    [x,y]=meshgrid(1:N(1),1:N(2));
    x_center=round((N(1)-1)/2)+1;
    y_center=round((N(2)-1)/2)+1;
    k=sqrt((x-x_center).^2+(y-y_center).^2);
%     figure;imshow(k,[]); figure;
    f=k.*0;
    ratt=(N(1)/N(2));
    kc=k1.*sqrt(((x-x_center).^2+(y-y_center).^2)./((x-x_center).^2+(y-y_center).^2*ratt));
%     figure;imshow(kc,[]); figure;
    mask=(k<=kc);
    mask(k==0)=1;
%     figure;imshow(mask,[]); figure;
    f(mask)=1;
    mask=0;
    mask=(k>kc ) &(k<kc+w);
    ff=(cos(pi*(k-kc)/2/w).^2);

    f(mask)=ff(mask);
%     figure;imshow(f,[]); 

    f(k>kc+w)=0;
end