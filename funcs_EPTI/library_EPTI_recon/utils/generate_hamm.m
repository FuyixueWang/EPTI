function window = generate_hamm(N,Nfull)

center = (1 + N)/2;
size = round(N*sqrt(2));
window = zeros(N,N);

for cc=1:N
    
    for dd=1:N
        
        r = sqrt((cc-center)^2 + (dd-center)^2);

        
        window(cc,dd) = 0.54+0.46*cos(2*pi*r/size);

    end
            
end

 window = (window-min(window(:)))/(max(window(:))-min(window(:)));
 window = (window-min(window(:)))/(max(window(:))-min(window(:)));
 if nargin > 1
 window = padarray(window,[round((Nfull-N)/2) round((Nfull-N)/2)]);
 end
 

end