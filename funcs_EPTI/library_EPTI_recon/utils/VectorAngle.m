function [ angle ] = VectorAngle( a,b )
%Calculate the angle between to vectors

angle=rad2deg(acos(dot(a,b)./(norm(a).*norm(b))));


end