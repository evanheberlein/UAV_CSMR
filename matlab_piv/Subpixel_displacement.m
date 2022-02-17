function [Xdis,Ydis] = Pixel_displacement(a)
%Subpixel_displacement - Determine the subpixel part of the displacement
%
%[Xdis,Ydis] = Pixel_displacement(a)
%a - the vector of amplitudes found in the subwindow:
%    a(1) min, a(2)=a_x+1, a(3)=a_x-1, a(4)=a_y+1, a(5)=a_y-1,
%Xdis, Ydis - the x and y displacements subpixel parts found

% For MQD
% a=-a; %MQD is a local min so now it is a local max
if (a(2)<a(1) & a(3)<a(1)& a(2)>0 & a(3)>0 )
   A=log(a(3)/a(1));
   C=log(a(1)/a(2));
   Xdis=(C+A)/(2*(A-C));
else
   Xdis=2000.
end

if (a(4)<a(1) & a(5)<a(1) & a(4)>0 & a(5)>0 )
   A=log(a(5)/a(1));
   C=log(a(1)/a(4));
   Ydis=(C+A)/(2*(A-C));
else
   Ydis=2000.
end