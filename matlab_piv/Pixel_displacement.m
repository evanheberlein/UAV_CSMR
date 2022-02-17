function [a,nXdis,nYdis] = Pixel_displacement(D,Nx,Ny)
%Pixel_displacement - Determine the integer size of the displacement
%(nXdis,nYdis) for a subwinodow of size (Nx,Ny).
%
%[min_ampl,nXdis,nYdis] = Pixel_displacement(D,Nx,Ny)
%a - the vector of amplitudes found in the subwindow:
%    a(1) min, a(2)=a_x+1, a(3)=a_x-1, a(4)=a_y+1, a(5)=a_y-1,
%nXdis, nYdis - the x and y displacements found (corrected for phase)
%D - the positive definite array to search
%Nx,Ny - the x and y size of the array

%For MQD find minimum
%[a(1),nYdis]=min(min(D,[],2));
%[a(1),nXdis]=min(min(D,[],1));

%For PIV maximum
[a(1),nYdis]=max(max(D,[],2));
[a(1),nXdis]=max(max(D,[],1));

%{
figure(12)
 imagesc(D)
 axis image
 colorbar
 sprintf('a, nXdis, nYdis: %6.2f %d %d',a,nXdis,nYdis)
pause
%}
if nXdis+1>Nx
    indXp1=1;
else
    indXp1=nXdis+1;
end
if nXdis-1<1
    indXm1=Nx;
else
    indXm1=nXdis-1;
end
if nYdis+1>Ny
    indYp1=1;
else
    indYp1=nYdis+1;
end
if nYdis-1<1
    indYm1=Ny;
else
    indYm1=nYdis-1;
end

a(2)=D(nYdis,indXp1);
a(3)=D(nYdis,indXm1);
a(4)=D(indYp1,nXdis);
a(5)=D(indYm1,nXdis);

nYdis=nYdis-1; %correct for lack of 0 lag
nXdis=nXdis-1; %correct for lack of 0 lag
      
if nYdis>=Ny/2
   nYdis=Ny-nYdis;
else
   nYdis=-nYdis;
end
      
if nXdis>=Nx/2
   nXdis=nXdis-Nx;
end