% Spring 2015 - Computational Engineering - Patanjali Parimi

%***********************************************************************
% 1-D FDTD code with simple radiation boundary conditions
%***********************************************************************
%
% Implement the finite-difference time-domain
% solution of Maxwell's curl equations over a one-dimensional space
% comprised of uniform grid cells.
%
% Model a sinusoidal wave (1GHz) propagating in a lossy medium (epsr=1.0, sigma=5.0e-3 S/m)
% The simplified finite difference system for nonpermeable media is implemented.
%
% Set grid resolution (dx = 1.5 cm) to provide 20
% samples per wavelength. Set Courant factor S=c*dt/dx
% stability limit: S=1.
%
% Use the simplest radiation boundary condition for wave propagation in free space:
%
% Ez(imax,n+1) = Ez(imax-1,n)
% Change stability criterion and see what happens to the solution
%**************************************************************************
% SOLUTION
% To execute this M-file, type "fdtd1D" at the MATLAB prompt. ???
% M-file displays the FDTD-computed Ez and Hy fields at every
% time step, and records those frames in a movie matrix, M, which is
% played at the end of the simulation using the "movie" command.
%
%***********************************************************************

clear

%***********************************************************************
% Fundamental constants
%***********************************************************************

cc=2.99792458e8; %speed of light in free space
muz=4.0*pi*1.0e-7; %permeability of free space
epsz=1.0/(cc*cc*muz); %permittivity of free space
freq=1.0e+9; %frequency of source excitation
lambda=cc/freq; %wavelength of source excitation
omega=2.0*pi*freq; 

%***********************************************************************
% Grid parameters
%***********************************************************************

ie=200; %number of grid cells in x-direction
ib=ie+1;
dx=lambda/20.0; %space increment of 1-D lattice
dt=dx/cc; %time step
omegadt=omega*dt;
nmax=round(12.0e-9/dt); %total number of time steps

%***********************************************************************
% Material parameters
%***********************************************************************

eps=2.0;
sig=5.0e-3;

%***********************************************************************
% Updating coefficients for space region with nonpermeable media
%***********************************************************************

scfact=dt/muz/dx;
ca=(1.0-(dt*sig)/(2.0*epsz*eps))/(1.0+(dt*sig)/(2.0*epsz*eps));
cb=scfact*(dt/epsz/eps/dx)/(1.0+(dt*sig)/(2.0*epsz*eps));

%***********************************************************************
% Field arrays
%***********************************************************************

ez(1:ib)=0.0;
hy(1:ie)=0.0;

%***********************************************************************
% Movie initialization
%***********************************************************************

x=linspace(dx,ie*dx,ie);
subplot(2,1,1),plot(x,ez(1:ie)/scfact,'r'),axis([0 3 -1 1]);
ylabel('EZ');
subplot(2,1,2),plot(x,hy,'b'),axis([0 3 -3.0e-3 3.0e-3]);
xlabel('x (meters)');ylabel('HY'); 
rect=get(gcf,'Position');
rect(1:2)=[0 0];
M=moviein(nmax/2,gcf,rect);

%***********************************************************************
% BEGIN TIME-STEPPING LOOP
%***********************************************************************

for n=1:nmax

%***********************************************************************
% Update electric fields
%***********************************************************************

ez(1)=scfact*sin(omegadt*n);
rbc=ez(ie);
ez(2:ie)=ca*ez(2:ie)+cb*(hy(2:ie)-hy(1:ie-1));
ez(ib)=rbc;

%***********************************************************************
% Update magnetic fields
%***********************************************************************

hy(1:ie)=hy(1:ie)+ez(2:ib)-ez(1:ie);

%***********************************************************************
% Visualize fields
%***********************************************************************

if mod(n,2)==0;
rtime=num2str(round(n*dt/1.0e-9));
subplot(2,1,1),plot(x,ez(1:ie)/scfact,'r'),axis([0 3 -1 1]);
title(['time = ',rtime,' ns']);
ylabel('EZ');
subplot(2,1,2),plot(x,hy,'b'),axis([0 3 -3.0e-3 3.0e-3]);
title(['time = ',rtime,' ns']);
xlabel('x (meters)');ylabel('HY');
M(:,n/2)=getframe(gcf,rect);
end

%***********************************************************************
% END TIME-STEPPING LOOP
%***********************************************************************

end
movie(gcf,M,0,10,rect);
