%**********************************************************
% Kyle Bashford - Spring 2015 - Computational Engineering
% 1-D FDTD Algorithm
%**********************************************************

%**********************************************************
% Physical Constants
%**********************************************************

eps0 = 8.854e-12;      % Permittivity of free space (F/m)
mu0 = 4*pi*1e-7;       % Permeability of free space (H/m)
c0 = 1/sqrt(mu0*eps0); % Speed of light in free space (m/s)

%**********************************************************
% Problem Definition
%**********************************************************

epsr1 = 1;             % Relative permittivity of region 1
epsr2 = 2;             % Relative permittivity of region 2
% c1 = c0/sqrt(epsr1);   % Speed of light in region 1
% c2 = c0/sqrt(epsr2);   % Speed of light in region 2

mu1 = 1;
mu2 = 1;
c1 = c0/sqrt(mu1/epsr1);   % Speed of light in region 1
c2 = c0/sqrt(mu2/epsr2);   % Speed of light in region 2

%**********************************************************
% Simulation Domain
%**********************************************************

D = 1;                 % Length of simulation domain (m)
tmax = 4.5e-9;         % Simulation time (s)

%**********************************************************
% Grid Generation
%**********************************************************

r = 1;                 % FDTD grid parameter
NX = 100;              % Number of spatial grid points
dx = D/(NX-1);         % Spatial grid step size (m)
dt = r*dx/c0;          % Temporal step size (s)
NT = ceil(tmax/dt);    % Number of time steps
x = 0:dx:D;            % Spatial grid points

%**********************************************************
% Constants in source
%**********************************************************

t0 = 60*dt;            % Time at pulse center (s)
s = 10*dt;             % Pulse width (s)

%**********************************************************
% Initialization
%**********************************************************

u1 = zeros(NX,1);      % Field at past time step
u2 = zeros(NX,1);      % Field at current time step
u3 = zeros(NX,1);      % Field at future time step

%**********************************************************
% Constant in update equation
%**********************************************************

a1 = (c1*dt/dx)^2*ones(NX,1);
idx2 = find(x > D/2);
a1(idx2) = (c2*dt/dx)^2;

%**********************************************************
% Constant in absorbing boundary condition
%**********************************************************

a2 = (sqrt(a1(NX))-1)/(sqrt(a1(NX)+1));

%**********************************************************
% Loop over time
%**********************************************************

figure(1)
for n = 1:NT,

    % Update equation for interior grid points
    u3(2:NX-1) = a1(2:NX-1).*(u2(3:NX)-2*u2(2:NX-1) + ...
        u2(1:NX-2)) + 2*u2(2:NX-1) - u1(2:NX-1);
    % Source at left-hand side
    t = n*dt;
    u3(1) = exp(-(t-t0)^2/(2*s^2));
    
    % Absorbing boundary condition at right-hand side
    u3(NX) = u2(NX-1) + a2*(u3(NX-1)-u2(NX));
    
    % Plot field
    plot(x,u3)
    ylim([-2 2])
    xlabel('x (m)')
    ylabel('u')
    drawnow
    
    % Shift field variables in time
    u1 = u2;
    u2 = u3;
end

%**********************************************************
% Postprocessing
%**********************************************************

idx1 = find(x < D/2);

%**********************************************************
% Get numerical reflection coefficient
%**********************************************************

[tmp,idx] = max(abs(u3(idx1)));

R1 = u3(idx)

%**********************************************************
%Exact reflection coefficient
%**********************************************************

R2 = (sqrt(1/epsr2) - sqrt(1/epsr1))/(sqrt(1/epsr2)+ sqrt(1/epsr1))