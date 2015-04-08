% Kyle Bashford - Computational Engineering - Spring 2015

% Tranformation Rule - Integrand

% Define the integrand

f = inline('x.^(-1/2)');
a = 0;
b = 1;

% Reference value of integral

IO = quadl(f,a,b);

% Direct numerical integration with midpoint rule

dx = 0.3;
x = (a+dx/2):dx:(b-dx/2);
I1 = sum(f(x))*dx;

% Numerical integration using the same midpoint rule
% but for the transformed variable

gamma = 1/2;
f2 = inline('t.^(gamma/(1-gamma)).*((t.^(1/(1/gamma))).^(-gamma))/(1-gamma)');
N = length(x);
t1 = 0;
t2 = (b-a)^(1-gamma);
dt = (t2-t1)/N; 
t = (t1+dt/2):dt:(t2-dx/2);
I2 = sum(f2(gamma,t))*dt;