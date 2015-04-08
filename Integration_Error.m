% Kyle Bashford - Computational Engineering - Spring 2015

% Midpoint rule integration error

% Define integrand

f = inline('exp(j*k*x)');

% Define region of integration
 
a = 0;
b = 1;

% Number of integration points and step size

N = 10;
h = (b-a)/N;

% Integration Points

x = (a+h/2):h:(b-h/2);

list = 0:2:40;
for n = 1:length(list)
    k = list(n);
    
    % Reference value of integral
    
    IO = quadl(@(x) f(x,k),a,b);
    
    % Approximate value of integral using midpoint rule
    
    I = sum(f(k,x))*h;
    
    % Compute relative error
    
    err(n) = abs(I - IO)/(abs(IO));
    
    err_theory(n) = abs(1-1./sinc(k*h/(2*pi)));
    
end % loop over n1

figure(1)
plot(list,err*100,'.')
hold on
plot(list,err_theory*100,'-')
hold off
xlabel('k')
ylabel('Relative Error ( % )')
    