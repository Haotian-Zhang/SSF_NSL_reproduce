%% Fig. 4

N = 30;
q = 2.0;
dt = 0.09;
M = 400;
J = 1;
L = 16;

h = L/N;                       % Space step
n = [-N/2:1:N/2-1]';           % Indices
x = n*h;                       % Grid points

% u = eval(u0); 
u = [0.5*(1+0.1*(1+x(1:N/2)/8)); 0.5*(1+0.1*(1-x(N/2+1:N)/8))];
U = u;           % Compute initial condition; save it in U
e = -4*n.*n*pi*pi/L/L;         % Squares of wavenumbers.
Uf = fftshift(fft(u));
for m = 1:1:M                  % Start time evolution

u = exp(dt*i*q*(abs(u).*abs(u))).*u;  % Solve nonlinear part of NLS
c = fftshift(fft(u));                 % Take Fourier transform
c = exp(dt*i*e).*c;                   % Advance in Fourier space
u = ifft(fftshift(c));                % Return to physical space 

if rem(m,J) == 0                      % Save output every J steps.
U = [U u];    
Uf = [Uf c];
% disp([m M])
end
                                      % Solution is stored in U.
end                                   % Display with mesh(abs(U))

figure(1)
[X,Y]=meshgrid(x,(0:M)*dt);
h = surf(X',Y',abs(U))
set(h,'LineStyle','none')
colorbar
xlabel('x')
ylabel('t')
zlabel('$|U_j^m|$','Interpreter','latex')
% mesh(abs(U))

%%
figure(2)
[X,Y]=meshgrid(0:N/2-1,(0:M)*dt);
h = surf(X',Y', abs(Uf(16:N,:)))
set(h,'LineStyle','none')
colorbar
% xlim([16,30])
xlabel('n')
% ylim([0,36])
ylabel('t')
zlabel('$|\hat{U}_n^m|$','Interpreter','latex')