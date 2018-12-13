%% Fig. 1
clear;clc;
N = 20;
q = -4.0;
dt = 0.036;
M = 5000;
J = 1;
L = 1;

h = L/N;                       % Space step
n = [-N/2:1:N/2-1]';           % Indices
x = n*h;                       % Grid points
r = dt/h^2;
theta= 0.50;
e1=ones(N,1);
S=spdiags([e1 -2*e1 e1],[-1 0 1],N,N);S(1,N)=1; S(N,1)=1;
I = speye(N);
% u = eval(u0); 
% u = [0.5*(1+0.1*(1+x(1:N/2)/8)); 0.5*(1+0.1*(1-x(N/2+1:N-1)/8))];
u = ones(N,1)*0.5;
% u = 0.5*(1+0.1*cos(pi*x/8));
U = u;           % Compute initial condition; save it in U
Uf = fftshift(fft(u));
for m = 1:1:M                  % Start time evolution

    vm = exp(1i*dt*q*(abs(u).*abs(u))).*u;
    
    u = (I-1i*r*theta*S)^-1*(I+1i*r*(1-theta)*S)*vm;
    c = fftshift(fft(u));                 % Take Fourier transform
    if rem(m,1) == 0
        U = [U u];    
        Uf = [Uf c];
    end
end
                                      % Solution is stored in U.

figure(1)
[X,Y]=meshgrid(x,(0:M)*dt);
mesh(X',Y',abs(U))
xlabel('x')
ylabel('t')
zlabel('$|U_j^m|$','Interpreter','latex')
mesh(abs(U))
figure(2)
[X,Y]=meshgrid(0:N/2-1,(0:M)*dt);
mesh(X',Y', abs(Uf(11:N,:)))
% xlim([16,30])
xlabel('n')
% ylim([0,36])
ylabel('t')
zlabel('$|\hat{U}_n^m|$','Interpreter','latex')