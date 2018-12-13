%% Amath 581 Final Project
% Author: Haotian Zhang
% Date: 11/15/2018

%% Figure 1
N = 50;
q = 2.0;
dt = 0.005;
M = 4800;
J = 1;
L = 16;

h = L/N;                       % Space step
n = [-N/2:1:N/2-1]';           % Indices
x = n*h;                       % Grid points
r = dt/h^2;
theta=0.5;
e1=ones(N,1);
S=spdiags([e1 -2*e1 e1],[-1 0 1],N,N);S(1,N)=1; S(N,1)=1;
I = speye(N);

u = 0.5*(1+0.1*cos(pi*x/8));
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
h = surfc(X',Y',abs(U));
set(h,'LineStyle','none')
colorbar
xlabel('x')
ylabel('t')
zlabel('$|U_j^m|$','Interpreter','latex')
% mesh(abs(U))


%% Figure 2
figure(2)
[X,Y]=meshgrid(0:N/2-1,(0:M)*dt);
h = surf(X',Y', abs(Uf(26:N,:)));
set(h,'LineStyle','none')
colorbar
xlabel('n')
ylabel('t')
zlabel('$|\hat{U}_n^m|$','Interpreter','latex')

%% Figure 3
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

u = ones(N,1)*0.5;

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

figure(3)
[X,Y]=meshgrid(x,(0:M)*dt);
h = surf(X',Y',abs(U));
set(h,'LineStyle','none')
colorbar

xlabel('x')
ylabel('t')
zlabel('$|U_j^m|$','Interpreter','latex')
mesh(abs(U))

%% Figure 4
figure(4)
[X,Y]=meshgrid(0:N/2-1,(0:M)*dt);
h = surf(X',Y', abs(Uf(11:N,:)))
set(h,'LineStyle','none')
colorbar

xlabel('n')
ylabel('t')
zlabel('$|\hat{U}_n^m|$','Interpreter','latex')

%% 

