function [Y,T]=ode_navier_fourier
% Navier stokes spectral Galerkin routine
% ode obtained
% before running run comp_alphas.m and save as temp10.mat and the KL_grid.m
% which saves as grid_10
% Also note that temp10 can be swapped to differnt case 
% Author: N. Kantas, July 2012

tic
load temp10 indx X Y L
global dim;
dim=indx^2;
Tfinal=50;

dx=2*pi/L;
dy=2*pi/L;

border1=.5;
border2=.7;


border1y=.1;
border2y=1.1;

% different initial conditions

% initial=exp(-((X*dx-pi).^2+(Y*dy-pi+pi/4).^2)/(0.2))...
%     +exp(-((X*dx-pi).^2+(Y*dy-pi-pi/4).^2)/(0.2))...
%     -0.5*exp(-((X*dx-pi-pi/4).^2+(Y*dy-pi-pi/4).^2)/(0.4));

% init=(exp(-i*2*pi/L*X*border)-ones(size(X))).*(exp(-i*2*pi/L*Y*border)-ones(size(Y)))/(4*pi*pi);

init=(exp(-i*2*pi/L*X*border2)-exp(-i*2*pi/L*X*border1)).*(exp(-i*2*pi/L*Y*border2y)-exp(-i*2*pi/L*Y*border1y))/(4*pi*pi);
init=init./X./Y;

init(isnan(init))=0;

initial = reshape(init,dim,1);     %initial conditions
disp('Solving system of odes');

options=odeset('MaxStep',0.001*Tfinal);
% options=odeset('RelTol',1e-6);

[T,Y] = ode23(@yprime1, [0, Tfinal],initial,options);
toc
subplot(121)
plot(T,abs(Y));
subplot(122)
plot(T,angle(Y));
% ylim([-pi,pi])


function dy = yprime1(t,y)
global dim;
load temp10 alpha X Y
load grid_10 I_kl 
alpha13=alpha;
Lambda_k=X.^2+Y.^2;
Lambda_k=reshape(Lambda_k,dim,1);
% dy=ones(dim,1);
nu=.01;
L=2;
A=nu*((2*pi/L)^2)*diag(Lambda_k);
dy= -A*y-i*conv_grid(y,'temp10.mat',alpha13,I_kl);


