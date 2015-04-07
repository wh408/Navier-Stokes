function [U,T]=rk4_systems(a, b, N)
%function rk4_systems() approximates the solutions of systems of m
%differential equations that are written in the form
%dy1/dt = f1(y1,y2,...,ym)
%dy2/dt = f2(y1,y2,...,ym)
%.
%.
%.
%dym/dt = fm(y1,y2,...,ym)
%with t in the interval [a; b] and the initial conditions are in the
%m-dimensional vector initial
%as with function runge_kutta4(), the inputs are the endpoints a and b, the
%number of subdivisions N in the interval [a; b], and the initial
%conditions - but this time, the initial condition is a vector
% Author: N. Kantas, July 2012


load temp10 indx X Y L
global dim;
dim=indx^2;

dx=2*pi/L;
dy=2*pi/L;
% initial=exp(-((X*dx-pi).^2+(Y*dy-pi+pi/4).^2)/(0.2))...
%     +exp(-((X*dx-pi).^2+(Y*dy-pi-pi/4).^2)/(0.2))...
%     -0.5*exp(-((X*dx-pi-pi/4).^2+(Y*dy-pi-pi/4).^2)/(0.4));

border=.1;
init=(exp(-i*2*pi/L*X*border)-ones(size(X))).*(exp(-i*2*pi/L*Y*border)-ones(size(Y)))/(4*pi*pi);
init=init./X./Y;

init(isnan(init))=0;

% m = size(initial,1);
% if m == 1
%    initial = initial';
% elseif m~=dim
%     disp('dimension mismatch in initial condition')
% end

h = (b-a)/N;        %the step size
t = a:h:b;
w(:,1) = reshape(init,dim,1);     %initial conditions


disp('Solving system of odes');
% Runge kutta order 4
for ii = 1:N
   k1 = h*f( w(:,ii));
   k2 = h*f( w(:,ii)+0.5*k1);
   k3 = h*f( w(:,ii)+0.5*k2); 
   k4 = h*f( w(:,ii)+k3);
   w(:,ii+1) = w(:,ii) + (k1 + 2*k2 + 2*k3 + k4)/6;  
end


T=t';
U=w';

subplot(121)
plot(T,abs(U));
subplot(122)
plot(T,angle(U));
ylim([-pi,pi])



%function relating the right-hand side of the differential equation
function dy = f(y)
global dim;
load temp10 alpha X Y
load grid_10 I_kl 
alpha13=alpha;
clear alpha
Lambda_k=X.^2+Y.^2;
Lambda_k=reshape(Lambda_k,dim,1);
nu=.01;
L=2;
A=nu*((2*pi/L)^2)*diag(Lambda_k);
dy= -A*y-i*conv_grid(y,'temp10.mat',alpha13,I_kl);
   