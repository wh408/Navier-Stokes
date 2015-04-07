% function to compute alpha_{k,l} constants to be used in conv_grid.m
% script saves as temp10.mat so if this changes change accordingly rest of
% scripts
% Author: N. Kantas, July 2012

clear
tic 

magnitude=inline('sqrt(x(1)^2+x(2)^2)','x');

magnitude_mesh=inline('sqrt(x1.^2+x2.^2)','x1','x2');

lambda_k=inline('((2*pi/L)^2)*(x(1)^2+x(2)^2)','x','L');

lambda_k_mesh=inline('((2*pi/L)^2)*(x1.^2+x2.^2)','x1','x2','L');


dot_kl=inline('abs(x(1)^2+x(2)^2 - x(1)*y(1)-x(2)*y(2))','x','y'); % x, y 2-D vectors

cross_kl=inline('abs(x(1)*y(2)-x(2)*y(1))','x','y'); % x, y 2-D vectors

% torus size
L=2;

% discretisation of k=[k_x;k_y] with max k_x = max k_y = k_max  
k_max=10;

lambda2=k_max^2*4*(pi^2)/(L^2);

% X and Y components of k on grid
[X,Y]=meshgrid(-k_max:1:k_max,-k_max:1:k_max);

% indicator whether point on grid is in circle
Z=[magnitude_mesh(X,Y)<=k_max];

% number of elements in X,Y,Z
indx=2*k_max+1;

% vector of nonzero elements (in Z)
active_grid=find(Z>0)';
size_active=size(active_grid,2);

% init alpha
alpha=zeros(size_active,size_active);

% brute force computing of alpha
for i1=1:size_active
    ii=active_grid(i1);
    for i2=1:size_active
        ll=active_grid(i2);
        
        
        k=[X(ii);Y(ii)];
        l=[X(ll);Y(ll)]; 
        
        m1=magnitude(k-l);
        m2=magnitude(l);
        m3=magnitude(k);
        
        if m3>0 && m2>0 && m1>0
        	alpha(i1,i2)=dot_kl(k,l)*cross_kl(k,l)/(m3*m2*m1);      
        
        	% if Z(ii)==0 || Z(ll)==0
            %         disp(['....' num2str([Z(ii),Z(ll)])])
            %         pause(.1)
            % disp(['..neeed check..'])
            
        end
    end    
end

alpha=2*pi*alpha/L;

simtime=toc;

ZZZ=magnitude_mesh(X(active_grid),Y(active_grid));

figure
mesh(ZZZ,ZZZ,alpha)
xlabel('|k|')
ylabel('|l|')
zlabel('Im(\alpha_{k}^{l,k-l})')

% saveas(gcf,'alphas10.jpg')
% 
save temp10









