% Basic plotting routine
% before running load u_k and T produced by ode_navier_fourier.m or
% rk4_system.m Also note that temp10 can be swapped to differnt case 
% Author: N. Kantas, July 2012

record=0;
comp_new=1;
%     close all
if comp_new==1
    tic
    load temp10 X Y indx k_max L

    magnitude_mesh=inline('sqrt(x1.^2+x2.^2)','x1','x2');

     [u_k,T]=ode_navier_fourier;
    % load u_temp_10 u_k T

    X_vec=reshape(X,indx^2,1);
    Y_vec=reshape(Y,indx^2,1);
    Z_vec=reshape(magnitude_mesh(X,Y),indx^2,1);

    Total_time=size(T,1);

    grid_L=.05;
    grid_plot=0:grid_L:L;
    size_grid=size(grid_plot,2);
    U_1=zeros(size_grid,size_grid,Total_time);
    U_2=U_1;

    for i1=1:size_grid
        for i2=1:size_grid
            i_x=grid_plot(i1);
            i_y=grid_plot(i2);

            Vec_exp=exp(i*2*pi/L*(i_x*X_vec'+i_y*Y_vec'));
            Matrx_exp=ones(Total_time,1)*Vec_exp;
            U_dot=u_k.*Matrx_exp./(ones(Total_time,1)*Z_vec');
            U_dot_1=U_dot.*(ones(Total_time,1)*Y_vec');
            U_dot_2=U_dot.*(-ones(Total_time,1)*X_vec');
            U_1(i1,i2,:)=nansum(U_dot_1,2);
            U_2(i1,i2,:)=nansum(U_dot_2,2);
        end
    end

    toc

    [Ix,Iy]=meshgrid(grid_plot,grid_plot);

end

% comment out stuff from now on to produce and avi
% doesnt work in mac for some reason

fig=figure(33);
 set(fig,'DoubleBuffer','on');
 set(gca,'xlim',[-80 80],'ylim',[-80 80],...
     'NextPlot','replace','Visible','off')

dt=10;
if record == 1
     fig1=figure(13);
     winsize = get(fig1,'Position');
     winsize(1:2) = [0 0];
     A=moviein(ceil(Total_time/dt),fig1);
    mov = avifile('example.avi');
end

for t=1:dt:Total_time
    quiver(Ix,Iy,reshape(U_1(:,:,t),size_grid,size_grid),reshape(U_2(:,:,t),size_grid,size_grid))
    title(['Time is ' num2str(T(t)) ' of '  num2str(T(end))])
    xlim([-0.1 2.1])
    ylim([-0.1 2.1])
    xlabel('u-direction')
    ylabel('v-direction')
         Amov(:,t)=getframe(fig1,winsize);
    if record == 1
        frame = getframe(gca);
        mov = addframe(mov,frame);
    end    
    pause
end

if record == 1
    mov = close(mov);

     movie(Amov)
     mpgwrite(A,jet,'movie1.mpg');
end