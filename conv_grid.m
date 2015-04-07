function [Out]=conv_grid(U,filename,alpha13,I_kl)
% function to compute convolution sum for the gradient in the system of
% odes
% It is used within ode_navier_fourier or rk4_systems.
% NOTE this computation is very naive and can be speeded up using FFT etc
% Author: N. Kantas, July 2012

dim=length(U);

filevariables = { 'X' ,'Y' ,'active_grid' ,'size_active','indx'};
load(filename,filevariables{:});

if indx^2~=dim
    disp('lengths of vectors not same as file')
end

Out=zeros(dim,1);

for i1=1:size_active
    
    ii=active_grid(i1);    
    k=[X(ii);Y(ii)];
    temp=0;
    
    for i2=1:size_active
        
        ll=active_grid(i2);
        l=[X(ll);Y(ll)];
        
        i_special=I_kl(i1,i2); 
        if i_special>0
            temp=temp+alpha13(i1,i2)*U(i2)*U(i_special); 
        
        end                    
        
    end
    Out(i1,1)=temp;
end
