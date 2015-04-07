function [I_kl]=KL_grid(filename)
% function to compute indices k-l to be used in conv_grid.m
% before running run comp_alphas.m and save as temp10.mat or copy paste
% from there the parts before the for loop.
% Also note that temp10 can be swapped to differnt case 
% Author: N. Kantas, July 2012

filename='temp10.mat';
% init alpha
filevariables = { 'X' ,'Y' ,'active_grid' ,'size_active','indx'};
load(filename,filevariables{:});

I_kl=zeros(size_active,size_active);

% brute force computing of I_kl
for i1=1:size_active
    
    ii=active_grid(i1);    
    k=[X(ii);Y(ii)];
    
    for i2=1:size_active
        
        ll=active_grid(i2);
        l=[X(ll);Y(ll)];
        
        vec=k-l;
        o1 =find(X==vec(1));
        o2 =find(Y==vec(2));        
        i_special=intersect(o1,o2);
        

        if ~isempty(i_special) && ismember(i_special,active_grid) && (vec(1)^2+vec(2)^2)>0
            I_kl(i1,i2)=i_special;
        end                    
        
    end
    
end
save grid_10 I_kl
