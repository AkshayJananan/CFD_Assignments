%% AKSHAY J_21105012
function [phi] = gauss_seidel(nxp2,nyp2,dom_size,s,phi_old,phi_exact)
%% Initialising
delta_x=dom_size/(nxp2);
delta_y=dom_size/(nyp2);
phi(1:nxp2,1:nyp2)=phi_old;
x_array=linspace(delta_x/2,1-delta_x/2,nxp2-2);
y_array=linspace(delta_y/2,1-delta_y/2,nyp2-2);
error_tol=10^(-6);%tolerance level
        pc=2*((1/(delta_x^2))+(1/(delta_y^2)));
        pcx=1/(delta_x^2);
        pcy=1/(delta_y^2);
%% GAUSS-SEIDEL METHOD
    for j=2:nyp2-1
        for i=2:nxp2-1
            phi(j,i)=(-s(j-1,i-1)+phi(j-1,i)*pcy+phi(j,i-1)*pcx+phi_old(j+1,i)*pcy+phi_old(j,i+1)*pcx)/pc;
        end
    end
end

