%% AKSHAY J_21105012
function [phi] = jacobi(nxp2,nyp2,dom_size,s,phi_old)
%% initialising 
delta_x=dom_size/(nxp2);
delta_y=dom_size/(nyp2);
phi(1:nxp2,1:nyp2)=phi_old;
for j=2:nyp2-1
    for i=1:nxp2
        pc(j,i)=2*((1/(delta_x^2))+(1/(delta_y^2)));
        pcx(j,i)=1/(delta_x^2);
        pcy(j,i)=1/(delta_y^2);
    end
end 
%%JACOBI METHOD
    for j=2:nyp2-1
        for i=2:nxp2-1
            phi(j,i)=(-s(j-1,i-1)+phi_old(j,i-1)*pcx(j,i)+phi_old(j-1,i)*pcy(j,i)+phi_old(j+1,i)*pcy(j,i)+phi_old(j,i+1)*pcx(j,i))/pc(j,i);
        end
    end
end
   