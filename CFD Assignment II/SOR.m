%% AKSHAY J_21105012
function [phi] = SOR(nxp2,nyp2,dom_size,s,phi_old)
%% Intialising
delta_x=dom_size/(nxp2);
delta_y=dom_size/(nyp2);
v=1.2;
phi(1:nxp2,1:nyp2)=phi_old;
for j=2:nyp2-1
    for i=1:nxp2-1
        pc(j,i)=2*((1/(delta_x^2))+(1/(delta_y^2)));
        pcx(j,i)=1/(delta_x^2);
        pcy(j,i)=1/(delta_y^2);
    end
end
%% SUCCESSIVE OVER RELAXATION METHOD
    for j=2:nyp2-1
        for i=2:nxp2-1
             phi(j,i)=(v*(-s(j-1,i-1)+phi(j,i-1)*pcx(j,i)+phi(j-1,i)*pcy(j,i)+phi_old(j+1,i)*pcy(j,i)+phi_old(j,i+1)*pcx(j,i))/pc(j,i))+(1-v)*phi_old(j,i);
        end
    end
end

