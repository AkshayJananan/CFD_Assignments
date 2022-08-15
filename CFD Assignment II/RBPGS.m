%% AKSHAY J_21105012
function [phi] = RBPGS(nxp2,nyp2,dom_size,s,phi_old)
%% Intialsing
delta_x=dom_size/(nxp2);
delta_y=dom_size/(nyp2);
phi(1:nxp2,1:nyp2)=phi_old;
        pc=2*((1/(delta_x^2))+(1/(delta_y^2)));
        pcx=1/(delta_x^2);
        pcy=1/(delta_y^2);
 %% RED BLACK POINTS GAUSS SEIDEL METHOD 
    for j=2:nyp2-1
        u=mod(j,2)+2;
        for i=u:2:nxp2-1
            phi(j,i)=(-s(j-1,i-1)+phi_old(j,i-1)*pcx+phi_old(j-1,i)*pcy+phi_old(j+1,i)*pcy+phi_old(j,i+1)*pcx)/pc;
        end
    end
    phi_old=phi;
     for j=2:nyp2-1
         u=3-mod(j,2);
         for i=u:2:nxp2-1
             phi(j,i)=(-s(j-1,i-1)+phi(j-1,i)*pcy+phi(j,i-1)*pcx+phi(j+1,i)*pcy+phi(j,i+1)*pcx)/pc;
         end
     end
     
end

