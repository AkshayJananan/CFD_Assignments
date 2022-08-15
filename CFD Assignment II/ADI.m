%% AKSHAY J_21105012
function [phi] = ADI(nxp2,nyp2,dom_size,s,phi_old)
%% Intialising
delta_x=dom_size/(nxp2);
delta_y=dom_size/(nyp2);
x_array=linspace(delta_x/2,1-delta_x/2,nxp2-2);
y_array=linspace(delta_y/2,1-delta_y/2,nyp2-2);
pc=2*((1/delta_x)^2+(1/delta_y)^2);
pcx=(1/delta_x)^2;
pcy=(1/delta_y)^2;
phi(1:nxp2,1:nyp2)=phi_old;
%% ALTERNATING DIRECTION IMPLICIT
%% ROW SWEEP
       a(1:nxp2-2)=-pcx;
       b(1:nxp2-2)=pc;
       c(1:nxp2-2)=-pcx;
       d(1:nxp2-2)=0;
        for j=2:nyp2-1
          for i=2:nxp2-1
              d(i)=pcy*phi(j-1,i)+-pcy*phi_old(j+1,i)-s(j-1,i-1);
          end
           phi_old(j,2:nxp2-1)=TDMA(a,b,c,d,nxp2-2);
        end
        phi_old=phi;
        %UPDATING BOUNDARY CONDITION
        for i=2:nxp2-1
          phi_old(1,i)=2*(0.25*sinh(-5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+1)-phi_old(2,i);
          phi_old(nyp2,i)=2*(0.25*sinh(5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+exp(2*x_array(i-1)))-phi_old(nyp2-1,i);
        end
       for j=2:nyp2-1
         phi_old(j,1)=2*(0.25*sinh(-5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+1)-phi_old(j,2);
         phi_old(j,nxp2)=2*(0.25*sinh(5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+exp(2*y_array(j-1)))-phi_old(j,nxp2-1);
       end
%% COLUMN SWEEP
       a(1:nyp2-2)=-pcy;
       b(1:nyp2-2)=pc;
       c(1:nyp2-2)=-pcy;
       d(1:nyp2-2)=0;
        for i=2:nyp2-1
          for j=2:nxp2-1
              d(j)=pcx*phi(j-1,i)+pcx*phi_old(j+1,i)-s(j-1,i-1);
          end
           phi(2:nxp2-1,i)=TDMA(a,b,c,d,nyp2-2);
        end
      %UPDATING BOUNDARY CONDITION
        for i=2:nxp2-1
         phi_old(1,i)=2*(0.25*sinh(-5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+1)-phi_old(2,i);
         phi_old(nyp2,i)=2*(0.25*sinh(5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+exp(2*x_array(i-1)))-phi_old(nyp2-1,i);
        end
       for j=2:nyp2-1
         phi_old(j,1)=2*(0.25*sinh(-5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+1)-phi_old(j,2);
         phi_old(j,nxp2)=2*(0.25*sinh(5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+exp(2*y_array(j-1)))-phi_old(j,nxp2-1);
end

