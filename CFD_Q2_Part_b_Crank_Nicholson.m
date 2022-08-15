clear all;clc;
close all;
%%
% AKSHAY J_21105012
%%
%defining the grid
T=42;
step_size_array=[0.1,0.6,2.1];
%initial values
phi_1(1)=1;
phi_2(1)=1;
phi_3(1)=1;
u(1)=0;%u(t) is assumed to be the first derivative of phi wrt t
%%
%Crank Nicholson Calculation
N=T/1;
i=0;
 for t=0:0.1:N
    phi_exact(i+1)=cos(t);
    i=i+1;
 end
for i=step_size_array(1:end)
 N=T/i;
for t=1:N
    if i==step_size_array(1)
          u(t+1)=(u(t).*(1-0.25.*i.^2)-(i.*phi_1(t)))./(1+0.25.*i.^2);
          phi_1(t+1)=phi_1(t)+(0.5.*i.*(u(t)+u(t+1)));
    elseif i==step_size_array(2)
         u(t+1)=(u(t).*(1-0.25.*i.^2)-(i.*phi_2(t)))./(1+0.25.*i.^2);
          phi_2(t+1)=phi_2(t)+(0.5.*i.*(u(t)+u(t+1)));
    elseif i==step_size_array(3)
         u(t+1)=(u(t).*(1-0.25.*i.^2)-(i.*phi_3(t)))./(1+0.25.*i.^2);
          phi_3(t+1)=phi_3(t)+(0.5.*i.*(u(t)+u(t+1)));
      end
 
end
end
%%
%plotting of Crank Nicholson
x1=0:step_size_array(1):42;
x2=0:step_size_array(2):42;
x3=0:step_size_array(3):42;
x=0:0.1:42;
subplot(4,1,1)
plot(x1,phi_1);
title('Step Size=0.1');
subplot(4,1,2)
plot(x2,phi_2);
title('Step Size=0.6');
subplot(4,1,3)
plot(x3,phi_3);
title('Step Size=2.1');
subplot(4,1,4)
plot(x,phi_exact);
title('Exact Graph');
suptitle('Crank Nicholson Subplot')
figure;
plot(x1,phi_1,'.',x2,phi_2,'--',x3,phi_3,'-.',x,phi_exact,'-')
legend('step size=0.1','step size=0.6','step size=2.1','exact graph')
title('Crank Nicholson Graph')
