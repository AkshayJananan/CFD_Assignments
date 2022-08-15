clear all;
close all;
clc;
%%
% AKSHAY J_21105012
%%
%defining the grid
T=42;
step_size_1=0.1;%delta(t) values
step_size_2=0.6;
step_size_3=2.1;
%initial values
phi_1(1)=1;
phi_2(1)=1;
phi_3(1)=1;
%%
%calculation of implicit euler
N=T/1;
for n=0:N
    phi_exact(n+1)=exp(-n);
end
N=T/step_size_1;
for n=1:N
    phi_1(n+1)=phi_1(n)./(1+step_size_1);

end
N=T/step_size_2;
for n=1:N
    phi_2(n+1)=phi_1(n)./(1+step_size_2);
end
N=T/step_size_3;
for n=1:N
     phi_3(n+1)=phi_3(n)./(1+step_size_3);
end
%Plot of Implicit Euler
figure;
subplot(4,1,1)
x1=0:step_size_1:42;
plot(x1,phi_1);
title('step size=0.1');
subplot(4,1,2)
x2=0:step_size_2:42;
plot(x2,phi_2);
title('step size=0.6');
subplot(4,1,3)
x3=0:step_size_3:42;
plot(x3,phi_3);
title('step size=2.1');
subplot(4,1,4)
x4=0:42;
plot(x4,phi_exact);
title('exact graph');
suptitle('Implicit Euler Subplot')
figure;
plot(x1,phi_1,'-.',x2,phi_2,'.',x3,phi_3,'-',x4,phi_exact,'^')
legend('step size=0.1','step size=0.6','step size=2.1','exact graph')
title('Implicit Euler Graph');
