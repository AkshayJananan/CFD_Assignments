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
phi_exact(1)=1;
%%
%RK-3 Method Calculation
N=T/1;
for n=0:N
    phi_exact(n+1)=exp(-n);
end
N=T/step_size_1;
for n=1:N
    %phi* is provided as phi_x
phi_x=phi_1(n);%1st RK
k1=-phi_x;
phi_x=phi_x+((step_size_1)/(3)).*k1;%2nd RK
k1=(-5.*k1/9)-phi_x;
phi_x=phi_x+(15.*step_size_1.*k1)/(16);%3rd RK
k1=(-153.*k1/128)-phi_x;
phi_1(n+1)=phi_x+(8.*step_size_1.*k1)/(15);
end
N=T/step_size_2;
for n=1:N
    %phi* is provided as phi_x
phi_x=phi_2(n);%1st RK
k1=-phi_x;
phi_x=phi_x+((step_size_2)/(3)).*k1;%2nd RK
k1=(-5.*k1/9)-phi_x;
phi_x=phi_x+(15.*step_size_2.*k1)/(16);%3rd RK
k1=(-153.*k1/128)-phi_x;
phi_2(n+1)=phi_x+(8.*step_size_2.*k1)/(15);
end
N=T/step_size_3;
for n=1:N
    %phi* is provided as phi_x
phi_x=phi_3(n);%1st RK
k1=-phi_x;
phi_x=phi_x+((step_size_3)/(3)).*k1;%2nd RK
k1=(-5.*k1/9)-phi_x;
phi_x=phi_x+(15.*step_size_3.*k1)/(16);%3rd RK
k1=(-153.*k1/128)-phi_x;
phi_3(n+1)=phi_x+(8.*step_size_3.*k1)/(15);
end
%Plot of RK-3
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
suptitle('RK3 Subplot')
figure;
plot(x1,phi_1,'-.',x2,phi_2,'.',x3,phi_3,'-',x4,phi_exact,'^')
legend('step size=0.1','step size=0.6','step size=2.1','exact graph')
title('RK3 Graph');
