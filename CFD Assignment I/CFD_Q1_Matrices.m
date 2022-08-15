clear all;close all;clc;
%%
% AKSHAY J_21105012
%%
%AX=B,MATRIX FORM
x=linspace(0.1,1,10);
%%
%exact solution
y_exact=(2.7183).^x;
%%
%N=1
y1=1+2.*x;
error1=y1-y_exact;
%%
%N=2
A=[1/2 2/3;1/6 5/12];
B=[1;1/2];
X=A\B;
y2=1+X(1,1).*x+X(2,1).*x.^2;
error2=y2-y_exact;
%%
%N=3
A=[1/2 2/3 3/4;1/6 5/12 11/20;1/12 3/10 13/30];
B=[1;1/2;1/3];
X=A\B;
y3=1+X(1,1).*x+X(2,1).*x.^2+X(3,1).*x.^3;
error3=y3-y_exact;
%%
%N=4
A=[1/2 2/3 3/4 4/5;1/6 5/12 11/20 19/30;1/12 3/10 13/30 11/21;1/20 7/30 15/42 25/56];
B=[1;1/2;1/3;1/4];
X=A\B;
y4=1+X(1,1).*x+X(2,1).*x.^2+X(3,1).*x.^3+X(4,1).*x.^4;
error4=y4-y_exact;
%%
%N=5
A=[1/2 2/3 3/4 4/5 5/6;1/6 5/12 11/20 19/30 29/42;1/12 3/10 13/30 11/21 33/56;1/20 7/30 15/42 25/56 37/72;1/30 4/21 17/56 7/18 41/90];
B=[1;1/2;1/3;1/4;1/5];
X=A\B;
y5=1+X(1,1).*x+X(2,1).*x.^2+X(3,1).*x.^3+X(4,1).*x.^4+X(5,1).*x.^5;
error5=y5-y_exact;
%%
%N=7
A=[1/2 2/3 3/4 4/5 5/6 6/7 7/8;1/6 5/12 11/20 19/30 29/42 41/56 55/72;1/12 3/10 13/30 11/21 33/56 23/36 61/90;1/20 7/30 15/42 25/56 37/72 51/90 67/110;1/30 4/21 17/56 7/18 41/90 28/55 73/132;1/42 9/56 19/72 31/90 45/110 61/132 79/156;1/56 10/72 21/90 34/110 49/132 66/156 85/182];
B=[1;1/2;1/3;1/4;1/5;1/6;1/7];
X=A\B;
y7=1+X(1,1).*x+X(2,1).*x.^2+X(3,1).*x.^3+X(4,1).*x.^4+X(5,1).*x.^5+X(6,1).*x.^6+X(7,1).*x.^7;
error7=y7-y_exact;
%%
%plotting
figure;
subplot(3,2,1)
plot(x,error1)
title('N=1');
subplot(3,2,2)
plot(x,error2)
title('N=2');
subplot(3,2,3)
plot(x,error3)
title('N=3');
subplot(3,2,4)
plot(x,error4)
title('N=4');
subplot(3,2,5)
plot(x,error5)
title('N=5');
subplot(3,2,6)
plot(x,error7)
title('N=7');
figure;
plot(x,error1,'.',x,error2,'-^',x,error3,'--',x,error4,'-.',x,error5,'*',x,error7,'^');
legend('error1','error2','error3','error4','error5','error6','error7')
title('Error vs x');
xlabel('x values');
ylabel('error');
