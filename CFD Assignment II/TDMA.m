%% AKSHAY J_21105012
function [phi] = TDMA(a,b,c,d,N)
%% TDMA
%% setting initial value
c(1)=c(1)/b(1);
d(1)=d(1)/b(1);
for i=2:N
    c(i)=c(i)/(b(i)-a(i)*c(i-1));
    d(i)=(d(i)-a(i)*d(i-1))/(b(i)-a(i)*c(i-1));
end
%% Back substitution
phi(N)=d(N);
for i=  N-1:-1:1
    phi(i)=d(i)-phi(i+1)*c(i);
end
end

