%% AKSHAY J_21105012
%MAIN CODE

clear all;clc;
close all;
%% defining the mesh
for g=1:2
iterations=0;
nxp2=[43,83];%no.of cells +ghost cells 
nyp2=[43,83];
delta_x=1/(nxp2(g)-1);%step_size in x
delta_y=1/(nyp2(g)-1);%step_size in y
phi1(1:nyp2(g),1:nxp2(g))=0;%initial guesses
phi2(1:nyp2(g),1:nxp2(g))=0;
phi3(1:nyp2(g),1:nxp2(g))=0;
phi4(1:nyp2(g),1:nxp2(g))=0;
phi5(1:nyp2(g),1:nxp2(g))=0;
phi_old1(1:nyp2(g),1:nxp2(g))=0;
phi_old2(1:nyp2(g),1:nxp2(g))=0;
phi_old3(1:nyp2(g),1:nxp2(g))=0;
phi_old4(1:nyp2(g),1:nxp2(g))=0;
phi_old5(1:nyp2(g),1:nxp2(g))=0;
error_tol=10^(-6);%tolerance level
x_array=linspace(delta_x/2,1-delta_x/2,nxp2(g)-2);
y_array=linspace(delta_y/2,1-delta_y/2,nyp2(g)-2);
L1(:)=1;
L2(:)=1;
L3(:)=1;
L4(:)=1;
L5(:)=1;
L6(:)=1;
dom_size=1;%domain range
iterations=0;
%% setting the source term and analytical colution
j=1;
for y=y_array(2:end)
    j=j+1;
    i=1;
    for x=x_array(2:end)
        i=i+1;
        %source term
        s(j,i)=2*sinh(10*(x-0.5))+40*(x-0.5)*cosh(10*(x-0.5))+100*((x-0.5)^2)*sinh(10*(x-0.5))+2*sinh(10*(y-0.5))+40*(y-0.5)*cosh(10*(y-0.5))+100*((y-0.5)^2)*sinh(10*(y-0.5))+4*(x^2+y^2)*exp(2*x*y);
    end
end
j=1;
for y=y_array(1:end)
    j=j+1;
    i=1;
    for x=x_array(1:end)
        i=i+1;
        %exact value of phi
        phi_exact(j,i)=((x-0.5)^2)*sinh(10*(x-0.5))+((y-0.5)^2)*sinh(10*(y-0.5))+exp(2*x*y);
    end
end
%% UPDATING GHOST CELLS
iterations=0;
while L1>error_tol| L2>error_tol|L3>error_tol|L4>error_tol|L5>error_tol
if L1>error_tol 
 for i=2:nxp2(g)-1
     phi_old1(1,i)=2*(0.25*sinh(-5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+1)-phi_old1(2,i);
     phi_old(nyp2(g),i)=2*(0.25*sinh(5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+exp(2*x_array(i-1)))-phi_old1(nyp2(g)-1,i);
 end
for j=2:nyp2(g)-1
     phi_old1(j,1)=2*(0.25*sinh(-5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+1)-phi_old1(j,2);
     phi_old1(j,nxp2(g))=2*(0.25*sinh(5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+exp(2*y_array(j-1)))-phi_old1(j,nxp2(g)-1);
end
end
if L2>error_tol
 %UPDATING GHOST CELLS
 for i=2:nxp2(g)-1
     phi_old2(1,i)=2*(0.25*sinh(-5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+1)-phi_old2(2,i);
     phi_old2(nyp2(g),i)=2*(0.25*sinh(5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+exp(2*x_array(i-1)))-phi_old2(nyp2(g)-1,i);
 end
for j=2:nyp2(g)-1
     phi_old2(j,1)=2*(0.25*sinh(-5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+1)-phi_old2(j,2);
     phi_old2(j,nxp2(g))=2*(0.25*sinh(5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+exp(2*y_array(j-1)))-phi_old2(j,nxp2(g)-1);
end
end
if L3>error_tol
 %UPDATING GHOST CELLS
 for i=2:nxp2(g)-1
     phi_old3(1,i)=2*(0.25*sinh(-5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+1)-phi_old3(2,i);
     phi_old3(nyp2(g),i)=2*(0.25*sinh(5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+exp(2*x_array(i-1)))-phi_old3(nyp2(g)-1,i);
 end
for j=2:nyp2(g)-1
     phi_old3(j,1)=2*(0.25*sinh(-5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+1)-phi_old3(j,2);
     phi_old3(j,nxp2(g))=2*(0.25*sinh(5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+exp(2*y_array(j-1)))-phi_old3(j,nxp2(g)-1);
end
end
if L4>error_tol
 %UPDATING GHOST CELLS
 for i=2:nxp2(g)-1
     phi_old4(1,i)=2*(0.25*sinh(-5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+1)-phi_old4(2,i);
     phi_old4(nyp2(g),i)=2*(0.25*sinh(5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+exp(2*x_array(i-1)))-phi_old4(nyp2(g)-1,i);
 end
for j=2:nyp2(g)-1
     phi_old4(j,1)=2*(0.25*sinh(-5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+1)-phi_old4(j,2);
     phi_old4(j,nxp2(g))=2*(0.25*sinh(5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+exp(2*y_array(j-1)))-phi_old4(j,nxp2(g)-1);
end
end
if L5>error_tol
 %UPDATING GHOST CELLS
 for i=2:nxp2(g)-1
     phi_old5(1,i)=2*(0.25*sinh(-5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+1)-phi_old5(2,i);
     phi_old5(nyp2(g),i)=2*(0.25*sinh(5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+exp(2*x_array(i-1)))-phi_old5(nyp2(g)-1,i);
 end
for j=2:nyp2(g)-1
     phi_old5(j,1)=2*(0.25*sinh(-5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+1)-phi_old5(j,2);
     phi_old5(j,nxp2(g))=2*(0.25*sinh(5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+exp(2*y_array(j-1)))-phi_old5(j,nxp2(g)-1);
end
end
%% calling out subroutine
if L1>error_tol
phi1=jacobi(nxp2(g),nyp2(g),dom_size,s,phi_old1); 
end
 for i=2:nxp2(g)-1
     phi_old1(1,i)=2*(0.25*sinh(-5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+1)-phi_old1(2,i);
     phi_old(nyp2(g),i)=2*(0.25*sinh(5)+((x_array(i-1)-0.5)^2)*sinh(10*(x_array(i-1)-0.5))+exp(2*x_array(i-1)))-phi_old1(nyp2(g)-1,i);
 end
for j=2:nyp2(g)-1
     phi_old1(j,1)=2*(0.25*sinh(-5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+1)-phi_old1(j,2);
     phi_old1(j,nxp2(g))=2*(0.25*sinh(5)+((y_array(j-1)-0.5)^2)*sinh(10*(y_array(j-1)-0.5))+exp(2*y_array(j-1)))-phi_old1(j,nxp2(g)-1);
end
if L2>error_tol
    phi2=gauss_seidel(nxp2(g),nyp2(g),dom_size,s,phi_old2); 
end
if L3>error_tol
    phi3=SOR(nxp2(g),nyp2(g),dom_size,s,phi_old3);
end
if L4>error_tol
    phi4=RBPGS(nxp2(g),nyp2(g),dom_size,s,phi_old4);
end
if L5>error_tol
    phi5=ADI(nxp2(g),nyp2(g),dom_size,s,phi_old5);
end
%% Residual Error
r1=0;
r2=0;
r3=0;
r4=0;
r5=0;
iterations=iterations+1;
for j=2:nyp2(g)-1
        for i=2:nxp2(g)-1
            r1=r1+((phi1(j,i)-phi_old1(j,i)).^2);
            r2=r2+((phi2(j,i)-phi_old2(j,i)).^2);
            r3=r3+((phi3(j,i)-phi_old3(j,i)).^2);
            r4=r4+((phi4(j,i)-phi_old4(j,i)).^2);
            r5=r5+((phi5(j,i)-phi_old5(j,i)).^2);
        end
end
L1(iterations)=abs(sqrt(r1));
L2(iterations)=abs(sqrt(r2));
L3(iterations)=abs(sqrt(r3));
L4(iterations)=abs(sqrt(r4));
L5(iterations)=abs(sqrt(r5));
if L1> error_tol
    phi_old1(:,:)=phi1(:,:);
else
    break
end
if L2> error_tol
    phi_old2(:,:)=phi2(:,:);
else
    break
end
if L3> error_tol
    phi_old3(:,:)=phi3(:,:);
else
    break
end
if L4> error_tol
    phi_old4(:,:)=phi4(:,:);
else
    break
end
if L5> error_tol
    phi_old5(:,:)=phi5(:,:);
else
    break
end
end
%Error=phi 
error1=phi1(2:nyp2(g)-1,2:nxp2(g)-1)-phi_exact(2:end,2:end);
error2=phi2(2:nyp2(g)-1,2:nxp2(g)-1)-phi_exact(2:end,2:end);
error3=phi3(2:nyp2(g)-1,2:nxp2(g)-1)-phi_exact(2:end,2:end);
error4=phi4(2:nyp2(g)-1,2:nxp2(g)-1)-phi_exact(2:end,2:end);
error5=phi5(2:nyp2(g)-1,2:nxp2(g)-1)-phi_exact(2:end,2:end);
%% Plotting
%Residual Error
figure;
hold on
plot((1:20),L1(1:20),'b',(1:20),L2(1:20),'r',(1:20),L3(1:20),'k',(1:20),L4(1:20),'g',(1:20),L5(1:20),'m');
title("Residuals of different methods(Grid size=" +nxp2(g)+ ")");
legend('Jacobi','GAUSS-SEIDEL','SOR','RBPGS','ADI');
xlabel('iterations-->');
ylabel('residual error-->');
hold off
[X,Y]=meshgrid(x_array,y_array);
figure;
%ERROR contours
subplot(3,2,1)
contour(X,Y,error1);
title('Jacobi')
subplot(3,2,2)
contour(X,Y,error2);
title('Gauss-Seidel')
subplot(3,2,3)
contour(X,Y,error3);
title('SOR')
subplot(3,2,4)
contour(X,Y,error4);
title('RBPGS')
subplot(3,2,5)
contour(X,Y,error5);
title('ADI')
suptitle("ERROR CONTOUR (Grid size=" +nxp2(g)+ ")");
%Plotting-Jacobi
[X,Y]=meshgrid(x_array,y_array);
figure;subplot(3,2,1)
contour(X,Y,phi1(2:end-1,2:end-1),50)
colorbar
title('JACOBI')

%Plotting-GS
subplot(3,2,2)
contour(X,Y,phi2(2:end-1,2:end-1),50)
colorbar
title('GAUSS-SEIDEL')

%Plotting-SOR
subplot(3,2,3)
contour(X,Y,phi3(2:end-1,2:end-1),50)
colorbar
title('SUCCESSIVE OVER RELAXATION')

%Plotting-RBPGS
subplot(3,2,4)
contour(X,Y,phi4(2:end-1,2:end-1),50)
colorbar
title('RED BLACK PGS')

%Plotting-ADI
subplot(3,2,5)
contour(X,Y,phi5(2:end-1,2:end-1),50)
colorbar
title('(ADI)')

%%Plotting-Analytical solution
subplot(3,2,6)
contour(X,Y,phi_exact(1:end-1,1:end-1),50);
colorbar
title('EXACT')

suptitle("Phi Contours (Grid size=" +nxp2(g)+ ")");
end