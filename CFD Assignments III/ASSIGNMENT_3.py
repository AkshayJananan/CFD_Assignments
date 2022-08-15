# AKSHAY J
#ROLL NO:21105012
#ASSIGNMENT-3
import numpy
from matplotlib import pyplot
from matplotlib.animation import FuncAnimation
#Setting the cell
nx=101
nt=1500
dx=2/(nx-1)
Re1=10
Re2=50
dt=0.05*dx
dt1=0.10*dx
Xmax=1
x1=numpy.array(numpy.linspace(-Xmax,Xmax,nx))
x=x1.tolist()
u_e1=numpy.empty((nx,nt))
u_e2=numpy.empty((nx,nt))
u_i1=numpy.empty((nx,nt))
u_i2=numpy.empty((nx,nt))
####### Explicit
# Intital Condition
u_e1[x.index(-Xmax):x.index(0), 0] = 1
u_e1[x.index(0) + 1:x.index(Xmax), 0] = 0
u_e2[x.index(-Xmax):x.index(0), 0] = 1
u_e2[x.index(0) + 1:x.index(Xmax), 0] = 0
#Calculation
for n in range(nt-1):
   for i in range(1, nx - 1):
      u_e1[i,n+1] = u_e1[i,n] - (u_e1[i,n] * dt / dx) * (u_e1[i,n] - u_e1[i - 1,n]) + (dt / (Re1 * dx ** 2)) * (
                    u_e1[i + 1,n] - 2 * u_e1[i,n] + u_e1[i - 1,n])
      u_e2[i, n + 1] = u_e2[i, n] - (u_e2[i, n] * dt / dx) * (u_e2[i, n] - u_e2[i - 1, n]) + (dt / (Re2 * dx ** 2)) * (
              u_e2[i + 1, n] - 2 * u_e2[i, n] + u_e2[i - 1, n])
      # Boundary Condition
      u_e1[x.index(-Xmax),:] = 1
      u_e1[x.index(Xmax),:] = 0
      u_e2[x.index(-Xmax), :] = 1
      u_e2[x.index(Xmax), :] = 0
#Graph of Explicit for time=0.36,0.47
pyplot.subplot(2,1,1)
pyplot.plot(x,u_e1[:,360],label='Re=10',color='r')
pyplot.plot(x,u_e2[:,360],label='Re=50',color='g')
pyplot.ylabel('U Values-->')
pyplot.title('Explicit Solution at t=0.36')
pyplot.legend()
pyplot.subplot(2,1,2)
pyplot.plot(x,u_e1[:,470],label='Re=10',color='r')
pyplot.plot(x,u_e2[:,470],label='Re=50',color='g')
pyplot.xlabel('X Values-->')
pyplot.ylabel('U Values-->')
pyplot.title('Explicit Solution at t=0.47')
pyplot.legend()
pyplot.show()
#Animation of Explicit plot
fig=pyplot.figure()
axis=pyplot.axes(xlim=(-1.0,1.0),ylim=(0,1.5))
axis2=pyplot.axes(xlim=(-1.0,1.0),ylim=(0,1.5))
line1,=axis.plot([],[],lw=3,label='Re=10',color='r')
line2,=axis2.plot([],[],lw=3,label='Re=50',color='g')
axis2.legend([line1,line2],[line1.get_label(),line2.get_label()])
def init():
     line1.set_data([],[])
     line2.set_data([],[])
     return line1,line2,
def animate1(i):
      line1.set_data(x,u_e1[:,i])
      line2.set_data(x,u_e2[:,i])
      return line1,line2,
anim=FuncAnimation(fig,animate1,init_func=init,frames=nt,interval=1,blit=True)
pyplot.xlabel('X Values-->')
pyplot.ylabel('U Values-->')
pyplot.title('Explicit Method=FTBS')
pyplot.show()
#######Implicit Solution
a=numpy.zeros(nx)
b=numpy.zeros(nx)
c=numpy.zeros(nx)
d=numpy.zeros(nx)
s1 = dt / (Re1 * dx ** 2)
s2 = dt / (Re2 * dx ** 2)
# Intital Condition
u_i1[x.index(-Xmax):x.index(0),0] = 1
u_i1[x.index(0) + 1:x.index(Xmax),0] = 0
u_i2[x.index(-Xmax):x.index(0),0] = 1
u_i2[x.index(0) + 1:x.index(Xmax),0] = 0
r1=dt/(2*Re1*dx**2)
r2=dt/(2*Re2*dx**2)
p=dt/(4*dx)
#Calculation
for n in range(nt-1):
      for i in range(0,nx-1):
         a[i]=-(r1+p*u_i1[i,n])
         b[i]=2*r1+1+p*u_i1[i+1,n]-p*u_i1[i-1,n]
         c[i]=-r1+p*u_i1[i,n]
         d[i]=r1*u_i1[i+1,n]+(1-2*r1)*u_i1[i,n]+r1*u_i1[i-1,n]
      from MATRIX_Solver import TDMA
      u_i1[:,n+1]=TDMA(a,b,c,d,nx,u_i1[:,n])
      # Boundary Condition
      u_i1[x.index(-Xmax), :] = 1
      u_i1[x.index(Xmax), :] = 0
for n in range(nt-1):
      for i in range(1,nx-1):
          a[i] = -(r2 + p * u_i2[i, n])
          b[i] = 2 * r2 + 1 + p * u_i2[i + 1, n] - p * u_i2[i - 1, n]
          c[i] = -r2 + p * u_i2[i, n]
          d[i] = r2 * u_i2[i + 1, n] + (1 - 2 * r2) * u_i2[i, n] + r2 * u_i2[i - 1, n]
      from MATRIX_Solver import TDMA

      u_i2[:, n + 1] = TDMA(a, b, c, d, nx, u_i2[:, n])
      # Boundary Condition
      u_i2[x.index(-Xmax), :] = 1
      u_i2[x.index(Xmax), :] = 0
#Graph of Implicit Solution for time=0.36,0.47
pyplot.subplot(2,1,1)
pyplot.plot(x,u_i1[:,360],label='Re=10',color='r')
pyplot.plot(x,u_i2[:,360],label='Re=50',color='g')
pyplot.ylabel('U Values-->')
pyplot.title('Implicit Solution at t=0.36')
pyplot.legend()
pyplot.subplot(2,1,2)
pyplot.plot(x,u_i1[:,470],label='Re=10',color='r')
pyplot.plot(x,u_i2[:,470],label='Re=50',color='g')
pyplot.xlabel('X Values-->')
pyplot.ylabel('U Values-->')
pyplot.title('Implicit Solution at t=0.47')
pyplot.legend()
pyplot.show()
#Animation of Implicit
fig=pyplot.figure()
axis=pyplot.axes(xlim=(-1.0,1.0),ylim=(0,1.5))
axis2=pyplot.axes(xlim=(-1.0,1.0),ylim=(0,1.5))
line1,=axis.plot([],[],lw=3,label='Re=10',color='r')
line2,=axis2.plot([],[],lw=3,label='Re=50',color='g')
axis2.legend([line1,line2],[line1.get_label(),line2.get_label()])
def init():
     line1.set_data([],[])
     line2.set_data([],[])
     return line1,line2,
def animate1(i):
      line1.set_data(x,u_i1[:,i])
      line2.set_data(x,u_i2[:,i])
      return line1,line2,
anim=FuncAnimation(fig,animate1,init_func=init,frames=nt,interval=20,blit=True)
pyplot.xlabel('X Values-->')
pyplot.ylabel('U Values-->')
pyplot.title('Implicit Method')
pyplot.show()
#Graph of Explicit vs Implicit
pyplot.subplot(2,1,1)
pyplot.plot(x,u_e1[:,360],label='Explicit for Re=10,t=0.36',color='r')
pyplot.plot(x,u_i1[:,360],label='Implicit for Re=10,t=0.36',color='g')
pyplot.title('For t=0.36')
pyplot.legend()
pyplot.subplot(2,1,2)
pyplot.plot(x,u_e1[:,470],label='Explicit for Re=10,t=0.47',color='r')
pyplot.plot(x,u_i1[:,470],label='Implicit for Re=10,t=0.47',color='g')
pyplot.title('For t=0.47')
pyplot.xlabel('X Values-->')
pyplot.ylabel('U Values-->')
pyplot.legend()
pyplot.suptitle('Explicit vs Implicit')
pyplot.show()
#Animation of Explicit vs Implicit
fig=pyplot.figure()
axis=pyplot.axes(xlim=(-1.0,1.0),ylim=(0,1.5))
axis2=pyplot.axes(xlim=(-1.0,1.0),ylim=(0,1.5))
line1,=axis.plot([],[],lw=3,label='Explicit',color='r')
line2,=axis2.plot([],[],lw=3,label='Implicit',color='k')
axis2.legend([line1,line2],[line1.get_label(),line2.get_label()])
def init():
     line1.set_data([],[])
     line2.set_data([],[])
     return line1,line2,
def animate1(i):
      line1.set_data(x,u_e1[:,i])
      line2.set_data(x,u_i1[:,i])
      return line1,line2,
anim=FuncAnimation(fig,animate1,init_func=init,frames=nt,interval=1,blit=True)
pyplot.xlabel('X Values-->')
pyplot.ylabel('U Values-->')
pyplot.title('Explicit vs Implicit for Re=10')
pyplot.show()
####CHECKING FOR STABILITY USING DIFFERNENT CFL VALUES
u_e1=numpy.empty((nx,nt))
u_e2=numpy.empty((nx,nt))
u_i1=numpy.empty((nx,nt))
u_i2=numpy.empty((nx,nt))
####### Explicit
# Intital Condition
u_e1[x.index(-Xmax):x.index(0), 0] = 1
u_e1[x.index(0) + 1:x.index(Xmax), 0] = 0
#Calculation
for n in range(nt-1):
   for i in range(1, nx - 1):
      u_e1[i,n+1] = u_e1[i,n] - (u_e1[i,n] * dt1 / dx) * (u_e1[i,n] - u_e1[i - 1,n]) + (dt1 / (Re1 * dx ** 2)) * (
                    u_e1[i + 1,n] - 2 * u_e1[i,n] + u_e1[i - 1,n])
      # Boundary Condition
      u_e1[x.index(-Xmax),:] = 1
      u_e1[x.index(Xmax),:] = 0
#######Implicit Solution
a=numpy.zeros(nx)
b=numpy.zeros(nx)
c=numpy.zeros(nx)
d=numpy.zeros(nx)
s1 = dt1 / (Re1 * dx ** 2)
# Intital Condition
u_i1[x.index(-Xmax):x.index(0),0] = 1
u_i1[x.index(0) + 1:x.index(Xmax),0] = 0
r1=dt1/(2*Re1*dx**2)
p=dt1/(4*dx)
#Calculation
for n in range(nt-1):
      for i in range(0,nx-1):
         a[i]=-(r1+p*u_i1[i,n])
         b[i]=2*r1+1+p*u_i1[i+1,n]-p*u_i1[i-1,n]
         c[i]=-r1+p*u_i1[i,n]
         d[i]=r1*u_i1[i+1,n]+(1-2*r1)*u_i1[i,n]+r1*u_i1[i-1,n]
      from MATRIX_Solver import TDMA
      u_i1[:,n+1]=TDMA(a,b,c,d,nx,u_i1[:,n])
      # Boundary Condition
      u_i1[x.index(-Xmax), :] = 1
      u_i1[x.index(Xmax), :] = 0
pyplot.plot(x,u_e1[:,360],label='Explicit',color='r')
pyplot.plot(x,u_i1[:,360],label='Implicit',color='g')
pyplot.xlabel('X Values-->')
pyplot.ylabel('U Values-->')
pyplot.title('CFL VALUE=0.1*dx')
pyplot.legend()
pyplot.show()