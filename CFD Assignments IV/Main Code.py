####AKSHAY J
##21105012
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
###SETTING THE CELL
nx=64
ny=64
nt=60
dx=1/(nx-1)
dy=1/(ny-1)
dt=0.0001
Re1=10000
error_tol=10**-2
###INITIALSING
psi=np.random.uniform(-1,1,(ny,nx,nt))
u=np.zeros((ny,nx,nt))
v=np.zeros((ny,nx,nt))
w=np.random.uniform(-1,1,(ny,nx,nt))
from Velocity import velocity
u[:,:,0], v[:,:,0] = velocity(psi[:,:,0], u[:,:,0], v[:,:,0], nx, ny)
for j in range(0,ny-1):
   for i in range(0,nx-1):
      w[j,i,0]=((v[j,i+1,0]-v[j,i,0])/dx)-((u[j+1,i,0]-u[j,i,0])/dy)
for n in range(1,nt-1):
   ###SOLVING FOR w(n+1)
   from Runge_Kutta import RKW3
   w[:,:,n+1]=RKW3(u,v,w,nx,ny,n,Re1,dt)
   #Updating BC
   w[0,:,:]=w[ny-2,:,:]
   w[ny-1,:,:]=w[1,:,:]
   w[:,0,:]=w[:,nx-2,:]
   w[:,nx-1,:]=w[:,1,:]
   ###SOLVING FOR psi(n+1)
   from Gauss_Seidel import GS
   psi[:,:,n+1]=GS(psi[:,:,n],w[:,:,n+1],nx,ny,error_tol)
   ###FINDING VELOCITY FROM PSI
   u[:,:,n+1],v[:,:,n+1]=velocity(psi[:,:,n+1],u[:,:,n],v[:,:,n],nx,ny)
   ##UPDATING BC
   u[0, :,:] = u[ny - 2, :,:]
   u[ny - 1, :,:] = u[1, :,:]
   u[:, 0,:] = u[:, nx - 2,:]
   u[:, nx - 1,:] = u[:, 1,:]
   v[0, :,:] = v[ny - 2, :,:]
   v[ny - 1, :,:] = v[1, :,:]
   v[:, 0,:] = v[:, nx - 2,:]
   v[:, nx - 1,:] = v[:, 1,:]
   psi[0, :, :] = psi[ny - 2, :, :]
   psi[ny - 1, :, :] = psi[1, :, :]
   psi[:, 0, :] = psi[:, nx - 2, :]
   psi[:, nx - 1, :] = psi[:, 1, :]
###ANIMATION OF TURBULENT FLOW
x=np.linspace(0,1,nx)
y=np.linspace(0,1,ny)
[X,Y]=np.meshgrid(x,y)
fig, ax = plt.subplots()
ax=plt.axes(xlim=(0,1.0),ylim=(0,1))
def animate(i):
   ax.clear()
   ax.contourf(X,Y,w[:,:,i])
ani = animation.FuncAnimation(fig, animate, frames=nt-1, interval=100, blit=False)
plt.show()