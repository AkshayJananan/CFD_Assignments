import numpy as np
def discrtn(u,v,w,nx,ny,Re):
    dx = 1 / (nx - 1)
    dy = 1 / (ny - 1)
    f=np.zeros((ny,nx))
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            f[j,i]=(-0.5*u[j,i]*(w[j,i+1]-w[j,i-1])/dx)+(-0.5*v[j,i]*(w[j+1,i]-w[j-1,i])/dy)+(
                (1/Re)*((w[j,i+1]-2*w[j,i]+w[j,i-1]/dx**2)+(w[j+1,i]-2*w[j,i]+w[j-1,i]/dy**2)))
    return f