def velocity(psi,u,v,nx,ny):
    dx=1/(nx-1)
    dy=1/(ny-1)
    for i in range(1,nx-1):
        for j  in range(1,ny-1):
            u[j,i]=0.5*(psi[j+1,i]-psi[j-1,i])/dy
            v[j, i] = -0.5 * (psi[j, i+1] - psi[j, i-1]) / dx

    return u,v
