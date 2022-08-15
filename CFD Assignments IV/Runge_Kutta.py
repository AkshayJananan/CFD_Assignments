def RKW3(u,v,w,nx,ny,n,Re,dt):
    dx = 1 / (nx - 1)
    dy = 1 / (ny - 1)
    n=n-1
    from Vorticity import discrtn
    w_1=w.copy()
    k1=dt*discrtn(u[:,:,n],v[:,:,n],w_1[:,:,n],nx,ny,Re)
    n1=int(n+8/15)
    w_1[:, :, n1] = w_1[:, :, n] + (8 / 15) * dt * k1  # 2nd RK
    k2=dt*discrtn(u[:,:,n1],v[:,:,n1],w_1[:,:,n1],nx,ny,Re)
    n2 = int(n1 + 2 / 3)
    w_1[:,:,n2]=w_1[:,:,n1]+(1/4)*dt*k1+(8/15)*dt*k2 #3rd RK
    k3=dt*discrtn(u[:,:,n2],v[:,:,n2],w_1[:,:,n2],nx,ny,Re)
    w[:,:,n+1]=w[:,:,n]+(1/4)*k1+(3/4)*k3
    return w[:,:,n+1]

