def TDMA(a,b,c,d,nx,u):
    c[0]=c[0]/b[0]
    d[0]=d[0]/b[0]
    for i in range(1,nx-1):
        c[i]=c[i]/(b[i]-a[i]*c[i-1])
        d[i]=(d[i]-a[i]*d[i-1])/(b[i]-a[i]*c[i-1])
    #Back Substitution
    u[nx-1]=d[nx-1]
    for i in range(nx-2,1,-1):
        u[i]=d[i]-u[i+1]*c[i]
    return(u)
