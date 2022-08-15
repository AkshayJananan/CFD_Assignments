import numpy as np
def GS(psi_old,w,nx,ny,error_tol):
    dx = 1 / (nx - 1)
    dy = 1 / (ny - 1)
    l=1
    psi=psi_old.copy()
    while l>error_tol:
       sum = 0
       for i in range(1,nx-1):
         for j in range(1,ny-1):
             psi[j,i]=0.25*(psi_old[j , i+1]+psi_old[j+1 , i]+psi[j , i-1]+psi[j-1 , i]+w[j,i]*dx**2)
             sum = sum + ((psi[j, i] - psi_old[j, i]) ** 2)
       ##ERROR CALCULATION
       l=np.sqrt(sum)
       psi_old = psi.copy()
    return psi


