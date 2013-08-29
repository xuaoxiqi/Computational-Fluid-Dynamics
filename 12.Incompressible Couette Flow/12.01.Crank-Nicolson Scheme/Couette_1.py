
#This program sloves imcopressible Couette Flow using Crank-Nicolson scheme
#Dirichelt boundary condition    
#This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
#Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

import numpy as np
import matplotlib.pyplot as plt
  
def Couette():
    n = 50   
    Re = 5000   # Reynold Number
    dy = 1.0/n
    dt = 12.5

    e=dt/dy/dy/Re
    temp_a = -0.5*e
    temp_b = 1+e

    Y = np.arange(0,1+dy,dy)
    A = np.zeros(n-1)
    B = np.zeros(n-1)
    C = np.zeros(n-1)
    F = np.zeros(n-1)
    K = np.zeros(n)
    
    U = np.zeros(n+1)   #Numerical solution
    V = np.zeros(n+1)   #Exact solution
#Dirichlet B.C.
    U[n] = 1.0    

    for nt in range(25):   #Time develop    
        for i in range(1,n):
            K[i]=(1-e)*U[i]+0.5*e*(U[i+1]+U[i-1])

        A[0] = 0.0
        B[0] = temp_b
        C[0] = temp_a

        for i in range(1,n-1):
            A[i] = temp_a
            B[i] = temp_b
            C[i] = temp_a
            F[i-1] = K[i]

        A[n-2] = temp_a
        B[n-2] = temp_b
        C[n-2] = 0.0
        F[n-2] = K[n-1] -temp_a

        U[1:n] = Thomas(A,B,C,F,n-1)

#Exact solution   
    for i in range(n+1):
        temp = 0
        for  m in range(1,151):
            temp = temp + 2.0*(-1)**m/np.pi/m*np.sin(m*np.pi*Y[i])*np.exp(-m*m*np.pi*np.pi*nt*dt/Re)
        V[i] = Y[i] + temp
        
#Show results
    plt.plot(Y,U,'ro',label="Crank-Nicolson Scheme",linewidth=1)
    plt.plot(Y,V,'k--',label="Analytical Solution")
    plt.xlabel(r'$\frac{y}{D}$')
    plt.ylabel(r'$u/u_e$')   
    plt.title('Incompressible Couette Flow'+'   nt='+str(nt)+'$\Delta t$')
    plt.ylim(-0.2,1.3)
    plt.legend(prop={"size":12})
    plt.savefig('1_Couette')
    plt.show()
    
    return 

        
# Thomas algorithm solve a trigonal matrix
def Thomas(A,B,C,X,n):

    C[0] = C[0]/B[0]
    X[0] = X[0]/B[0]

    for k in range(1,n):
        t = B[k] - C[k-1]*A[k]
        C[k] = C[k]/t
        X[k] = ( X[k]-X[k-1]*A[k] )/t

    for k in range(n-2, 0, -1):
            X[k] = X[k] - C[k]*X[k+1]

    return X


def main():
    print 'Incompressible Couette Flow:'
    print 'Crank-Nicolson Scheme'
    Couette()
    print 'Program End!'

if __name__=='__main__':
    main()
    
