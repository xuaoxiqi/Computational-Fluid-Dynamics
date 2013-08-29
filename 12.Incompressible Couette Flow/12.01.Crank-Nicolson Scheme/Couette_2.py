
#This program sloves imcopressible Couette Flow using Crank-Nicolson scheme
#Neumann boundary condition
#This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
#Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

import numpy as np
import matplotlib.pyplot as plt
  
def Couette():
    n = 50   
    Re = 5000.0   # Reynold Number
    dy = 1.0/n
    dt = 12.5

    e=dt/dy/dy/Re
    temp_a = -0.5*e
    temp_b = 1+e
    bt = 1

    Y = np.arange(0,1+dy,dy)
    A = np.zeros(n)
    B = np.zeros(n)
    C = np.zeros(n)
    F = np.zeros(n)
    K = np.zeros(n+1)
    
    U = np.zeros(n+2)   #Numerical solution
    V = np.zeros(n+2)   #Exact solution
    

    for nt in range(25):   #Time develop    
        for i in range(1,n):
            K[i]=(1-e)*U[i]+0.5*e*(U[i+1]+U[i-1])

        A[0] = 0.0
        B[0] = temp_b
        C[0] = temp_a

        for i in range(1,n):
            A[i] = temp_a
            B[i] = temp_b
            C[i] = temp_a
            F[i-1] = K[i]
#Neumann B.C.
        A[n-1] = -2.*e
        B[n-1] = 2.*e
        C[n-1] = 0.
        F[n-1] = 2.*e*dy*bt

        U[1:n+1] = Thomas(A,B,C,F,n)

        
    

#Exact solution   
    for i in range(n+1):
        temp = 0.0
        for m in range(1,301):
                temp = temp + (-1)**m/(2.*m-1)/(2.*m-1)*np.exp(-(m-0.5)*(m-0.5)*np.pi*np.pi*nt*dt/Re)*np.sin((m-0.5)*np.pi*Y[i])        
        V[i] = (Y[i] + 8./np.pi/np.pi*temp)*bt

#Show results
    plt.plot(Y,U[:n+1],'ro',label="Crank-Nicolson Scheme",linewidth=1)
    plt.plot(Y,V[:n+1],'k--',label="Analytical Solution")
    plt.xlabel(r'$\frac{y}{D}$')
    plt.ylabel(r'$u/u_e$')
    plt.ylim(0.0,1.0)
    plt.title('Incompressible Couette Flow'+'   nt='+str(nt)+'$\Delta t$')
    plt.legend(prop={"size":12})
    plt.savefig('2_Couette')
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
    print 'Crank - Nicolson scheme'
    Couette()
    print 'Program End!'

if __name__=='__main__':
    main()
    
