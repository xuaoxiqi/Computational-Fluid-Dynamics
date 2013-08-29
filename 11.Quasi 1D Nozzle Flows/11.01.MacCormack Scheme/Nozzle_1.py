
#This program solves quasi one dimensional nozzle flows using MacCormack Scheme
#Subsonic-Supersonic Isentropic Flow
#Governing Equations in Nonconservation Form
#This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
#Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

import numpy as np
import matplotlib.pyplot as plt

def mac():
    n = 61
    gamma = 1.4
    dx = 3.0/(n-1) 
    
    X = np.arange(0,3+dx,dx)
    S = np.zeros(n)
    rho = np.zeros(n)
    Te = np.zeros(n)
    P = np.zeros(n)
    V = np.zeros(n)
    Ma = np.zeros(n)

    drho = np.zeros(n)
    dV = np.zeros(n)
    dTe = np.zeros(n)
    
    Rrho = np.zeros(n)
    RV = np.zeros(n)
    RTe = np.zeros(n)
    
    for i in range(n):
        S[i] = 1+2.2*(X[i]-1.5)*(X[i]-1.5)
        rho[i] = 1-0.3146*X[i]
        Te[i] = 1-0.2314*X[i]
        V[i] = (0.1+1.09*X[i])*np.sqrt(Te[i])

    dt = 1.0
    for nt in range(1400):
    #dt
        for i in range(n-1):
            t_min = 0.5*dx/(np.sqrt(Te[i])+V[i])
            if t_min < dt:
                dt = t_min
            
    #step1
            drho[i] = -rho[i]*(V[i+1]-V[i])/dx-rho[i]*V[i]*(np.log(S[i+1])-np.log(S[i]))/dx-V[i]*(rho[i+1]-rho[i])/dx
            dV[i] = -V[i]*(V[i+1]-V[i])/dx-1.0/gamma*((Te[i+1]-Te[i])/dx + Te[i]/rho[i]*(rho[i+1]-rho[i])/dx)
            dTe[i] = -V[i]*(Te[i+1]-Te[i])/dx-(gamma-1)*Te[i]*((V[i+1]-V[i])/dx + V[i]*(np.log(S[i+1])-np.log(S[i]))/dx)
    #step2           
        for i in range(n-1):
            Rrho[i] = rho[i]+drho[i]*dt
            RV[i] = V[i]+dV[i]*dt
            RTe[i] = Te[i]+dTe[i]*dt
            
        for i in range(1,n-1):
    #step3
            dRrho = -Rrho[i]*(RV[i]-RV[i-1])/dx - Rrho[i]*RV[i]*(np.log(S[i])-np.log(S[i-1]))/dx - RV[i]*(Rrho[i]-Rrho[i-1])/dx
            dRV = -RV[i]*(RV[i]-RV[i-1])/dx - 1.0/gamma*((RTe[i]-RTe[i-1])/dx + RTe[i]/Rrho[i]*(Rrho[i]-Rrho[i-1])/dx)
            dRTe = -RV[i]*(RTe[i]-RTe[i-1])/dx - (gamma-1)*RTe[i]*((RV[i]-RV[i-1])/dx + RV[i]*(np.log(S[i])-np.log(S[i-1]))/dx)                                                            
    #step4
            drho_av = 0.5*(drho[i]+dRrho)
            dV_av = 0.5*(dV[i]+dRV)
            dTe_av = 0.5*(dTe[i]+dRTe)
    #step5
            rho[i] = rho[i]+drho_av*dt
            V[i] = V[i]+dV_av*dt
            Te[i] = Te[i]+dTe_av*dt
    #B.C.
        V[0] = 2*V[1]-V[2]
        rho[0] = 1.0
        Te[0] = 1.0

        rho[n-1] = 2*rho[n-2]-rho[n-3]
        V[n-1] = 2*V[n-2]-V[n-3]
        Te[n-1] = 2*Te[n-2]-Te[n-3]

        for i in range(n):
            P[i] = rho[i]*Te[i]
            Ma[i] = V[i]/np.sqrt(Te[i])
  
#Show results

    plt.plot(X,rho,'ro-',label="MacCormack Scheme",linewidth=1)
    #plt.plot(X,pre_rho,'k--',label="Exact Solution")
    plt.xlabel(r'$X$')
    plt.ylabel(r'$\rho$')   
    plt.title('Quasi 1D Nozzle Flows'+' (Density)')
    plt.legend(prop={"size":12})
    plt.savefig('1_rho')
    plt.show()
    
    plt.plot(X,Te,'ro-',label="MacCormack Scheme",linewidth=1)
    #plt.plot(X,pre_Te,'k--',label="Exact Solution")
    plt.xlabel(r'$X$')
    plt.ylabel(r'$T$')   
    plt.title('Quasi 1D Nozzle Flows'+' (Temperature)')
    plt.legend(prop={"size":12})
    plt.savefig('1_T')
    plt.show()

    plt.plot(X,P,'ro-',label="MacCormack Scheme",linewidth=1)
    #plt.plot(X,pre_Te,'k--',label="Exact Solution")
    plt.xlabel(r'$X$')
    plt.ylabel(r'$P$')   
    plt.title('Quasi 1D Nozzle Flows'+' (Pressure)')
    plt.legend(prop={"size":12})
    plt.savefig('1_P')
    plt.show()

    plt.plot(X,Ma,'ro-',label="MacCormack Scheme",linewidth=1)
    #plt.plot(X,pre_V,'k--',label="Exact Solution")
    plt.xlabel(r'$X$')
    plt.ylabel(r'$Ma$')   
    plt.title('Quasi 1D Nozzle Flows'+' (Mach number)')
    plt.legend(prop={"size":12})
    plt.savefig('1_Ma')
    plt.show()
    
    return

    
def main():
    print 'This program solves quasi one dimensional nozzle flows using MacCormack Scheme'
    print 'Subsonic-Supersonic Isentropic Flow'
    print 'Governing Equations in Nonconservation Form'
    mac()
    print 'End!'


if __name__=='__main__':
    main()
