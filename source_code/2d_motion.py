# tashihe nx , ny
import numpy as np
import matplotlib.pyplot as plt

def Results(u,v,p,nx,ny):
    
    x = np.linspace(0.,Lx,nx+1)
    y=np.linspace(0.,Ly,ny)
    X,Y=np.meshgrid(x,y)
    
    plt.figure(1)
    x = np.linspace(0.,Lx,nx+1)
    y=np.linspace(0.,Ly,ny)
    X,Y=np.meshgrid(x,y)
    plt.contourf(X,Y,u,50,cmap='jet')
    plt.contourf(X,Y,u,50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('u contours')
    plt.colorbar(orientation='horizontal')
    
    plt.figure(2)
    x = np.linspace(0.,Lx,nx)
    y=np.linspace(0.,Ly,ny+1)
    X,Y=np.meshgrid(x,y)
    plt.contourf(X,Y,v,50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('v contours')
    plt.colorbar(orientation='horizontal')
    
    plt.figure(3)
    x = np.linspace(0.,Lx,nx)
    y=np.linspace(0.,Ly,ny)
    X,Y=np.meshgrid(x,y)
    plt.contourf(X,Y,p,50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('p contours')
    plt.colorbar(orientation='horizontal')

    plt.show()


'''
def coefu(i,j,u,v,dx,dy,dt,Re,ns,nx,ny):
    I,J = i,j
    
    Fe = (u[J,i+1]+u[J-1,i+1])*dy/2.; Fw = (u[J,i]+u[J-1,i])*dy/2.
    Fn = (v[j,I]+v[j+1,I])*dx/2.; Fs = (v[j,I]+v[j-1,I])*dx/2.
    De = dy/(Re*dx);Dw = dy/(Re*dx)
    Dn = dx/(Re*dy);Ds = dx/(Re*dy)
    
    if j==0:    #bottom
        Fs = 0; Ds =  2.*dx/(Re*dy)
    elif j==ny-1:   #top
        Fn = 0; Dn = 2.*dx/(Re*dy)

    if ns == 1:
        aE = Fe/2.-De; aW = -Fw/2.-Dw
        aN = Fn/2.-Dn; aS = -Fs/2.-Ds
        aP = aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2
        b = aP0*v[J,i] 
        aP0 = dx*dy/dt
    if ns == 2:     #Hybrid
        aE = max(Fe/2.-De,Fe,0); aW = max(-Fw/2.-Dw,-Fw,0)
        aN = Fn/2.-Dn; aS = -Fs/2.-Ds
        aP = aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2
        b = aP0*v[J,i] 
        aP0 = dx*dy/dt
    
    return aP,aE,aW,aN,aS,b

def coefv(i,j,u,v,dx,dy,dt,Re,ns,nx,ny):
    I,J = i,j
    Fe = (u[J,i+1]+u[J-1,i+1])*dy/2.; Fw = (u[J,i]+u[J-1,i])*dy/2.
    Fn = (v[j,I]+v[j+1,I])*dx/2.; Fs = (v[j,I]+v[j-1,I])*dx/2.
    
    De = dy/(Re*dx);Dw = dy/(Re*dx)
    Dn = dx/(Re*dy);Ds = dx/(Re*dy)
    aP0 = dx*dy/dt
    b = aP0*v[J,i]
    
    if i == 0: #left
        Fw = 0; Dw = 2*dy/(Re*dx)
    if i == nx-1:   #right
        Fe = 0; De = 2*dy/(Re*dx)
    
    aE = Fe/2.-De; aW = -Fw/2.-Dw
    aN = Fn/2.-Dn; aS = -Fs/2.-Ds
    aP = aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2.
    
    
    return aP,aE,aW,aN,aS,b
'''
       
def xmomentom(nx,ny,dx,dy,dt,Re,us,u,v,ps,Uinlet):    #Solve x-momentom & return u_star
    ukp = us.copy()
    
    for i in range(1,nx):    
        uk = ukp.copy()
        
        for j in range(0,ny):  
            I,J = i,j
            
            Fe = (u[J,i]+u[J,i+1])*dy/2.; Fw = (u[J,i]+u[J,i-1])*dy/2.
            Fn = (v[j+1,I-1]+v[j+1,I])*dx/2.; Fs = (v[j,I-1]+v[j,I])*dx/2.
            De = dy/(Re*dx);Dw = dy/(Re*dx)
            Dn = dx/(Re*dy);Ds = dx/(Re*dy)
            aP0 = dx*dy/dt
            b = aP0*u[J,i]
            aE = Fe/2.-De; aW = -Fw/2.-Dw
            aN = Fn/2.-Dn; aS = -Fs/2.-Ds
            #print(aP0)
            if j==0:    #bottom
                Fs = 0; Ds =  2.*dx/(Re*dy)
                aS = -Fs/2.-Ds
                aP =  aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2.
                ukp[J,i] = -(aE*uk[J,i+1]+aW*ukp[J,i-1]+aN*uk[J+1,i]+(ps[J,I]-ps[J,I-1])*dy-b)/aP
            elif j==ny-1:   #top
                Fn = 0; Dn = 2.*dx/(Re*dy)
                aN = Fn/2.-Dn
                aP =  aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2.
                ukp[J,i] = -(aE*uk[J,i+1]+aW*ukp[J,i-1]+aS*ukp[J-1,i]+(ps[J,I]-ps[J,I-1])*dy-b)/aP
            else:
                aP = aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2.
                ukp[J,i] = -(aE*uk[J,i+1]+aW*ukp[J,i-1]+aN*uk[J+1,i]+aS*ukp[J-1,i]+(ps[J,I]-ps[J,I-1])*dy-b)/aP
            
            
            
            
            #if Fe/De>2 or Fw/Dw>2 or Fn/Dn>2 or Fs/Ds>2 :
                #print('wrong answer because peclect greater than two')
            
            
        ukp[:,0] = Uinlet   # left i = 1
        ukp[:,nx] = ukp[:,nx-1] # right
        
    us = ukp.copy() 
    return us

def ymomentom(nx,ny,dx,dy,dt,Re,vs,u,v,ps):    #Solve x-momentom & return u_star
    vkp = vs.copy()
    
    for i in range(0,nx):    
        
        vk = vkp.copy()
        
        for j in range(1,ny):  
            I,J = i,j
            
            Fe = (u[J,i+1]+u[J-1,i+1])*dy/2.; Fw = (u[J,i]+u[J-1,i])*dy/2.
            Fn = (v[j,I]+v[j+1,I])*dx/2.; Fs = (v[j,I]+v[j-1,I])*dx/2.
            
            De = dy/(Re*dx);Dw = dy/(Re*dx)
            Dn = dx/(Re*dy);Ds = dx/(Re*dy)
            aP0 = dx*dy/dt
            b = aP0*v[J,i]
            
            aE = Fe/2.-De; aW = -Fw/2.-Dw
            aN = Fn/2.-Dn; aS = -Fs/2.-Ds
            aP = aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2.
            if i==0 : #left
                Fw = 0; Dw = 2*dy/(Re*dx)
                aW = -Fw/2-Dw
                aP =  aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2.
                
                vkp[J,i] = -(aE*vk[J,i+1]+aS*vk[J-1,i]+aN*vk[J+1,i]+(ps[J,I]-ps[J-1,I])*dx-b)/aP
            elif i==nx-1:   #right
                Fe = 0; De = 2*dy/(Re*dx)
                aE = Fe/2-De
                aP =  aP0+De+Dw+Dn+Ds+(Fe-Fw+Fn-Fs)/2.
                vkp[J,i] = -(aN*vk[J+1,i]+aW*vk[J,i-1]+aS*vk[J-1,i]+(ps[J,I]-ps[J-1,I])*dx-b)/aP
            else:
                vkp[J,i] = -(aE*vk[J,i+1]+aW*vk[J,i-1]+aN*vk[J+1,i]+aS*vk[J-1,i]+(ps[J,I]-ps[J-1,I])*dx-b)/aP
        vkp[0,:] = 0    # bottom , j=0
        vkp[ny,:] = 0   # top 

   
    vs = vkp.copy() 
    return vs
def pressure_prime(nx,ny,dx,dy,dt,us,vs,pp): #Calculate pressure_prime
    
    ppkp = pp.copy()

    
    for j in range(0,ny):
        
        
        ppk = ppkp.copy()
        
        for i in range(0,nx):
            I,J = i,j
                
            aE = -dt*dy/dx
            aW = -dt*dy/dx
            aS = -dt*dx/dy
            aN = -dt*dx/dy
            aP = -(aE+aW+aN+aS)
            b = dy*(us[J,i+1]-us[J,i])+dx*(vs[j+1,I]-vs[j,I])
            if j == 0 and 0<i<nx-1 : #bottom
                #print(j,i)
                aS = 0
                aP = -(aE+aW+aN+aS)
                ppkp[J,I] = (-b-aW*ppk[J,I-1]-aE*ppk[J,I+1]-aN*ppk[J+1,I])/aP
            elif j == ny-1 and 0<i<nx-1 : #top
                #print('j=',J)
                aN = 0
                aP = -(aE+aW+aN+aS)
                ppkp[J,I] = (-b-aE*ppk[J,I+1]-aW*ppk[J,I-1]-aS*ppk[J-1,I])/aP
            elif i == 0 and 0<j<ny-1 :  #left
                
                aW = 0
                aP = -(aE+aW+aN+aS)
                ppkp[J,I] = (-b-aE*ppk[J,I+1]-aN*ppk[J+1,I]-aS*ppk[J-1,I])/aP
            elif i == nx-1 and 0<j<ny-1 :   #right
               
                aE = 0
                aP = -(aE+aW+aN+aS)
                ppkp[J,I] = (-b-aW*ppk[J,I-1]-aN*ppk[J+1,I]-aS*ppk[J-1,I])/aP
            elif 0<i<nx-1 and 0<j<ny-1:
                #print(J,i)
                ppkp[J,I] = (-b-aE*ppk[J,I+1]-aW*ppk[J,I-1]-aN*ppk[J+1,I]-aS*ppk[J-1,I])/aP
                #
    pp = ppkp.copy()     
    return pp 
          
def velocity_prime(nx,ny,dx,dy,dt,pp,up,vp):      #Calculate u_prime & v_prime 
    for i in range(1,nx):
        for j in range(0,ny):
            I,J = i,j
            #print(J,i)
            up[J,I] = dt/dx*(pp[J,I-1]-pp[J,I])
            
    for i in range(1,nx):
        for j in range(1,ny):
            I,J = i,j         
            vp[J,I] = dt/dy*(pp[J-1,I]-pp[J,I])
            
    return up,vp
def rsm(nx,ny,u,v,up,vp):
    s1 = abs(np.linalg.norm(u*dy)/((nx+2)*(ny+2))-np.linalg.norm(v*dx)/((nx+2)*(ny+3)))
    s2 = np.linalg.norm(up)/((nx+2)*(ny+2))
    s3 = np.linalg.norm(vp)/((nx+2)*(ny+3))
    return s1,s2,s3

#Geometry
Lx = 20.   #length in x-direction
Ly = 4.    #length in y-direction

#Mesh
nx = 70    #number of node in x-direction
ny = 15    #number of node in x-direction

dx = Lx/nx  #size of step in x-direction
dy = Ly/ny  #size of step in y-direction

dt = 0.1       #size of time step
nt = 500      #number of time step
Uinlet = 0.1

#pressure & velocity field
u = np.zeros((ny,nx+1)) ; v = np.zeros((ny+1,nx)) ; p = np.zeros((ny,nx))
up = np.zeros((ny,nx+1)) ; vp = np.zeros((ny+1,nx)) ; pp = np.zeros((ny,nx))
us = np.zeros((ny,nx+1)) ; vs = np.zeros((ny+1,nx)) ; ps = np.zeros((ny,nx))    
u[:,0] = Uinlet



Re = 25.    #Renolds number    Re=rho*U*L/mu

alpha_p = 0.8   #Relaxation factor of pressure
alpha_u = 1
alpha_v = 1

tol = 5.4e-4

maxiteration = 200

for n in range(nt*5):
    
    for iteration in range(maxiteration):
        

        us = xmomentom(nx,ny,dx,dy,dt,Re,us,u,v,ps,Uinlet)
        vs = ymomentom(nx,ny,dx,dy,dt,Re,vs,u,v,ps)
        pp = np.zeros((ny,nx))
        #for iteration1 in range(30):
        pp = pressure_prime(nx,ny,dx,dy,dt,us,vs,pp)
        
        up ,vp = velocity_prime(nx,ny,dx,dy,dt,pp,up,vp)

        
        ps = ps+pp*alpha_p
        vs = vs+vp*alpha_v
        us = us+up*alpha_u

        s1,s2,s3 = rsm(nx,ny,u,v,up,vp)
     
        
        if s1<tol :
            break
    u = us.copy()
    v = vs.copy()
    p = ps.copy()
    
    if n%5 == 0:
        print(n, iteration, s1)
            
            
Results(u,v,p,nx,ny)

