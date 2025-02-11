# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 19:01:32 2025

@author: C10381149
"""

import numpy as np
import pandas as pd
from math import pi, factorial
from numpy.linalg import inv

# 1.1 --- Define function inputs and default values.
def slanted_grating_simulation_v21(

start_exp=0,# Start of exposure
end_exp=1e2, # End of exposure
total_time=1e2,# Total simulation time
lpmm=1e3,# Spatial frequency
I0=5,# Intensity of recording beam
slant_angle=1e-4,# Grating slant_angle
xi=0.3,# Scattering coefficient
n_m=1.55,# Monomer refractive index
rhom=1.15,# Monomer density
Dm=1.6e-7,# Monomer diffusion coefficient
Dp=6.35e-10,# Polymer diffusion coefficient
rhop=1.3,# Polymer density
n_p=1.56,# Oligomer refractive index
n_q=1.64,# Polymer refractive index
Gamma=1,# Rate of immobilization
wt_pc=5e-2,# Doping %
Dz=1e-10,# Nanodopant self-diffusion coefficient
epsilon_pz=13,# Cross-diffusion ratio
epsilon_qz=13,# Cross-diffusion ratio
rhoz=1.74,# Nanodopant mass density
n_z=1.366,# Nanodopant refractive index
b0=5.05,# Ratio of binder to monomer mass
n_b=1.5,# Binder refractive index
rhob=1.19, # Binder mass density
T0=50e-4,# Depth of photosensitive layer [cm]
zeta=139,# absorption coefficient [cm**-1]
lambda_probe=633e-7,# Wavelength of reconstruction beam
Delta_t=1/100,# Numerical scheme time step
Delta_x=1/20,# Numerical scheme spatial step
output_time_step=1):# Output time
  
    # 1.2 --- Define parameters
    Delta_Y=Delta_x
    Nx=int(1/Delta_x) + 1# Number of spatial points
    Ny=int(1/Delta_Y) + 1# Number of spatial points
    if Nx%2==0:
       return "Number of x mesh points must be an odd number."
    if Ny%2==0:
        return "Number of y mesh points must be an odd number."
    
    x=np.linspace(0,1,Nx)# Non-dimensional grating distance
    n_iterations = int(total_time/Delta_t)+1# Total number of iterations
    r=Delta_t/Delta_x/Delta_x# Ratio of finite time step to squared finite spatial step
    m0=1# Initial mass of monomer
    t0=1 #  Reference time [s]
    Lambda0=1/10/lpmm # Grating period [cm]
    Lambda1=Lambda0
    j_end_exp=end_exp/Delta_t # Iteration of exposure end
    z0 = wt_pc/(1 - wt_pc)*(m0 + b0)# Initial nanodopant to monomer
    
    # 1.3 --- Matrix initial conditions
    u1=1
    du_dt=0
    m1 = np.ones(Nx*Nx)# m at j=0
    p1 = np.zeros(Nx*Nx)# p at j=0
    q1 = np.zeros(Nx*Nx)# q at j=0
    z1 = np.ones(Nx*Nx)# z at j=0
    b1 = b0*np.ones(Nx*Nx)# b at j=0
    
    Volume0=m0/rhom + b0/rhob + z0/rhoz
    
    phi_m0=m0/rhom/Volume0
    phi_z0=z0/rhoz/Volume0
    phi_b0=b0/rhob/Volume0
    
    Lorentz_Lorenz_RHS = phi_m0*(n_m*n_m - 1)/(n_m*n_m + 2) + phi_b0*(n_b*n_b - 1)/(n_b*n_b + 2) + phi_z0*(n_z*n_z - 1)/(n_z*n_z + 2)
    Initial_RI = np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
    n1=Initial_RI*np.ones(Nx*Nx)
    if slant_angle==0:
        slant_angle=1e-05
    phi_0=np.arcsin(np.sin(slant_angle/180*pi)/Initial_RI)
    phi_1=phi_0
    theta_B0=np.arcsin(lambda_probe/2/Initial_RI/Lambda0) - phi_0
    y_hat0=Lambda0/np.sin(phi_0)
    y_hat1=y_hat0
    x_hat=Lambda0/np.cos(phi_0)
    
    # 1.5 --- Nondimensionalized parameters
    alpha_m_x=Dm*t0/x_hat/x_hat
    alpha_m_y=Dm*t0/T0/T0
    alpha_p_x=Dp*t0/x_hat/x_hat
    alpha_p_y=Dp*t0/T0/T0
    alpha_z_x=Dz*t0/x_hat/x_hat
    alpha_z_y=Dz*t0/T0/T0
    F0=0.1*I0**0.3
    beta=F0*t0
    gamma=Gamma*m0*t0
    zeta_star=zeta*T0
    xi_star=xi*z0
    interior_points = list(range(1,Nx-1))
    times_4 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2 != 0]
    times_2 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2 == 0]
    Y=np.arange(0,1+Delta_x, Delta_x)
    time=np.arange(1,total_time+output_time_step,output_time_step)
    Y1=[]
    x1=[]
    for i in range(Nx):
        for j in list(Y):
            Y1.append(j)
        for j in x:
            x1.append(j)
    Y1=np.sort(Y1)    
    indexDF=pd.DataFrame({'x':x1, 'y':Y1})
    
    spatial_profile_DF=pd.DataFrame({"x": x1, 'Y': Y1, "monomer": m1, "short_polymer": p1, "immobile_polymer": q1, 'nanoparticles': z1, 'binder':b1, 'refractive_index': n1, "time": np.zeros(len(n1)), 'N0':np.zeros(len(n1))+Initial_RI, 'Delta_n':np.zeros(len(n1)), 'd2':np.zeros(len(n1))})
    
    optical_properties_DF=pd.DataFrame({'Y': Y, "time": np.zeros(len(Y)), 'N0':np.zeros(len(Y))+Initial_RI,'Delta_n':np.zeros(len(Y)),'nu':np.zeros(len(Y)),'d2':np.zeros(len(Y))})
    
    shrinkage_DF=pd.DataFrame({'time':[0], 'actual_shrinkage':[0],'phi_t':[phi_0],'theta_B':[theta_B0],'apparent_shrinkage':[0]})
    
    def trapezoidal_rule_integration(arr):
        arr=arr.reshape(Nx,Nx).T
        return Delta_x*Delta_x/4*(arr[0,0] + arr[Nx-1,0] + arr[0,Nx-1] + arr[Nx-1,Nx-1] + np.sum(2*arr[0,1:(Nx-1)]) + np.sum(2*arr[1:(Nx-1),0]) + np.sum(2*arr[Nx-1,1:(Nx-1)]) + np.sum(2*arr[1:(Nx-1),Nx-1]) + np.sum(4*arr[1:(Nx-1),1:(Nx-1)]))
    
    
    # 1.6 --- Calculate each time step via implicit finite difference method
    for j in range(1,n_iterations):
    
        if j <= j_end_exp:
            Phi=1 
        else:
            Phi = 0# Phi=1 if illumination is on, 0 otherwise
        
        f = np.zeros(Nx*Nx).reshape(Nx,Nx)
        for i in range(Nx):
            matrix_z1=z1.reshape(Nx,Nx)
            z1_i=matrix_z1[:,i]
            f[:,i] = np.exp(-0.3*zeta_star*u1*(1-Y[i]))*(1 + np.exp(-xi_star*z1_i)*np.cos(2*pi*x - 2*pi*T0/y_hat1*u1*Y[i]))
        
        f = f.reshape(Nx*Nx,)
        
        MM2 = (2 + Phi*Delta_t*beta*f)*np.identity(Nx*Nx)
        MM1 = (2 - Phi*Delta_t*beta*f)*np.identity(Nx*Nx)
        PP2 = (2 + Delta_t*gamma*p1)*np.identity(Nx*Nx)
        PP1 = (2 - Delta_t*gamma*p1)*np.identity(Nx*Nx)
        PM2 = (+Phi*Delta_t*beta*f)*np.identity(Nx*Nx)
        PM1 = (+Phi*Delta_t*beta*f)*np.identity(Nx*Nx)
        QQ2 = 2*np.identity(Nx*Nx)
        QQ1 = 2*np.identity(Nx*Nx)
        QP2 = (+Phi*gamma*Delta_t*p1)*np.identity(Nx*Nx)
        QP1 = (+Phi*gamma*Delta_t*p1)*np.identity(Nx*Nx)
        ZZ2 = 2*np.identity(Nx*Nx)
        ZZ1 = 2*np.identity(Nx*Nx)
        BB2 = 2*np.identity(Nx*Nx)
        BB1 = 2*np.identity(Nx*Nx)
        
        for i in range(Nx*Nx):
            
            if i in list(indexDF.loc[indexDF.x==0,].index):
                i_minus_1 = i+Nx-2
            else:
                i_minus_1 = i-1
                
            if i in list(indexDF.loc[indexDF.x==1,].index):
                i_plus_1 = i-Nx+2
            else:
                i_plus_1 = i+1
                
            if i in list(indexDF.loc[indexDF.y==0,].index):
                j_minus_1 = i+Nx
            else:
                j_minus_1 = i-Nx
                
            if i in list(indexDF.loc[indexDF.y==1,].index):
                j_plus_1 = i-Nx
            else:
                j_plus_1 = i+Nx
            
            MM2[i, i_minus_1] = MM2[i, i_minus_1] - r*alpha_m_x
            MM2[i, j_minus_1] = MM2[i, j_minus_1] - r*alpha_m_y/u1/u1
            MM2[i, i] = MM2[i, i] + 2*r*alpha_m_x
            MM2[i, i] = MM2[i, i] + 2*r*alpha_m_y/u1/u1
            MM2[i, i_plus_1] = MM2[i, i_plus_1] - r*alpha_m_x
            MM2[i, j_plus_1] = MM2[i, j_plus_1] - r*alpha_m_y/u1/u1
        
            MM1[i, i_minus_1] = MM1[i, i_minus_1] + r*alpha_m_x
            MM1[i, j_minus_1] = MM1[i, j_minus_1] + r*alpha_m_y/u1/u1
            MM1[i, i] = MM1[i, i] - 2*r*alpha_m_x
            MM1[i, i] = MM1[i, i] - 2*r*alpha_m_y/u1/u1
            MM1[i, i_plus_1] = MM1[i, i_plus_1] + r*alpha_m_x
            MM1[i, j_plus_1] = MM1[i, j_plus_1] + r*alpha_m_y/u1/u1
            
            PP2[i, i_minus_1] = PP2[i, i_minus_1] - r*alpha_p_x*(1 + epsilon_pz*z0*z1[i_minus_1])
            PP2[i, j_minus_1] = PP2[i, j_minus_1] - r*alpha_p_y/u1/u1*(1 + epsilon_pz*z0*z1[j_minus_1])
            PP2[i, i] = PP2[i, i] + 2*r*alpha_p_x*(1 + epsilon_pz*z0*z1[i])
            PP2[i, i] = PP2[i, i] + 2*r*alpha_p_y/u1/u1*(1 + epsilon_pz*z0*z1[i])
            PP2[i, i_plus_1] = PP2[i, i_plus_1] - r*alpha_p_x*(1 + epsilon_pz*z0*z1[i_plus_1])
            PP2[i, j_plus_1] = PP2[i, j_plus_1] - r*alpha_p_y/u1/u1*(1 + epsilon_pz*z0*z1[j_plus_1])
            
            PP1[i, i_minus_1] = PP1[i, i_minus_1] + r*alpha_p_x*(1 + epsilon_pz*z0*z1[i_minus_1])
            PP1[i, j_minus_1] = PP1[i, j_minus_1] + r*alpha_p_y/u1/u1*(1 + epsilon_pz*z0*z1[j_minus_1])
            PP1[i, i] = PP1[i, i] - 2*r*alpha_p_x*(1 + epsilon_pz*z0*z1[i])
            PP1[i, i] = PP1[i, i] - 2*r*alpha_p_y/u1/u1*(1 + epsilon_pz*z0*z1[i])
            PP1[i, i_plus_1] = PP1[i, i_plus_1] + r*alpha_p_x*(1 + epsilon_pz*z0*z1[i_plus_1])
            PP1[i, j_plus_1] = PP1[i, j_plus_1] + r*alpha_p_y/u1/u1*(1 + epsilon_pz*z0*z1[j_plus_1])
            
            ZZ2[i, i_minus_1] = ZZ2[i, i_minus_1] - r*alpha_z_x*(1 + epsilon_qz*q1[i_minus_1] + epsilon_pz*p1[i_minus_1])
            ZZ2[i, j_minus_1] = ZZ2[i, j_minus_1] - r*alpha_z_y/u1/u1*(1 + epsilon_qz*q1[j_minus_1] + epsilon_pz*p1[j_minus_1])
            ZZ2[i, i] = ZZ2[i, i] + 2*r*alpha_z_x*(1 + epsilon_qz*q1[i] + epsilon_pz*p1[i])
            ZZ2[i, i] = ZZ2[i, i] + 2*r*alpha_z_y/u1/u1*(1 + epsilon_qz*q1[i] + epsilon_pz*p1[i])
            ZZ2[i, i_plus_1] = ZZ2[i, i_plus_1] - r*alpha_z_x*(1 + epsilon_qz*q1[i_plus_1] + epsilon_pz*p1[i_plus_1])
            ZZ2[i, j_plus_1] = ZZ2[i, j_plus_1] - r*alpha_z_y/u1/u1*(1 + epsilon_qz*q1[j_plus_1] + epsilon_pz*p1[j_plus_1])
            
            ZZ1[i, i_minus_1] = ZZ1[i, i_minus_1] + r*alpha_z_x*(1 + epsilon_qz*q1[i_minus_1] + epsilon_pz*p1[i_minus_1])
            ZZ1[i, j_minus_1] = ZZ1[i, j_minus_1] + r*alpha_z_y/u1/u1*(1 + epsilon_qz*q1[j_minus_1] + epsilon_pz*p1[j_minus_1])
            ZZ1[i, i] = ZZ1[i, i] - 2*r*alpha_z_x*(1 + epsilon_qz*q1[i] + epsilon_pz*p1[i])
            ZZ1[i, i] = ZZ1[i, i] - 2*r*alpha_z_y/u1/u1*(1 + epsilon_qz*q1[i] + epsilon_pz*p1[i])
            ZZ1[i, i_plus_1] = ZZ1[i, i_plus_1] + r*alpha_z_x*(1 + epsilon_qz*q1[i_plus_1] + epsilon_pz*p1[i_plus_1])
            ZZ1[i, j_plus_1] = ZZ1[i, j_plus_1] + r*alpha_z_y/u1/u1*(1 + epsilon_qz*q1[j_plus_1] + epsilon_pz*p1[j_plus_1])
            
            MM2[i, j_minus_1] = MM2[i, j_minus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            MM1[i, j_minus_1] = MM1[i, j_minus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            MM2[i, j_plus_1] = MM2[i, j_plus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            MM1[i, j_plus_1] = MM1[i, j_plus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            
            PP2[i, j_minus_1] = PP2[i, j_minus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            PP1[i, j_minus_1] = PP1[i, j_minus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            PP2[i, j_plus_1] = PP2[i, j_plus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            PP1[i, j_plus_1] = PP1[i, j_plus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            
            QQ2[i, j_minus_1] = QQ2[i, j_minus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            QQ1[i, j_minus_1] = QQ1[i, j_minus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            QQ2[i, j_plus_1] = QQ2[i, j_plus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            QQ1[i, j_plus_1] = QQ1[i, j_plus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            
            ZZ2[i, j_minus_1] = ZZ2[i, j_minus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            ZZ1[i, j_minus_1] = ZZ1[i, j_minus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            ZZ2[i, j_plus_1] = ZZ2[i, j_plus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            ZZ1[i, j_plus_1] = ZZ1[i, j_plus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            
            BB2[i, j_minus_1] = BB2[i, j_minus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            BB1[i, j_minus_1] = BB1[i, j_minus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            BB2[i, j_plus_1] = BB2[i, j_plus_1] - Y1[i]/u1*Delta_t/2/Delta_x*du_dt
            BB1[i, j_plus_1] = BB1[i, j_plus_1] + Y1[i]/u1*Delta_t/2/Delta_x*du_dt
        
        m2 = np.matmul(inv(MM2), np.matmul(MM1,m1))
        
        p2 = np.matmul(inv(PP2), np.matmul(PP1,p1) + np.matmul(PM2,m2) + np.matmul(PM1,m1))
        
        q2 = np.matmul(inv(QQ2), np.matmul(QQ1, q1) + np.matmul(QP2,p2) + np.matmul(QP1,p1))
        
        if z0==0:
            z2 = z1
        else:
            z2 = np.matmul(inv(ZZ2), np.matmul(ZZ1,z1))
        
        b2 = np.matmul(inv(BB2), np.matmul(BB1,b1))
        
        Vb = b2/rhob # cm**3
        Vm = m2*m0/rhom # cm**3
        Vp = p2*m0/rhop # cm**3
        Vq = q2*m0/rhop # cm**3
        Vz = z2*z0/rhoz # cm**3
        Vtotal=Vb+Vm+Vp+Vq+Vz # cm**3
        phi_m = Vm/Vtotal
        phi_b = Vb/Vtotal
        phi_p = Vp/Vtotal
        phi_q = Vq/Vtotal
        phi_z = Vz/Vtotal
        
        Lorentz_Lorenz_RHS = phi_m*(n_m*n_m - 1)/(n_m*n_m + 2) + phi_b*(n_b*n_b - 1)/(n_b*n_b + 2) + phi_p*(n_p*n_p - 1)/(n_p*n_p + 2) + phi_q*(n_q*n_q - 1)/(n_q*n_q + 2) + phi_z*(n_z*n_z - 1)/(n_z*n_z + 2)
        
        n2=np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
        n2=n2.reshape(Nx,Nx).T
        
        N0=[(Delta_x/3*(n2[0,i] + np.sum(2*n2[times_2,i]) + np.sum(4*n2[times_4,i]) + n2[(Nx-1),i])) for i in range(Nx)]
        
        n2_cos=np.zeros(n2.shape)
        for i in range(Nx):
            n2_cos[:,i] = n2[:,i]*np.cos(2*pi*x)
            
        n2_sin=np.zeros(n2.shape)
        for i in range(Nx):
            n2_sin[:,i] = n2[:,i]*np.sin(2*pi*x)
        
        N1_a=np.array([2*Delta_x/3*(n2_cos[0,i] + np.sum(2*n2_cos[times_2,i]) + np.sum(4*n2_cos[times_4,i]) + n2_cos[(Nx-1),i]) for i in range(Nx)])
        
        N1_b=np.array([2*Delta_x/3*(n2_sin[0,i] + np.sum(2*n2_sin[times_2,i]) + np.sum(4*n2_sin[times_4,i]) + n2_sin[(Nx-1),i]) for i in range(Nx)])
        
        n_tilde=np.ones((Nx,Nx))
        for i in range(Nx):
            n_tilde[:,i]=N0[i]*n_tilde[:,i] + N1_a[i]*np.cos(2*pi*x) + N1_b[i]*np.sin(2*pi*x)
        
        sq_diff=(n2 - n_tilde)**2
        
        d2=[(Delta_x/3*(sq_diff[0,i] + np.sum(2*sq_diff[times_2,i]) + np.sum(4*sq_diff[times_4,i]) + sq_diff[(Nx-1),i])) for i in range(Nx)]
        
        n2=n2.T.reshape(Nx*Nx,)
        
        Volume1=(trapezoidal_rule_integration(m2)/rhom + trapezoidal_rule_integration(p2)/rhop + trapezoidal_rule_integration(q2)/rhop + trapezoidal_rule_integration(z2)*z0/rhoz + trapezoidal_rule_integration(b2)/rhob)
        
        u1 = Volume1/Volume0
        
        phi_1=np.arctan(np.tan(phi_0)/u1)
        
        Lambda1=np.cos(phi_1)/np.cos(phi_0)*Lambda0
        
        y_hat1=Lambda1/np.sin(phi_1)
        
        Mean_RI = trapezoidal_rule_integration(n1)
        
        theta_B=np.arcsin(lambda_probe/2/Mean_RI/Lambda1)-phi_1
        
        Delta_theta_B=theta_B0-theta_B,
        
        Delta_n=np.sqrt(N1_a*N1_a+N1_b*N1_b)
        
        nu=pi*Delta_n*T0*u1/lambda_probe/np.cos(theta_B)
        
        m1 = m2
        p1 = p2
        q1 = q2
        z1 = z2
        b1 = b2
        n1 = n2
        
        if Delta_t*j in time:
            
            spatial_profile_DF=pd.concat([spatial_profile_DF, pd.DataFrame({"x": x1,                                                'Y': Y1,"monomer": m1,"short_polymer": p1,"immobile_polymer": q1,'nanoparticles': z1,'binder':b1,'refractive_index': n1,"time":j*Delta_t*np.ones(len(n1))})]).reset_index(drop=True)
        
            optical_properties_DF=pd.concat([optical_properties_DF,pd.DataFrame({'time':j*Delta_t,'Y':Y,'N0':N0,'Delta_n':Delta_n,'d2':d2, 'nu':nu})]).reset_index(drop=True)
            
            shrinkage_DF=pd.concat([shrinkage_DF,pd.DataFrame({'time':[Delta_t*j], 'actual_shrinkage':[1-u1],'phi_t':[phi_1], 'theta_B':[theta_B], 'apparent_shrinkage':[1 - np.tan(phi_0)/np.tan(phi_0 + Delta_theta_B)]})]).reset_index(drop=True)
    
    
    #Moharam_Young=lambda_probe*lambda_probe/Mean_RI/Delta_n/Lambda1/Lambda1/np.cos(phi_1)
    
    Klein_Cook=2*pi*lambda_probe*T0*u1/Mean_RI/Lambda1/Lambda1/np.cos(phi_1)
    
    if Klein_Cook < 10:
        Geometry='Planar'
        J0=optical_properties_DF.nu/2
        J1=J0
        for l in range(1,101):
            J1 = J1 + ((-1)**l)/factorial(l)/factorial(l+1)*(J0**(2*l + 1))
        eta=J1*J1
    else:
        Geometry='Volume'
        eta=np.sin(np.sqrt(optical_properties_DF.nu*optical_properties_DF.nu))**2
        
    optical_properties_DF['eta']=eta
    
    output={'spatial_profiles':spatial_profile_DF,'optical_properties': optical_properties_DF,'shrinkage':shrinkage_DF,'end_exp':end_exp,'total_time':total_time,'lpmm':lpmm,'T0':T0,'I0':I0,'xi':xi,'Dm':Dm,'n_m':n_m,'rhom':rhom,'Dp':Dp,'rhop':rhop,'n_p':n_p,'n_q':n_q,'Gamma':Gamma,'Dz':Dz,'epsilon_pz':epsilon_pz,'epsilon_qz':epsilon_qz,'wt_pc':wt_pc,'rhoz':rhoz,'n_z':n_z,'b0':b0,'n_b':n_b,'rhob':rhob,'lambda_probe':lambda_probe,'Delta_t':Delta_t,'Delta_x':Delta_x,'output_time_step':output_time_step,'Geometry':Geometry,'Model': "slanted_grating_simulation_v21",'slant_angle':slant_angle,'zeta':zeta}

    return output       
