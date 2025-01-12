# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 09:07:48 2025

@author: C10381149
"""
import numpy as np
import pandas as pd
from math import pi
from numpy.linalg import inv

def coupled_cross_diffusion_v7(

end_exp=1e2, # End of exposure
total_time=100,# Total simulation time
lpmm=1e3,# Spatial frequency
T0=50e-4,# cm
I0=5,# Intensity of first recording beam
xi=0.3,# Scattering coefficient
Dm=1.6e-7,# Monomer diffusion coefficient
n_m=1.55,# Monomer refractive index
rhom=1.15,# Monomer mass density
Dp=6.35e-10,# Oligomer diffusivity ratio
rhop=1.3,# Fractional van der Waals space loss
n_p=1.56,# Oligomer refractive index
n_q=1.64,# Polymer refractive index
Gamma=1,# Rate of immobilization
Dz=1e-10,# Nanoparticle cross-diffusion ratio
epsilon_mz=0,# Cross-diffusion
epsilon_pz=13,# Cross-diffusion
epsilon_qz=13,# Cross-diffusion
wt_pc=5e-2,# Doping %
rhoz=1.74,# Nanoparticle mass density
n_z=1.366,# Nanoparticle refractive index
b0=5.05,# Binder mass
n_b=1.5,# Binder refractive inedx
rhob=1.19, # Binder mass density
lambda_probe=633e-7,# Wavelength of reconstruction beam
Delta_t=1/100,# Numerical scheme time step
Delta_x=1/20,# Numerical scheme spatial step
output_time_step=1# Output time

):

    # 1.1 --- Define parameters
    iterations_per_second=1/Delta_t# Number of iterations each second
    Nx=int(1/Delta_x) + 1# Number of spatial points
    x = np.linspace(0,1,Nx)# Non-dimensional grating distance
    n_iterations = int(total_time*iterations_per_second) +1# Total number of iterations
    r=Delta_t/Delta_x/Delta_x# Ratio of finite time step to squared finite spatial step
    m0=1# Initial grams of monomer
    t0=1 #  Reference time [s]
    Lambda=1/10/lpmm # Grating period [cm]
    j_end_exp=end_exp/Delta_t # Iteration of exposure end
    z0 = wt_pc/(1 - wt_pc)*(m0 + b0)# Initial nanoparticle mass
    m1=np.ones(Nx); m2=m1
    p1=np.zeros(Nx); p2=p1
    q1=np.zeros(Nx); q2=q1
    z1=np.ones(Nx); z2=z1
    b1=b0*np.ones(Nx)
    
    Vb = b1/rhob # cm**3
    Vm = m1*m0/rhom # cm**3
    Vp = p1*m0/rhop # cm**3
    Vq = q1*m0/rhop # cm**3
    Vz = z1*z0/rhoz # cm**3
    Vtotal=Vb+Vm+Vp+Vq+Vz # cm**3
    
    phi_m = Vm/Vtotal
    phi_b = Vb/Vtotal
    phi_p = Vp/Vtotal
    phi_q = Vq/Vtotal
    phi_z = Vz/Vtotal
    
    Lorentz_Lorenz_RHS = phi_m*(n_m*n_m - 1)/(n_m*n_m + 2) + phi_b*(n_b*n_b - 1)/(n_b*n_b + 2) + phi_p*(n_p*n_p - 1)/(n_p*n_p + 2) + phi_q*(n_q*n_q - 1)/(n_q*n_q + 2) + phi_z*(n_z*n_z - 1)/(n_z*n_z + 2)
    
    n1=np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
    
    interior_points = list(range(1, Nx-1))
    
    times_4 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2 != 0]
    
    times_2 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2==0]
    
    N0=[Delta_x/3*(n1[0] + sum(2*n1[times_2]) + sum(4*n1[times_4]) + n1[Nx-1])]
    
    n1_cos=n1*np.cos(2*pi*x)
    
    N1=[2*Delta_x/3*(n1_cos[0] + sum(2*n1_cos[times_2]) + sum(4*n1_cos[times_4]) + n1_cos[Nx-1])]
        
    sq_diff=(n1 - N0[0] - N1[0]*np.cos(2*pi*x))**2
    
    d2 = [Delta_x/3*(sq_diff[0] + sum(2*sq_diff[times_2]) + sum(4*sq_diff[times_4]) + sq_diff[Nx-1])]
    
    # 1.2 --- Non-dimensional parameters
    alpha_m=Dm*t0/Lambda/Lambda# Monomer diffusion
    F0=0.1*I0**0.3
    beta=F0*t0# Monomer comsumption
    alpha_p=Dp*t0/Lambda/Lambda# Oligomer diffusion
    alpha_z=Dz*t0/Lambda/Lambda# Nondimensional self-nanoparticle diffusion
    gamma = m0*Gamma*t0# Immobilization
    if wt_pc==0:
        alpha_mz=0
        alpha_pz=0
        alpha_zm=0
        alpha_zp=0
        alpha_zq=0
    else:
        alpha_mz=z0*epsilon_mz*alpha_m
        alpha_pz=z0*epsilon_pz*alpha_p
        alpha_zm=epsilon_mz*alpha_z
        alpha_zp=epsilon_pz*alpha_z
        alpha_zq=epsilon_qz*alpha_z
    
    df=pd.DataFrame({"x": x, 
                     "monomer": m1,
                     "short_polymer": p1,
                     "immobile_polymer": q1,
                     'nanoparticles': z1,
                     'binder': b1,
                     'refractive_index': n1,
                     "time": np.zeros(Nx)})
    
    time_vals=np.arange(0, total_time+1, output_time_step)
    
    # 1.3 --- Calculate each time step via implicit finite difference method
    for j in range(1,n_iterations):
        # Phi=1 if illumination is on, 0 otherwise
        if ((j >= 0) and (j <= j_end_exp)):
            Phi=1
        else:
            Phi = 0
            
        # Illumination pattern
        f = 1 + np.exp(-xi*z0*z1)*np.cos(2*pi*x)
        MM2 = (2 + Phi*Delta_t*beta*f)*np.identity(Nx)
        MM1 = (2 - Phi*Delta_t*beta*f)*np.identity(Nx)
        PP2 = (2 + Phi*Delta_t*gamma*p1)*np.identity(Nx)
        PP1 = (2 - Phi*Delta_t*gamma*p1)*np.identity(Nx)
        PM2 = Phi*Delta_t*beta*f*np.identity(Nx)
        PM1 = Phi*Delta_t*beta*f*np.identity(Nx)
        QQ2=2*np.identity(Nx)
        QQ1=2*np.identity(Nx)
        QP2=(Phi*gamma*Delta_t*p1)*np.identity(Nx)
        QP1=(Phi*gamma*Delta_t*p1)*np.identity(Nx)
        ZZ2 = 2*np.identity(Nx)
        ZZ1 = 2*np.identity(Nx)
        
        for i in range(Nx):
            if i==0:
                i_minus_1=i+1
            else:
                i_minus_1=i-1
                
            if i==Nx-1:
                i_plus_1=i-1
            else:
                i_plus_1=i+1
            
            MM2[i, i_minus_1] = MM2[i, i_minus_1] - r*alpha_m - r*alpha_mz*z1[i_minus_1]
            MM2[i, i] = MM2[i, i] + 2*r*alpha_m + 2*r*alpha_mz*z1[i]
            MM2[i, i_plus_1] = MM2[i, i_plus_1] - r*alpha_m - r*alpha_mz*z1[i_plus_1]
    
            MM1[i, i_minus_1] = MM1[i, i_minus_1] + r*alpha_m + r*alpha_mz*z1[i_minus_1]
            MM1[i, i] = MM1[i, i] - 2*r*alpha_m - 2*r*alpha_mz*z1[i]
            MM1[i, i_plus_1] = MM1[i, i_plus_1] + r*alpha_m + r*alpha_mz*z1[i_plus_1]
    
            PP2[i, i_minus_1] = PP2[i, i_minus_1] - r*alpha_p -  r*alpha_pz*z1[i_minus_1]
            PP2[i, i] = PP2[i, i] + 2*r*alpha_p + 2*r*alpha_pz*z1[i]
            PP2[i, i_plus_1] = PP2[i, i_plus_1] - r*alpha_p - r*alpha_pz*z1[i_plus_1]
            
            PP1[i, i_minus_1] = PP1[i, i_minus_1] + r*alpha_p + r*alpha_pz*z1[i_minus_1]
            PP1[i, i] = PP1[i, i] - 2*r*alpha_p - 2*r*alpha_pz*z1[i]
            PP1[i, i_plus_1] = PP1[i, i_plus_1] + r*alpha_p + r*alpha_pz*z1[i_plus_1]
    
            ZZ2[i, i_minus_1] = ZZ2[i, i_minus_1] - r*alpha_z - r*alpha_zq*q1[i_minus_1] - r*alpha_zp*p1[i_minus_1]  - r*alpha_zm*m1[i_minus_1]
            ZZ2[i, i] = ZZ2[i, i] + 2*r*alpha_z + 2*r*alpha_zq*q1[i] + 2*r*alpha_zp*p1[i] + 2*r*alpha_zm*m1[i]
            ZZ2[i, i_plus_1] = ZZ2[i, i_plus_1] - r*alpha_z - r*alpha_zq*q1[i_plus_1]  - r*alpha_zp*p1[i_plus_1]  - r*alpha_zm*m1[i_plus_1]
    
            ZZ1[i, i_minus_1] = ZZ1[i, i_minus_1] + r*alpha_z + r*alpha_zq*q1[i_minus_1] + r*alpha_zp*p1[i_minus_1] + r*alpha_zm*m1[i_minus_1]
            ZZ1[i, i] = ZZ1[i, i] - 2*r*alpha_z - 2*r*alpha_zq*q1[i]  - 2*r*alpha_zp*p1[i]  - 2*r*alpha_zm*m1[i]
            ZZ1[i, i_plus_1] = ZZ1[i, i_plus_1] + r*alpha_z + r*alpha_zq*q1[i_plus_1] + r*alpha_zp*p1[i_plus_1] + r*alpha_zm*m1[i_plus_1]
            
        m2 = np.matmul(inv(MM2), np.matmul(MM1,m1))
        
        p2 = np.matmul(inv(PP2), (np.matmul(PP1,p1) + np.matmul(PM2,m2) + np.matmul(PM1,m1)))
        
        q2 = np.matmul(inv(QQ2), (np.matmul(QQ1,q1) + np.matmul(QP2,p2) + np.matmul(QP1,p1)))
        
        z2 = np.matmul(inv(ZZ2), np.matmul(ZZ1,z1))
        
        m1=m2; p1=p2; q1=q2; z1=z2
        
        Vb = b1/rhob # cm**3
        Vm = m1*m0/rhom # cm**3
        Vp = p1*m0/rhop # cm**3
        Vq = q1*m0/rhop # cm**3
        Vz = z1*z0/rhoz # cm**3
        Vtotal=Vb+Vm+Vp+Vq+Vz # cm**3
        phi_m = Vm/Vtotal
        phi_b = Vb/Vtotal
        phi_p = Vp/Vtotal
        phi_q = Vq/Vtotal
        phi_z = Vz/Vtotal
    
        Lorentz_Lorenz_RHS = phi_m*(n_m*n_m - 1)/(n_m*n_m + 2) + phi_b*(n_b*n_b - 1)/(n_b*n_b + 2) + phi_p*(n_p*n_p - 1)/(n_p*n_p + 2) + phi_q*(n_q*n_q - 1)/(n_q*n_q + 2) + phi_z*(n_z*n_z - 1)/(n_z*n_z + 2)
    
        n1=np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
    
        N0_new=Delta_x/3*(n1[0] + sum(2*n1[times_2]) + sum(4*n1[times_4]) + n1[Nx-1])
        
        n1_cos=n1*np.cos(2*pi*x)
    
        N1_new=2*Delta_x/3*(n1_cos[0] + sum(2*n1_cos[times_2]) + sum(4*n1_cos[times_4]) + n1_cos[Nx-1])
        
        sq_diff=(n1 - N0_new - N1_new*np.cos(2*pi*x))**2
    
        d2_new = [Delta_x/3*(sq_diff[0] + sum(2*sq_diff[times_2]) + sum(4*sq_diff[times_4]) + sq_diff[Nx-1])]
        
        if Delta_t*j in time_vals:
            
            new_df=pd.DataFrame({"x": x, 
                             "monomer": m1,
                             "short_polymer": p1,
                             "immobile_polymer": q1,
                             'nanoparticles': z1,
                             'binder':b1,
                             'refractive_index': n1,
                             "time": np.zeros(Nx) + Delta_t*j})
            
            df=pd.concat([df,new_df])
            
            N0.append(N0_new)
            
            N1.append(N1_new)
            
            d2.append(d2_new)
    
    optical_properties_df=pd.DataFrame({'time': time_vals, 'N0': N0, 'Delta_n': [2*i for i in N1]})              
    
    output = {'spatial_profiles':df,
              'optical_properties': optical_properties_df,
              'end_exp':end_exp,
              'total_time':total_time,
              'lpmm':lpmm,
              'T0':T0,
              'I0':I0,
              'xi':xi,
              'Dm':Dm,
              'n_m':n_m,
              'rhom':rhom,
              'Dp':Dp,
              'rhop':rhop,
              'n_p':n_p,
              'n_q':n_q,
              'Gamma':Gamma,
              'Dz':Dz,
              'epsilon_mz':epsilon_mz,
              'epsilon_pz':epsilon_pz,
              'epsilon_qz':epsilon_qz,
              'wt_pc':wt_pc,
              'rhoz':rhoz,
              'n_z':n_z,
              'b0':b0,
              'n_b':n_b,
              'rhob':rhob,
              'lambda_probe':lambda_probe,
              'Delta_t':Delta_t,
              'Delta_x':Delta_x,
              'output_time_step':output_time_step,
              'Model': "coupled_cross_diffusion_v7"}
    
    return output