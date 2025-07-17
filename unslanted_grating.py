# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 16:26:54 2025

@author: lyons
"""

import numpy as np
import pandas as pd
from math import pi
from numpy.linalg import inv

class unslanted_grating:
    
    def __init__(self, end_exp=1e2,total_time=100,lpmm=1e3,T0=50e-4,I0=5,xi=0.3,Dm=1.6e-7,
                 n_m=1.55,rhom=1.15,Dp=6.35e-10,rhop=1.3,n_p=1.56,n_q=1.64,Gamma=1,
                 Dz=1e-10,epsilon_mz=0,epsilon_pz=13,epsilon_qz=13,wt_pc=5e-2,rhoz=1.74,
                 n_z=1.366,b0=5.05,n_b=1.5,rhob=1.19,lambda_probe=633e-7,Delta_t=1/100,
                 Delta_x=1/20,output_time_step=1):
        
        self.total_time = total_time
        self.end_exp = end_exp
        self.lpmm = lpmm
        self.T0 = T0
        self.I0 = I0
        self.xi = xi
        self.Dm = Dm
        self.n_m = n_m
        self.rhom = rhom
        self.Dp = Dp
        self.rhop = rhop
        self.n_p = n_p
        self.n_q = n_q
        self.Gamma = Gamma
        self.Dz = Dz
        self.epsilon_mz = epsilon_mz
        self.epsilon_pz = epsilon_pz
        self.epsilon_qz = epsilon_qz
        self.wt_pc = wt_pc
        self.rhoz = rhoz
        self.n_z = n_z
        self.b0 = b0
        self.n_b = n_b
        self.rhob = rhob
        self.lambda_probe = lambda_probe
        self.Delta_x = Delta_x
        self.Delta_t = Delta_t
        self.output_time_step = output_time_step
        self.model = "coupled_cross_diffusion_v7"
        

    def coupled_cross_diffusion_v7(self):
        
        # 1.1 --- Define parameters
        iterations_per_second=1/self.Delta_t# Number of iterations each second
        Nx=int(1/self.Delta_x) + 1# Number of spatial points
        x = np.linspace(0,1,Nx)# Non-dimensional grating distance
        n_iterations = int(self.total_time*iterations_per_second) +1# Total number of iterations
        r=self.Delta_t/self.Delta_x/self.Delta_x# Ratio of finite time step to squared finite spatial step
        m0=1# Initial grams of monomer
        t0=1 #  Reference time [s]
        Lambda=1/10/self.lpmm # Grating period [cm]
        j_end_exp=self.end_exp/self.Delta_t # Iteration of exposure end
        z0 = self.wt_pc/(1 - self.wt_pc)*(m0 + self.b0)# Initial nanoparticle mass
        m1=np.ones(Nx); m2=m1
        p1=np.zeros(Nx); p2=p1
        q1=np.zeros(Nx); q2=q1
        z1=np.ones(Nx); z2=z1
        b1=self.b0*np.ones(Nx)
        
        Vb = b1/self.rhob # cm**3
        Vm = m1*m0/self.rhom # cm**3
        Vp = p1*m0/self.rhop # cm**3
        Vq = q1*m0/self.rhop # cm**3
        Vz = z1*z0/self.rhoz # cm**3
        Vtotal=Vb+Vm+Vp+Vq+Vz # cm**3
        
        phi_m = Vm/Vtotal
        phi_b = Vb/Vtotal
        phi_p = Vp/Vtotal
        phi_q = Vq/Vtotal
        phi_z = Vz/Vtotal
        
        Lorentz_Lorenz_RHS = phi_m*(self.n_m*self.n_m - 1)/(self.n_m*self.n_m + 2) + phi_b*(self.n_b*self.n_b - 1)/(self.n_b*self.n_b + 2) + phi_p*(self.n_p*self.n_p - 1)/(self.n_p*self.n_p + 2) + phi_q*(self.n_q*self.n_q - 1)/(self.n_q*self.n_q + 2) + phi_z*(self.n_z*self.n_z - 1)/(self.n_z*self.n_z + 2)
        
        n1=np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
        
        interior_points = list(range(1, Nx-1))
        
        times_4 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2 != 0]
        
        times_2 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2==0]
        
        N0=[self.Delta_x/3*(n1[0] + sum(2*n1[times_2]) + sum(4*n1[times_4]) + n1[Nx-1])]
        
        n1_cos=n1*np.cos(2*pi*x)
        
        N1=[2*self.Delta_x/3*(n1_cos[0] + sum(2*n1_cos[times_2]) + sum(4*n1_cos[times_4]) + n1_cos[Nx-1])]
            
        sq_diff=(n1 - N0[0] - N1[0]*np.cos(2*pi*x))**2
        
        d2 = [self.Delta_x/3*(sq_diff[0] + sum(2*sq_diff[times_2]) + sum(4*sq_diff[times_4]) + sq_diff[Nx-1])]
        
        # 1.2 --- Non-dimensional parameters
        alpha_m=self.Dm*t0/Lambda/Lambda# Monomer diffusion
        F0=0.1*self.I0**0.3
        beta=F0*t0# Monomer comsumption
        alpha_p=self.Dp*t0/Lambda/Lambda# Oligomer diffusion
        alpha_z=self.Dz*t0/Lambda/Lambda# Nondimensional self-nanoparticle diffusion
        gamma = m0*self.Gamma*t0# Immobilization
        if self.wt_pc==0:
            alpha_mz=0
            alpha_pz=0
            alpha_zm=0
            alpha_zp=0
            alpha_zq=0
        else:
            alpha_mz=z0*self.epsilon_mz*alpha_m
            alpha_pz=z0*self.epsilon_pz*alpha_p
            alpha_zm=self.epsilon_mz*alpha_z
            alpha_zp=self.epsilon_pz*alpha_z
            alpha_zq=self.epsilon_qz*alpha_z
        
        spatial_profiles_df=pd.DataFrame({"x": x, 
                         "monomer": m1,
                         "short_polymer": p1,
                         "immobile_polymer": q1,
                         'nanoparticles': z1,
                         'binder': b1,
                         'refractive_index': n1,
                         "time": np.zeros(Nx)})
        
        time_vals=np.arange(0, self.total_time+1, self.output_time_step)
        
        # 1.3 --- Calculate each time step via implicit finite difference method
        for j in range(1,n_iterations):
            # Phi=1 if illumination is on, 0 otherwise
            if ((j >= 0) and (j <= j_end_exp)):
                Phi=1
            else:
                Phi = 0
                
            # Illumination pattern
            f = 1 + np.exp(-self.xi*z0*z1)*np.cos(2*pi*x)
            MM2 = (2 + Phi*self.Delta_t*beta*f)*np.identity(Nx)
            MM1 = (2 - Phi*self.Delta_t*beta*f)*np.identity(Nx)
            PP2 = (2 + Phi*self.Delta_t*gamma*p1)*np.identity(Nx)
            PP1 = (2 - Phi*self.Delta_t*gamma*p1)*np.identity(Nx)
            PM2 = Phi*self.Delta_t*beta*f*np.identity(Nx)
            PM1 = Phi*self.Delta_t*beta*f*np.identity(Nx)
            QQ2=2*np.identity(Nx)
            QQ1=2*np.identity(Nx)
            QP2=(Phi*gamma*self.Delta_t*p1)*np.identity(Nx)
            QP1=(Phi*gamma*self.Delta_t*p1)*np.identity(Nx)
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
            
            Vb = b1/self.rhob # cm**3
            Vm = m1*m0/self.rhom # cm**3
            Vp = p1*m0/self.rhop # cm**3
            Vq = q1*m0/self.rhop # cm**3
            Vz = z1*z0/self.rhoz # cm**3
            Vtotal=Vb+Vm+Vp+Vq+Vz # cm**3
            phi_m = Vm/Vtotal
            phi_b = Vb/Vtotal
            phi_p = Vp/Vtotal
            phi_q = Vq/Vtotal
            phi_z = Vz/Vtotal
        
            Lorentz_Lorenz_RHS = phi_m*(self.n_m*self.n_m - 1)/(self.n_m*self.n_m + 2) + phi_b*(self.n_b*self.n_b - 1)/(self.n_b*self.n_b + 2) + phi_p*(self.n_p*self.n_p - 1)/(self.n_p*self.n_p + 2) + phi_q*(self.n_q*self.n_q - 1)/(self.n_q*self.n_q + 2) + phi_z*(self.n_z*self.n_z - 1)/(self.n_z*self.n_z + 2)
        
            n1=np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
        
            N0_new=self.Delta_x/3*(n1[0] + sum(2*n1[times_2]) + sum(4*n1[times_4]) + n1[Nx-1])
            
            n1_cos=n1*np.cos(2*pi*x)
        
            N1_new=2*self.Delta_x/3*(n1_cos[0] + sum(2*n1_cos[times_2]) + sum(4*n1_cos[times_4]) + n1_cos[Nx-1])
            
            sq_diff=(n1 - N0_new - N1_new*np.cos(2*pi*x))**2
        
            d2_new = [self.Delta_x/3*(sq_diff[0] + sum(2*sq_diff[times_2]) + sum(4*sq_diff[times_4]) + sq_diff[Nx-1])]
            
            if self.Delta_t*j in time_vals:
                
                new_df=pd.DataFrame({"x": x, 
                                 "monomer": m1,
                                 "short_polymer": p1,
                                 "immobile_polymer": q1,
                                 'nanoparticles': z1,
                                 'binder':b1,
                                 'refractive_index': n1,
                                 "time": np.zeros(Nx) + self.Delta_t*j})
                
                spatial_profiles_df=pd.concat([spatial_profiles_df,new_df])
                
                N0.append(N0_new)
                
                N1.append(N1_new)
                
                d2.append(d2_new)
                
        spatial_profiles_df=spatial_profiles_df.reset_index(drop=True)
        
        optical_properties_df=pd.DataFrame({'time': time_vals, 'N0': N0, 'Delta_n': [2*i for i in N1]})              
        
        return spatial_profiles_df, optical_properties_df