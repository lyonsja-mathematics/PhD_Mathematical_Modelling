#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 14:01:40 2025

@author: lyonsja
"""

# Fixed spatial domain, radiation BC

# 0.1 --- Import libraries
from warnings import filterwarnings, warn
filterwarnings("ignore")
import numpy as np
import pandas as pd
from math import pi, factorial
from numpy.linalg import inv
from time import time as gettime

def simpsons_rule_1D(array_1D, step):
    
    intpts=list(range(1,len(array_1D)-1))
    times_4=[i for i in intpts if i%2!=0]
    times_2=[i for i in intpts if i%2==0]
    
    return step/3*(array_1D[0] + 4*sum(array_1D[times_4]) + 2*sum(array_1D[times_2]) + array_1D[len(array_1D)-1])

def sensor(

	grating#,
	n_s=1.33#,
	n_a=2#,
	n_zs=1.46#,
	n_za=1.6#,
	tau_c_s=60#,# seconds
	tau_c_a=1#,# seconds
	tau_e_s=60#,# seconds
	tau_e_a=60#,# seconds
	a0=0.1*0.1#,# m0
	s0=0.1#,# m0
	rhoa=9#,# g/cm3
	rhos=1#,# g/cm3
	Da=2.3e-5#,# cm/2
	Ds=2.3e-5#,# cm2/s
	lambda_probe=633e-7#,# cm
	exposure_time=180#,# seconds
	output_time_step=5):
    
    if "spatial_profile_DF" not in dir(grating):
    
        return warn("Sensor needs a theoretically modelled holographic grating")


	if grating.slant_angle==0:
        
        grating.slant_angle=1e-4
	
	# 1.2 --- Define the holographic grating parameters
	Delta_x = grating.Delta_x
	Delta_t = grating.Delta_t
	r = Delta_t/Delta_x/Delta_x
	Nx = int(1/Delta_x + 1)# Number of spatial points
	
	interior_points = list(range(1,Nx-1))
    times_4 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2 != 0]
    times_2 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2 == 0]
	
    final_time = max(grating.spatial_profile_DF.time)
	pre_exposure_spatial_profile_DF = grating.spatial_profile_DF[grating.spatial_profile_DF.time == final_time].reset_index(drop=True)
	
    lpmm = grating.lpmm
    
    b0 = grating.b0
    
    final_grating_properties = grating.shrinkage_DF[grating.shrinkage_DF.time == max(grating.shrinkage_DF.time)].reset_index(drop=True)
    final_grating_properties = final_grating_properties[['time','phi_t','Lambda_t','theta_B','Thickness','Mean_RI']]
    
    Lambda_0 = final_grating_properties.Lambda_t[0]
    
    phi_0 = final_grating_properties.phi_t[0]
    
    theta_B_0 = final_grating_properties.theta_B[0]
    
    T_0 = final_grating_properties.Thickness[0]
    
    Mean_RI_0 = final_grating_properties.Mean_RI[0]
    
	x_hat = Lambda_0/np.cos(phi_0)
	
	y_hat = Lambda_0/np.sin(phi_0)
	
	lambda_r_0 = 2*Mean_RI_0*Lambda_0*np.sin(phi_0)
    
    n_ze = grating.n_z
		
	t0=1
	
	m0=1
	
	if grating.z0==0:
        grating.z0=1e-6
	
    if tau_c_s==0:
        tau_c_s=1e-4
	
    if tau_e_s==0:
        tau_e_s=1e-4
        
	if tau_c_a==0:
        tau_c_a=1e-4
        
	if tau_e_a==0:
        tau_e_a=1e-4
	
	gamma_a=1/grating.z0/tau_c_a
	gamma_s=1/grating.z0/tau_c_s
	omega_a=a0/grating.z0/tau_e_a
	omega_s=s0/grating.z0/tau_e_s
    
    tau_y_s=T_0*T_0/Ds
    tau_y_a=T_0*T_0/Da
    tau_x_s=x_hat*x_hat/Ds
    tau_x_a=x_hat*x_hat/Da
    
	n_before_exposure = np.array(pre_exposure_spatial_profile_DF.refractive_index)
	n_b = grating.n_b
	n_m = grating.n_m
	n_p = grating.n_p
	n_q = grating.n_q
	
	free_surface = list(pre_exposure_spatial_profile_DF[pre_exposure_spatial_profile_DF.Y == 1].index)
	fixed_surface = list(pre_exposure_spatial_profile_DF[pre_exposure_spatial_profile_DF.Y == 0].index)
	
	b1 = np.ones(Nx*Nx)
	m1 = np.array(pre_exposure_spatial_profile_DF.monomer)
	p1 = np.array(pre_exposure_spatial_profile_DF.short_polymer)
	q1 = np.array(pre_exposure_spatial_profile_DF.immobile_polymer)
	
    if grating.z0 == 0:
        ze1 = np.zeros(Nx*Nx)
    
    else:
        ze1 = np.array(pre_exposure_spatial_profile_DF.nanoparticles)
    
    s1 = np.zeros(Nx*Nx); s1[free_surface] = 1
    a1 = np.zeros(Nx*Nx); a1[free_surface] = 1
	zs1 = np.zeros(Nx*Nx)
	za1 = np.zeros(Nx*Nx)
	
	a = a1
	s = s1
	b = b1
	m = m1
	p = p1
	q = q1
	ze = ze1
	zs = zs1
	za = za1
    
    spatial_profile_DF = pre_exposure_spatial_profile_DF.rename(columns={'nanoparticles':'nanoparticles_vacant'})
    spatial_profile_DF.time = 0
    spatial_profile_DF['nanoparticles_analyte'] = za1
    spatial_profile_DF['nanoparticles_solvent'] = zs1
    spatial_profile_DF['analyte'] = a1
    spatial_profile_DF['solvent'] = s1
    
    def calculate_RI(solvent, analyte, vacant_nanoparticles, 
                     solvent_occupied_nanoparticles, 
                     analyte_occupied_nanoparticles):
        
        mass_distributions = [grating.b0*b1,m1,p1,q1,
                              grating.z0*vacant_nanoparticles,
                              grating.z0*solvent_occupied_nanoparticles,
                              grating.z0*analyte_occupied_nanoparticles,
                              s0*solvent,
                              a0*analyte]
        
        densities = [grating.rhob,grating.rhom,grating.rhop,grating.rhop,
                 grating.rhoz,grating.rhoz,grating.rhoz,
                 rhos,rhoa]
        
        RIs = [grating.n_b, grating.n_m, grating.n_p, grating.n_q, n_ze, n_zs, n_za, n_s, n_a]
        
        Vtotal = 0
        
        for md, rho in zip(mass_distributions, densities):
            
            Vtotal = Vtotal + md/rho
            
        Lorentz_Lorenz_RHS = 0
            
        for md, rho, ri in zip(mass_distributions, densities, RIs):
            
            Lorentz_Lorenz_RHS = Lorentz_Lorenz_RHS + md/rho/Vtotal*(ri*ri - 1)/(ri*ri + 2)
            
        return np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
    
    n1 = calculate_RI(s1, a1, ze1, zs1, za1)
    spatial_profile_DF['refractive_index'] = n1
    
    x1 = np.array(pre_exposure_spatial_profile_DF.x)
    Y1 = np.array(pre_exposure_spatial_profile_DF.Y)
		
	matrix_m = m1.reshape(Nx,Nx).T
	matrix_p = p1.reshape(Nx,Nx).T
	matrix_q = q1.reshape(Nx,Nx).T
	
	int_m_dx_dy=Delta_x*Delta_x/4*(matrix_m[0,0] + matrix_m[Nx-1,0] + matrix_m[0,Nx-1] + matrix_m[Nx-1,Nx-1] + np.sum(2*matrix_m[0,range(1, Nx-1)]) + np.sum(2*matrix_m[range(1, Nx-1),0]) + np.sum(2*matrix_m[Nx-1,range(1, Nx-1)]) + np.sum(2*matrix_m[range(1, Nx-1),Nx-1]) + np.sum(4*matrix_m[range(1, Nx-1),range(1, Nx-1)]) )
	
	int_p_dx_dy=Delta_x*Delta_x/4*(matrix_p[0,0] + matrix_p[Nx-1,0] + matrix_p[0,Nx-1] + matrix_p[Nx-1,Nx-1] + np.sum(2*matrix_p[0,range(1, Nx-1)]) + np.sum(2*matrix_p[range(1, Nx-1),0]) + np.sum(2*matrix_p[Nx-1,range(1, Nx-1)]) + np.sum(2*matrix_p[range(1, Nx-1),Nx-1]) + np.sum(4*matrix_p[range(1, Nx-1),range(1, Nx-1)]) )
	
	int_q_dx_dy=Delta_x*Delta_x/4*(matrix_q[0,0] + matrix_q[Nx-1,0] + matrix_q[0,Nx-1] + matrix_q[Nx-1,Nx-1] + np.sum(2*matrix_q[0,range(1, Nx-1)]) + np.sum(2*matrix_q[range(1, Nx-1),0]) + np.sum(2*matrix_q[Nx-1,range(1, Nx-1)]) + np.sum(2*matrix_q[range(1, Nx-1),Nx-1]) + np.sum(4*matrix_q[range(1, Nx-1),range(1, Nx-1)]) )
    
    def calculate_volume(solvent, analyte):
        
        matrix_a = analyte.reshape(Nx,Nx).T
        
        matrix_s = solvent.reshape(Nx,Nx).T
        
        int_a_dx_dy=Delta_x*Delta_x/4*(matrix_a[0,0] + matrix_a[Nx-1,0] + matrix_a[0,Nx-1] + matrix_a[Nx-1,Nx-1] + np.sum(2*matrix_a[0,range(1, Nx-1)]) + np.sum(2*matrix_a[range(1, Nx-1),0]) + np.sum(2*matrix_a[Nx-1,range(1, Nx-1)]) + np.sum(2*matrix_a[range(1, Nx-1),Nx-1]) + np.sum(4*matrix_a[range(1, Nx-1),range(1, Nx-1)]) )
        
        int_s_dx_dy=Delta_x*Delta_x/4*(matrix_s[0,0] + matrix_s[Nx-1,0] + matrix_s[0,Nx-1] + matrix_s[Nx-1,Nx-1] + np.sum(2*matrix_s[0,range(1, Nx-1)]) + np.sum(2*matrix_s[range(1, Nx-1),0]) + np.sum(2*matrix_s[Nx-1,range(1, Nx-1)]) + np.sum(2*matrix_s[range(1, Nx-1),Nx-1]) + np.sum(4*matrix_s[range(1, Nx-1),range(1, Nx-1)]) )
        
        return b0/grating.rhob + int_m_dx_dy/grating.rhom + int_p_dx_dy/grating.rhop + int_q_dx_dy/grating.rhop + grating.z0/grating.rhoz + a0*int_a_dx_dy/rhoa + s0*int_s_dx_dy/rhos
    
    Volume0 = calculate_volume(s1, a1)
    
    #@simpsons_rule_1D
    def calculate_optical_properties(RI, Nx = grating.Nx, Delta_x = grating.Delta_x):
        
        n2=RI.reshape(Nx,Nx).T
        
        midpoint = int((Nx-1)/2)
        inside_ri = simpsons_rule_1D(n2[midpoint], step=Delta_x)
        
        mean_ri = np.mean(n2)
        
        N0=[(Delta_x/3*(n2[0,i] + np.sum(2*n2[times_2,i]) + np.sum(4*n2[times_4,i]) + n2[(Nx-1),i])) for i in range(Nx)]
        
        n2_cos=np.zeros(n2.shape)
        x = np.linspace(0,1,Nx)
        for i in range(Nx):
            n2_cos[:,i] = n2[:,i]*np.cos(2*pi*x)
            
        n2_sin=np.zeros(n2.shape)
        for i in range(Nx):
            n2_sin[:,i] = n2[:,i]*np.sin(2*pi*x)
        
        N1_a=np.array([2*Delta_x/3*(n2_cos[0,i] + np.sum(2*n2_cos[times_2,i]) + np.sum(4*n2_cos[times_4,i]) + n2_cos[(Nx-1),i]) for i in range(Nx)])
        
        N1_b=np.array([2*Delta_x/3*(n2_sin[0,i] + np.sum(2*n2_sin[times_2,i]) + np.sum(4*n2_sin[times_4,i]) + n2_sin[(Nx-1),i]) for i in range(Nx)])
        
        Delta_n=np.sqrt(N1_a*N1_a+N1_b*N1_b)
        
        n_tilde=np.ones((Nx,Nx))
        for i in range(Nx):
            n_tilde[:,i]=N0[i]*n_tilde[:,i] + N1_a[i]*np.cos(2*pi*x) + N1_b[i]*np.sin(2*pi*x)
        
        sq_diff=(n2 - n_tilde)**2
        
        d2=[(Delta_x/3*(sq_diff[0,i] + np.sum(2*sq_diff[times_2,i]) + np.sum(4*sq_diff[times_4,i]) + sq_diff[(Nx-1),i])) for i in range(Nx)]
        
        u1 = calculate_volume(s1, a1)/Volume0
        
        actual_shrinkage = 1-u1
      
        phi_1=np.arctan(np.tan(phi_0)/u1)
      
        Lambda_1=np.cos(phi_1)/np.cos(phi_0)*Lambda_0
        
        y_hat_1=Lambda_1/np.sin(phi_1)
      
        theta_B_1=np.arcsin(grating.lambda_probe/2/inside_ri/Lambda_1)-phi_1
      
        Delta_theta_B = theta_B_0-theta_B_1
        
        nu=pi*Delta_n*T_0*u1/lambda_probe/np.cos(theta_B_1)
        
        lambda_r=2*mean_ri*Lambda_1*np.sin(theta_B_1 + phi_1)
        
        if grating.Klein_Cook < 10:
            
            xi=0
            J0=nu/2
            J1=J0
      	
            for l in range(1,101):
                
                J1 = J1 + ((-1)**l)/factorial(l)/factorial(l+1)*J0**(2*l + 1)
                  
            eta=J1*J1
            
        else:
            
            xi=pi*T_0*u1*(phi_1 - phi_0)/Lambda_1
            
            eta = np.sin(np.sqrt(xi*xi + nu*nu))**2/np.sqrt(1 + (xi*xi)/(nu*nu))
        
        return mean_ri, Delta_n, theta_B_1, nu, lambda_r, eta
        
    pre_exposure_Mean_RI, pre_exposure_Delta_n, pre_exposure_theta_B, pre_exposure_nu, pre_exposure_lambda_r, pre_exposure_eta = calculate_optical_properties(n1)
    
	optical_properties_DF = pd.DataFrame({
        'time':np.zeros(Nx),
        'Y': np.linspace(0, 1, Nx),
        'Delta_n': pre_exposure_Delta_n,
        'nu': pre_exposure_nu,
        'eta': pre_exposure_eta})
    
    Bragg_Angle_Wavelength_DF = pd.DataFrame({
        'time': [0],
        'theta_B': pre_exposure_theta_B,
        'Mean_RI': pre_exposure_Mean_RI,
        'lambda_r': pre_exposure_lambda_r})
	
    # 3.1 --- Diffusion of target a
	alpha_s_x=Ds*t0/x_hat/x_hat
	alpha_a_x=Da*t0/x_hat/x_hat
	alpha_s_y=Ds*t0/T_0/T_0
	alpha_a_y=Da*t0/T_0/T_0
	gamma_ss=gamma_s*grating.z0*t0
    
    if s0==0:
        omega_ss = 0
        
    else:
        omega_ss = omega_s*grating.z0*t0/s0
	
    gamma_sz=gamma_s*s0*t0
	omega_sz=omega_s*t0
	gamma_aa=gamma_a*grating.z0*t0

    if a0==0:
        omega_aa = 0

    else:
        omega_aa = omega_a*grating.z0*t0/a0
	
    gamma_az=gamma_a*a0*t0
	omega_az=omega_a*t0
	
    n_iterations = exposure_time/Delta_t + 1# Total number of iterations
    
    time_vals = list(range(0, exposure_time+1, output_time_step))
	
    if exposure_time < output_time_step:
        bind_iter = 0
         
    else:
        bind_iter = [i/Delta_t for i in time_vals]
        bind_iter = bind_iter[1::]
	
	# Simulation of holographic grating exposed to a loaded solvent
	for j in range(0, n_iterations):
        
        AA2 = (2 + 2*r*alpha_a_x + 2*r*alpha_a_y + gamma_aa*Delta_t*ze1)*np.identity(Nx*Nx)
        AA1 = (2 - 2*r*alpha_a_x - 2*r*alpha_a_y - gamma_aa*Delta_t*ze1)*np.identity(Nx*Nx)
        
        SS2 = (2 + 2*r*alpha_s_x + 2*r*alpha_s_y + gamma_ss*Delta_t*ze1)*np.identity(Nx*Nx)
        SS1 = (2 - 2*r*alpha_s_x - 2*r*alpha_s_y - gamma_ss*Delta_t*ze1)*np.identity(Nx*Nx)
        
        ZEZE2 = (2 + gamma_az*Delta_t*a1 + gamma_sz*Delta_t*s1)*np.identity(Nx*Nx)
        ZEZE1 = (2 - gamma_az*Delta_t*a1 - gamma_sz*Delta_t*s1)*np.identity(Nx*Nx)
        
        ZAZA2 = (2 + omega_az*Delta_t)*np.identity(Nx*Nx)
        ZAZA1 = (2 - omega_az*Delta_t)*np.identity(Nx*Nx)
        
        ZSZS2 = (2 + omega_sz*Delta_t)*np.identity(Nx*Nx)
        ZSZS1 = (2 - omega_sz*Delta_t)*np.identity(Nx*Nx)
        
        BB2 = 2*np.identity(Nx*Nx)
  		BB1 = 2*np.identity(Nx*Nx)
  		
  		MM2 = 2*np.identity(Nx*Nx)
  		MM1 = 2*np.identity(Nx*Nx)
  		
  		PP2 = 2*np.identity(Nx*Nx)
  		PP1 = 2*np.identity(Nx*Nx)
  		
  		QQ2 = 2*np.identity(Nx*Nx)
  		QQ1 = 2*np.identity(Nx*Nx)
        
        for i in range(0, Nx*Nx):
            
            if i in list(pre_exposure_spatial_profile_DF[pre_exposure_spatial_profile_DF.x == 0].index):
                i_minus_1 = i + Nx - 2
            else:
                i_minus_1 = i - 1
                
            if i in list(pre_exposure_spatial_profile_DF[pre_exposure_spatial_profile_DF.x == 1].index):
                i_plus_1 = i - Nx + 2
            else:
                i_plus_1 =  i + 1
    
            if i in list(pre_exposure_spatial_profile_DF[pre_exposure_spatial_profile_DF.Y == 0].index):
                j_minus_1 = i + Nx
            else:
                j_minus_1 = i - Nx
                
            if i in list(pre_exposure_spatial_profile_DF[pre_exposure_spatial_profile_DF.Y == 1].index):
                j_plus_1 = i - Nx
            else:
                j_plus_1 = i + Nx
	    
            AA2[i, i_minus_1] = AA2[i, i_minus_1] - r*alpha_a_x
			AA2[i, j_minus_1] = AA2[i, j_minus_1] - r*alpha_a_y
			AA2[i, i_plus_1] = AA2[i, i_plus_1] - r*alpha_a_x
			AA2[i, j_plus_1] = AA2[i, j_plus_1] - r*alpha_a_y

			SS2[i, i_minus_1] = SS2[i, i_minus_1] - r*alpha_s_x
			SS2[i, j_minus_1] = SS2[i, j_minus_1] - r*alpha_s_y
			SS2[i, i_plus_1] = SS2[i, i_plus_1] - r*alpha_s_x
			SS2[i, j_plus_1] = SS2[i, j_plus_1] - r*alpha_s_y
			
			AA1[i, i_minus_1] = AA1[i, i_minus_1] + r*alpha_a_x
			AA1[i, j_minus_1] = AA1[i, j_minus_1] + r*alpha_a_y
			AA1[i, i_plus_1] = AA1[i, i_plus_1] + r*alpha_a_x
			AA1[i, j_plus_1] = AA1[i, j_plus_1] + r*alpha_a_y

			SS1[i, i_minus_1] = SS1[i, i_minus_1] + r*alpha_s_x
			SS1[i, j_minus_1] = SS1[i, j_minus_1] + r*alpha_s_y
			SS1[i, i_plus_1] = SS1[i, i_plus_1] + r*alpha_s_x
			SS1[i, j_plus_1] = SS1[i, j_plus_1] + r*alpha_s_y
				    	
		
        a2 = np.matmul(inv(AA2), np.matmul(AA1, a1) + 2*Delta_t*omega_aa*za1)
        
        s2 = np.matmul(inv(SS2), np.matmul(SS1, s1) + 2*Delta_t*omega_ss*zs1)
        
        ze2 = np.matmul(inv(ZEZE2), np.matmul(ZEZE1, ze1) + 2*omega_az*Delta_t*za1 + 2*omega_sz*Delta_t*zs1)
        
        za2 = np.matmul(inv(ZAZA2), np.matmul(ZAZA1, za1) + 2*Delta_t*gamma_az*a1*ze1)
        
        zs2 = np.matmul(inv(ZSZS2), np.matmul(ZSZS1, zs1) + 2*Delta_t*gamma_sz*s1*ze1)
        
        b2 = np.matmul(inv(BB2), np.matmul(BB1, b1))
        
        m2 = np.matmul(inv(MM2), np.matmul(MM1, m1))
        
        p2 = np.matmul(inv(PP2), np.matmul(PP1, p1))
        
        q2 = np.matmul(inv(QQ2), np.matmul(QQ1, q1))
        	  
	  a2[free_surface]=1
	  s2[free_surface]=1
	  
	  a1 = a2
	  s1 = s2
	  ze1 = ze2
	  za1 = za2
	  zs1 = zs2
	  b1 = b2
	  m1 = m2
	  p1 = p2
	  q1 = q2
      
      n1 = calculate_RI(s1, a1, ze1, zs1, za1)
      
      new_Mean_RI, new_Delta_n, new_theta_B, new_nu, new_lambda_r, new_eta = calculate_optical_properties(n1)
      
      if j*Delta_t in time_vals:
          
          spatial_profile_DF = pd.concat([spatial_profile_DF, 
                             pd.DataFrame({
                                 'x':x, 'Y': Y, 
                                 'time': np.zeros(Nx) + j*Delta_t,
                                 'binder': b1,
                                 'monomer': m1,
                                 'short_polymer': p1,
                                 'immobile_polymer': q1,
                                 'refractive_index': n1,
                                 'nanoparticles_vacant': ze1,
                                 'nanoparticles_solvent': zs1,
                                 'nanoparticles_analyte': za1})])
          
          optical_properties_DF = pd.concat([optical_properties_DF,
                                             pd.DataFrame({
                                                 'time':np.zeros(Nx) + j*Delta_t,
                                                 'Y': np.linspace(0, 1, Nx),
                                                 'Delta_n': new_Delta_n,
                                                 'nu': new_nu,
                                                 'eta': new_eta})])
          
          Bragg_Angle_Wavelength_DF = pd.concat([Bragg_Angle_Wavelength_DF,
                                                 pd.DataFrame({
                                                     'time': [j*Delta_t],
                                                     'theta_B': new_theta_B,
                                                     'Mean_RI': new_Mean_RI,
                                                     'lambda_r': new_lambda_r})
                                                 
    return spatial_profile_DF, optical_properties_DF, Bragg_Angle_Wavelength_DF
