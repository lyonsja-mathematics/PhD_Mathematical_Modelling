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

# 1.1 --- Define function inputs and default values.
#def sensor_2D_model_v2(

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
	output_time_step=5# seconds
	
#){

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
    
    theta_B_0 = sdf0.theta_B[0]
    
    T_0 = final_grating_properties.Thickness[0]
    
    Mean_RI_0 = final_grating_properties.Mean_RI[0]
    
	x_hat = Lambda_1/np.cos(phi_0)
	
	y_hat = Lambda_1/np.sin(phi_0)
	
	lambda_r_0 = 2*Mean_RI_0*Lambda_0*np.sin(phi_0)
    
    n_ze = grating.n_z
	
	#LLRHS = phi_interior*((ns*ns - 1)/(ns*ns + 2)) + ((nze*nze - 1)/(nze*nze + 2))
	
	#nzs = sqrt((2*LLRHS + 1)/(1 - LLRHS))
	
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
        
        return mean_ri, Delta_n, theta_B_1, nu, lambda_r
        
    pre_exposure_Mean_RI, pre_exposure_Delta_n, pre_exposure_theta_B, pre_exposure_nu, pre_exposure_lambda_r = calculate_optical_properties(n1)
    
	optical_properties_DF = pd.DataFrame({
        'time':np.zeros(Nx),
        'Y': np.linspace(0, 1, Nx),
        'Delta_n': pre_exposure_Delta_n,
        'nu': pre_exposure_nu})
    
    Bragg_Angle_Wavelength_DF = pd.DataFrame({
        'time': [0],
        'theta_B': pre_exposure_theta_B,
        'Mean_RI': pre_exposure_Mean_RI,
        'lambda_r': pre_exposure_lambda_r})
	
    # 3.1 --- Diffusion of target a
	alpha_s_x=Ds*t0/x_hat/x_hat
	alpha_a_x=Da*t0/x_hat/x_hat
	alpha_s_y=Ds*t0/T_1/T_1
	alpha_a_y=Da*t0/T_1/T_1
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
      
      new_Mean_RI, new_Delta_n, new_theta_B, new_nu, new_lambda_r = calculate_optical_properties(n1)
      
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
                                                 'Delta_n': pre_exposure_Delta_n,
                                                 'nu': pre_exposure_nu})])
          
          Bragg_Angle_Wavelength_DF = pd.concat([Bragg_Angle_Wavelength_DF,
                                                 pd.DataFrame({
                                                     'time': [j*Delta_t],
                                                     'theta_B': new_theta_B,
                                                     'Mean_RI': new_Mean_RI,
                                                     'lambda_r': new_lambda_r})
      	
          
          
          
          
          
          
          
          
	  	
	  	matrix(cbind(a, a1), nrow=Nx*Nx) -> a
	  	matrix(cbind(s, s1), nrow=Nx*Nx) -> s
	  	matrix(cbind(b, b1), nrow=Nx*Nx) -> b
	  	matrix(cbind(m, m1), nrow=Nx*Nx) -> m
	  	matrix(cbind(p, p1), nrow=Nx*Nx) -> p
	  	matrix(cbind(q, q1), nrow=Nx*Nx) -> q
	  	matrix(cbind(ze, ze1), nrow=Nx*Nx) -> ze
	  	matrix(cbind(za, za1), nrow=Nx*Nx) -> za
	  	matrix(cbind(zs, zs1), nrow=Nx*Nx) -> zs
	  	
  	}
  
	}

	# Mass concentration data frame
	x=seq(0,1,Delta_x)
	time=seq(0,exposure_time,output_time_step)
	Nt=length(time)
	
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),a) -> df_a
	names(df_a)=c("x", "Y", paste0("t",time_vals[1:ncol(a)]))# Rename columns
	df_a %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, a=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""), time=time), a=a0*a) %>% arrange(time,Y,x) -> df_a
	
	# Solvent mass data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),s) -> df_s
	names(df_s)=c("x", "Y", paste0("t",time_vals[1:ncol(s)]))# Rename columns
	df_s %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, s=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=time), s=s0*s) %>% arrange(time,Y,x) -> df_s
	
	# Empty nanozeolite mass data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),ze) -> df_ze
	names(df_ze)=c("x", "Y", paste0("t",time_vals[1:ncol(ze)]))# Rename columns
	df_ze %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, ze=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=time), ze=z0*ze) %>% arrange(time,Y,x) -> df_ze
	
	# Analyte-filled nanozeolite mass data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),za) -> df_za
	names(df_za)=c("x", "Y", paste0("t",time_vals[1:ncol(za)]))# Rename columns
	df_za %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, za=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=time), za=z0*za) %>% arrange(time,Y,x) -> df_za
	
	# Solvent-filled nanozeolite mass data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),zs) -> df_zs
	names(df_zs)=c("x", "Y", paste0("t",time_vals[1:ncol(zs)]))# Rename columns
	df_zs %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, zs=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=time), zs=z0*zs) %>% arrange(time,Y,x) -> df_zs
	
	# Binder data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),b) -> df_b
	names(df_b)=c("x", "Y", paste0("t",time_vals[1:ncol(b)]))# Rename columns
	df_b %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, b=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=time),b=b0*b) %>% arrange(time,Y,x) -> df_b
	
	# Monomer data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),m) -> df_m
	names(df_m)=c("x", "Y", paste0("t",time_vals[1:ncol(m)]))# Rename columns
	df_m %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, m=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=time),m=m) %>% arrange(time,Y,x) -> df_m
	
	# Short polymer data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),p) -> df_p
	names(df_p)=c("x", "Y", paste0("t",time_vals[1:ncol(p)]))# Rename columns
	df_p %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, p=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=time),p=p) %>% arrange(time,Y,x) -> df_p
	
	# Cross-linked polymer data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),q) -> df_q
	names(df_q)=c("x", "Y", paste0("t",time_vals[1:ncol(q)]))# Rename columns
	df_q %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, q=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=time),q=q) %>% arrange(time,Y,x) -> df_q
	
	# 3.1 --- Refractive index matrix
  Vb = b*b0/rhob # cm**3
  Vm = m*m0/rhom # cm**3
  Vp = p*m0/rhop # cm**3
  Vq = q*m0/rhop # cm**3
  Vze = ze*grating.z0/rhoz # cm**3
  Vza = za*grating.z0/rhoz # cm**3
  Vzs = zs*grating.z0/rhoz # cm**3
  Va = a*a0/rhoa # cm**3
  Vs = s*s0/rhos # cm**3
  Vtotal=Vb+Vm+Vp+Vq+Vze+Vza+Vzs+Va+Vs# cm**3
  
  phi.b=Vb/Vtotal
  phi.m=Vm/Vtotal
  phi.p=Vp/Vtotal
  phi.q=Vq/Vtotal
  phi.ze=Vze/Vtotal
  phi.zs=Vzs/Vtotal
  phi.za=Vza/Vtotal
  phi.a=Va/Vtotal
  phi.s=Vs/Vtotal
  
  Lorentz.Lorenz.RHS = phi.b*(nb*nb - 1)/(nb*nb + 2) + phi.m*(nm*nm - 1)/(nm*nm + 2) + phi.p*(np*np - 1)/(np*np + 2) + phi.q*(nq*nq - 1)/(nq*nq + 2) + phi.ze*(nze*nze - 1)/(nze*nze + 2) + phi.zs*(nzs*nzs - 1)/(nzs*nzs + 2) + phi.za*(nza*nza - 1)/(nza*nza + 2) + phi.a*(na*na - 1)/(na*na + 2) + phi.s*(ns*ns - 1)/(ns*ns + 2)
  
  n=sqrt((2*Lorentz.Lorenz.RHS + 1)/(1 - Lorentz.Lorenz.RHS))
  
  data.frame(x=rep(seq(0,1,Delta_x), Nx), Y=sort(rep(seq(0,1,Delta_x), Nx)), n) -> df_n
  names(df_n)=c("x","Y", paste0("t",time_vals[1:ncol(n)]))
  df_n %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, n_after_exposure=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=time)) %>% arrange(time,Y,x) -> df_n# Additional columns containing info on input
  
  cbind(df_b,m=df_m$m,p=df_p$p, q=df_q$q, ze=df_ze$ze, za=df_za$za, zs=df_zs$zs, a=df_a$a, s=df_s$s, n_after_exposure=df_n$n_after_exposure) -> df1
	
  Nt = length(time)
  N0=matrix(0, nrow=Nt, ncol=Nx)
  N1a=matrix(0, nrow=Nt, ncol=Nx)
  N1b=matrix(0, nrow=Nt, ncol=Nx)
  d2=matrix(0, nrow=Nt, ncol=Nx)
  int_n2=matrix(0, nrow=Nt, ncol=Nx)
  Mean_RI=rep(0,Nt)
  
  for(k in 1:Nt){

matrix(n[, k], nrow=Nx) -> n1
  	
  	n2=n1**2
  	
  	Mean_RI[k]=Delta_x*Delta_x/4*(n1[0,0] + n1[Nx-1,0] + n1[0,Nx-1] + n1[Nx-1,Nx-1] + np.sum(2*n1[0,range(1, Nx-1)]) + np.sum(2*n1[range(1, Nx-1),0]) + np.sum(2*n1[Nx-1,range(1, Nx-1)]) + np.sum(2*n1[range(1, Nx-1),Nx-1]) + np.sum(4*n1[range(1, Nx-1),range(1, Nx-1)]) )

for(i in 1:Nx){
  
  n1[, i] -> n1i
  N0[k, i]=Delta_x/3*(n1i[1] + np.sum(2*n1i[times_2]) + np.sum(4*n1i[times_4]) + n1i[Nx])
  
  n1[, i]*cos(2*pi*x) -> n1cos
  N1a[k, i]=2*Delta_x/3*(n1cos[1] + np.sum(2*n1cos[times_2]) + np.sum(4*n1cos[times_4]) + n1cos[Nx])
  
  n1[, i]*sin(2*pi*x) -> n1sin
  N1b[k, i]=2*Delta_x/3*(n1sin[1] + np.sum(2*n1sin[times_2]) + np.sum(4*n1sin[times_4]) + n1sin[Nx])
  
}


n.Fourier=matrix(0, nrow=Nx, ncol=Nx)

for(i in 1:Nx){
  
  n.Fourier[,i] = N0[k, i] + N1a[k, i]*cos(2*pi*x) + N1b[k, i]*sin(2*pi*x)
  
}

diff = abs(n1 - n.Fourier)**2

for(i in 1:Nx){
	d2[k,i] = Delta_x/3*(diff[1, i] + np.sum(2*diff[times_2, i]) + np.sum(4*diff[times_4, i]) + diff[Nx, i])
	int_n2[k,i] = Delta_x/3*(n2[1, i] + np.sum(2*n2[times_2, i]) + np.sum(4*n2[times_4, i]) + n2[Nx, i])
}

  }# End for loop
  
  delta=d2/int_n2
  Delta_n=2*sqrt(N1a*N1a + N1b*N1b)
  
  data.frame(
time=rep(time, Nx),
Y=sort(rep(seq(0,1,Delta_x), Nt)),
Delta_n=melt(Delta_n)$value,
delta=melt(delta)$value
  ) -> df2
  
  df1 %>% mutate(
  	
  	Delta_n=0, 
  	delta=0, 
  	Mean_RI=0
  	
  ) %>% arrange(time,Y,x) -> df1
  
  for(i in time){
  	for(j in seq(0,1,Delta_x)){
  		df1[df1$time==i & df1$Y==j, "Delta_n"] = df2[df2$time==i & df2$Y==j, "Delta_n"]
  		df1[df1$time==i & df1$Y==j, "delta"] = df2[df2$time==i & df2$Y==j, "delta"]
  	}
  }
  
  for(i in 1:length(time)){
  	df1[df1$time==time[i], "Mean_RI"]=Mean_RI[i]
  }
  
  df1 %>% mutate(
  	theta_B = asin(lambda_probe/2/Mean_RI/Lambda_1) - phi_0,
  	Delta_theta_B = theta_B1 - theta_B,
  	Delta_phi_r = 0
  ) -> df1
  
  if df0 %>% subset(Y==0) %>% distinct(Klein.Cook) %>% pull(Klein.Cook) < 10){
	
		df1 %>% mutate(nu=pi*Delta_n*T_1/lambda_probe, xi=0, J0=nu/2, J1=J0) -> df1
	
		for(l in 1:100){
			
			df1$J1 = df1$J1 + ((-1)**l)/factorial(l)/factorial(l+1)*df1$J0**(2*l + 1)
		
		}
			
		df1 %>% mutate(
			
			eta=J1*J1, 
			Geometry="Planar"
			
		) %>% select(-J1,-J0) -> df1
		
	} else {
			
			df1 %>% mutate(
				
				nu=pi*Delta_n*T_1/lambda_probe/cos(theta_B), 
				xi=pi*T_1*Delta_phi_r/Lambda_1, 
				eta=sin(sqrt(xi*xi + nu*nu))**2/sqrt(1 + (xi*xi)/(nu*nu)), 
				Geometry="Volume"
				
			) -> df1
			
	}
	
	df1 %>% subset(time==0) %>% rename(pre_exposure_eta=eta) %>% distinct(Y,pre_exposure_eta) -> df4
	
	df1 %>% mutate(pre_exposure_eta=0) -> df1
	
	for(i in seq(0,1,Delta_x)){
		df1[df1$Y==i, "pre_exposure_eta"] = df4[df4$Y==i,"pre_exposure_eta"]
	}
	
	df1 %>% mutate(
		
		normalized_eta=eta/pre_exposure_eta,
		lambda_r_t=2*Mean_RI*Lambda_1*sin(theta_B),
		nze=nze,
		nza=nza,
		nzs=nzs,
		na=na,
		ns=ns,
		z0=z0,
		a0=a0,
		s0=s0,
		b0=b0,
		tau_c_s=tau_c_s,
		tau_c_a=tau_c_a,
		tau_e_s=tau_e_s,
		tau_e_a=tau_e_a,
		tau_y_s=T_1*T_1/Ds,
		tau_y_a=T_1*T_1/Da,
		tau_x_s=x_hat*x_hat/Ds,
		tau_x_a=x_hat*x_hat/Da,
		gamma_a=gamma_a,
		gamma_s=gamma_s,
		omega_a=omega_a,
		omega_s=omega_s,
		rhoz=rhoz,
		rhoa=rhoa,
		rhos=rhos,
		Ds=Ds,
		Da=Da,
		wt_pc=wt_pc,
		lambda_probe=lambda_probe,
		T_1=T_1,
		pre_exposure_Bragg_angle=theta_B1,
		lpmm=lpmm,
		exposure_time=exposure_time,
		output_time_step=output_time_step,
		model="sensor_2D_model_v2"
	) %>% arrange(time,Y,x) -> df1
	
	return(df1)
	
}
