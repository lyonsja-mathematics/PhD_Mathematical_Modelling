# Holographic_Sensor_1D_Model_v1.R
# Author: Jack Lyons

# A function to run a numerical simulation of a mathematical model to describe operation of SRG based holographic sensor. This script differs from earlier versions in that instead of filling up the troughs with nanoparticles, a specified quantity of nanoparticles are coated upon a undoped holographic grating and will undergo redistribution governed by the same equation for hybrid holographic grating formation.

# 1 spatial dimension, 1 analyte

# 0.1 --- Import libraries
library(dplyr)
library(reshape2)
library(ggplot2)

# 0.2 --- Increase memory limit
memory.limit(size=1e6)

# 0.3 --- Set working directory
setwd("~/PhD Documents/Holographic Sensors/Modelling")

# 0.4 --- Import plot aesthetics
load("~/PhD Documents/plt.theme1.RData")
load("~/PhD Documents/plt.theme2.RData")

# 1.1 --- Define function inputs and default values.
sensor_simulation_v11 = function(
    
	sim_holo_grat=NULL,
	ns=1.33,
	na=2,
	nzs=1.46,
	nza=1.6,
	tau_c_s=60,# seconds
	tau_c_a=1,# seconds
	tau_e_s=60,# seconds
	tau_e_a=60,# seconds
	a0=0.1*0.1,# m0
	s0=0.1,# m0
	rhoa=9,# g/cm3
	rhos=1,# g/cm3
	Da=2.3e-5,# g/cm3
	Ds=2.3e-5,# cm2/s
	lambda_probe=633e-7,# cm
	exposure_time=180,# seconds
	output_time_step=5# seconds

	
){
	
	if(is.null(sim_holo_grat)){
		return(warning("Function needs a simulation of a holographic grating"))
	}
	
	# if(!(unique(sim_holo_grat$Model)=="slanted_grating_simulation_v16")){
	# 	return(warning("Theoretically modelled nanocomposite needs to be simulated with slanted_grating_simulation_v16()."))
	# }

	# 1.2 --- Define the holographic grating parameters
	sim_holo_grat %>% pull(Delta.x) %>% unique() -> Delta_x
	sim_holo_grat %>% pull(Delta.t) %>% unique() -> Delta_t
	r = Delta_t/Delta_x/Delta_x
	Nx=1/Delta_x + 1# Number of spatial points
	
	interior_points = 2:(Nx-1)
	times_4 = interior_points[interior_points %% 2 == 0]
	times_2 = interior_points[interior_points %% 2 != 0]
	
	df0=subset(sim_holo_grat, time==max(time))
	
	df0 %>% pull(T.t) %>% unique() -> T1
	
	df0 %>% pull(lpmm) %>% unique() -> lpmm
	
	df0 %>% pull(b0) %>% unique() -> b0
	
	Lambda=1/10/lpmm
	
	df0 %>% pull(nz) %>% unique() -> nze
	
	#LLRHS = phi_interior*((ns*ns - 1)/(ns*ns + 2)) + ((nze*nze - 1)/(nze*nze + 2))
	
	#nzs = sqrt((2*LLRHS + 1)/(1 - LLRHS))
	
	df0 %>% pull(wt.pc) %>% unique() -> wt_pc	
	
	t0=1
	
	m0=1
	
	sim_holo_grat %>% pull(z0) %>% unique() -> z0
	
	if(z0==0){z0=1e-6}
	if(tau_c_s==0){tau_c_s=1e-4}
	if(tau_e_s==0){tau_e_s=1e-4}
	
	gamma_s=1/z0/tau_c_s
	gamma_a=1/z0/tau_c_a
	omega_s=s0/z0/tau_e_s
	omega_a=a0/z0/tau_e_a

	df0 %>% arrange(x) %>% pull(n) -> n_before_exposure
	df0 %>% pull(nb) %>% unique -> nb
	df0 %>% pull(nm) %>% unique -> nm
	df0 %>% pull(np) %>% unique -> np
	df0 %>% pull(nq) %>% unique -> nq
	df0 %>% pull(nz) %>% unique -> nze
	
	rep(1, Nx) -> b1
	df0 %>% arrange(x) %>% pull(mass.m) -> m1
	df0 %>% arrange(x) %>% pull(mass.p) -> p1
	df0 %>% arrange(x) %>% pull(mass.q) -> q1
	df0 %>% arrange(x) %>% mutate(z1=ifelse(z0==0, mass.z, mass.z/z0)) %>% pull(z1) -> ze1
	rep(1, Nx) -> s1
	rep(1, Nx) -> a1
	rep(0, Nx) -> zs1
	rep(0, Nx) -> za1
	
	matrix(a1, ncol=1) -> a
	matrix(s1, ncol=1) -> s
	matrix(b1, ncol=1) -> b
	matrix(m1, ncol=1) -> m
	matrix(p1, ncol=1) -> p
	matrix(q1, ncol=1) -> q
	matrix(ze1, ncol=1) -> ze
	matrix(zs1, ncol=1) -> zs
	matrix(za1, ncol=1) -> za
	
	int_m_dx=Delta_x/3*(m1[1] + sum(2*m1[times_2]) + sum(4*m1[times_4]) + m1[Nx])
	int_p_dx=Delta_x/3*(p1[1] + sum(2*p1[times_2]) + sum(4*p1[times_4]) + p1[Nx])
	int_q_dx=Delta_x/3*(q1[1] + sum(2*q1[times_2]) + sum(4*q1[times_4]) + q1[Nx])
	
	sim_holo_grat %>% pull(rhom) %>% unique() -> rhom
	sim_holo_grat %>% pull(rhop) %>% unique() -> rhop
	sim_holo_grat %>% pull(rhob) %>% unique() -> rhob
	sim_holo_grat %>% pull(rhoz) %>% unique() -> rhoz
	
	v1=b0/rhob + int_m_dx/rhom + int_p_dx/rhop + int_q_dx/rhop + z0/rhoz
	v2=s0/rhos + a0/rhoa + v1
	T2=c(T1*v2/v1)
	
	alpha_s_x=Ds*t0/Lambda/Lambda
	alpha_a_x=Da*t0/Lambda/Lambda
	gamma_ss=gamma_s*z0*t0
	omega_ss=ifelse(s0==0, 0, omega_s*z0*t0/s0)
	gamma_sz=gamma_s*s0*t0
	omega_sz=omega_s*t0
	gamma_aa=gamma_a*z0*t0
	omega_aa=ifelse(a0==0, 0, omega_a*z0*t0/a0)
	gamma_az=gamma_a*a0*t0
	omega_az=omega_a*t0
	
	Nt = exposure_time/Delta_t+1# Total number of iterations
	
	seq(0, exposure_time, by=output_time_step) -> time_vals
	
	if(exposure_time < output_time_step){
		bind_iter="0"
	} else {
		time_vals/Delta_t + 1 -> bind_iter
		bind_iter[-1] -> bind_iter
	}
	
	# Simulation of holographic grating exposed to a loaded solvent
	for(j in 1:Nt){
			
			AA2 = diag(2 + 2*r*alpha_a_x + gamma_aa*Delta_t*as.numeric(ze1), Nx)
			AA1 = diag(2 - 2*r*alpha_a_x - gamma_aa*Delta_t*as.numeric(ze1), Nx)
			  
		  SS2 = diag(2 + 2*r*alpha_s_x + gamma_ss*Delta_t*as.numeric(ze1), Nx)
		  SS1 = diag(2 - 2*r*alpha_s_x - gamma_ss*Delta_t*as.numeric(ze1), Nx)
		  
		  ZEZE2 = diag(2 + gamma_sz*Delta_t*as.numeric(s1) + gamma_az*Delta_t*as.numeric(a1), Nx)
		  ZEZE1 = diag(2 - gamma_sz*Delta_t*as.numeric(s1) - gamma_az*Delta_t*as.numeric(a1), Nx)
		  
		  ZSZS2 = diag(2 + omega_sz*Delta_t, Nx)
		  ZSZS1 = diag(2 - omega_sz*Delta_t, Nx)
		  
		  ZAZA2 = diag(2 + omega_az*Delta_t, Nx)
		  ZAZA1 = diag(2 - omega_az*Delta_t, Nx)
		  
		  for(i in 1:Nx){
		    
		    i_minus_1=ifelse(i == 1, i + 1, i - 1)
	      i_plus_1=ifelse(i == Nx, i - 1, i + 1)
	          
		    SS2[i, i_minus_1] = SS2[i, i_minus_1] - r*alpha_s_x
				SS2[i, i_plus_1] = SS2[i, i_plus_1] - r*alpha_s_x
				
				SS1[i, i_minus_1] = SS1[i, i_minus_1] + r*alpha_s_x
				SS1[i, i_plus_1] = SS1[i, i_plus_1] + r*alpha_s_x
				
				AA2[i, i_minus_1] = AA2[i, i_minus_1] - r*alpha_a_x
				AA2[i, i_plus_1] = AA2[i, i_plus_1] - r*alpha_a_x
				
				AA1[i, i_minus_1] = AA1[i, i_minus_1] + r*alpha_a_x
				AA1[i, i_plus_1] = AA1[i, i_plus_1] + r*alpha_a_x
					    	
	    }
		  
		  s2 = solve(SS2) %*% ((SS1 %*% matrix(s1, ncol=1)) + matrix(2*Delta_t*omega_ss*zs1, ncol=1))
		  a2 = solve(AA2) %*% ((AA1 %*% matrix(a1, ncol=1)) + matrix(2*Delta_t*omega_aa*za1, ncol=1))
		  ze2 = solve(ZEZE2) %*% ((ZEZE1 %*% matrix(ze1, ncol=1)) + matrix(2*omega_sz*Delta_t*zs1, ncol=1) + matrix(2*omega_az*Delta_t*za1, ncol=1))
		  zs2 = solve(ZSZS2) %*% ((ZSZS1 %*% matrix(zs1, ncol=1)) + matrix(2*Delta_t*gamma_sz*s1*ze1, ncol=1))
		  za2 = solve(ZAZA2) %*% ((ZAZA1 %*% matrix(za1, ncol=1)) + matrix(2*Delta_t*gamma_az*a1*ze1, ncol=1))
		  
		  s2 -> s1
		  a2 -> a1
		  ze2 -> ze1
		  zs2 -> zs1
		  za2 -> za1
		  
		  int_s_dx=Delta_x/3*(s1[1] + sum(2*s1[times_2]) + sum(4*s1[times_4]) + s1[Nx])
		  int_a_dx=Delta_x/3*(a1[1] + sum(2*a1[times_2]) + sum(4*a1[times_4]) + a1[Nx])
		  
		  v2 = v1 + s0*int_s_dx/rhos + a0*int_a_dx/rhoa
		  
		  if(as.character(j) %in% as.character(bind_iter)){
		  	
			  	matrix(cbind(s, s1), nrow=Nx) -> s
		  		matrix(cbind(a, a1), nrow=Nx) -> a
			  	matrix(cbind(ze, ze1), nrow=Nx) -> ze
			  	matrix(cbind(zs, zs1), nrow=Nx) -> zs
			  	matrix(cbind(za, za1), nrow=Nx) -> za
			  	T2 = c(T2, T1*v2/v1)
			  	
		  }
		  
	}

	# Mass concentration data frame
	x=seq(0,1,Delta_x)
	time=seq(0,exposure_time,output_time_step)
	Nt=length(time)
	
	data.frame(x=rep(x, Nt), time=sort(rep(time, Nx)), s=s0*melt(s)$value, a=a0*melt(a)$value, ze=z0*melt(ze)$value, zs=z0*melt(zs)$value, za=z0*melt(za)$value, m=rep(m1, Nt), p=rep(p1, Nt), q=rep(q1, Nt), b=rep(b0, Nt), n_before_exposure=rep(n_before_exposure, Nt)) -> df2
	
	df2 %>% 
		
		mutate(
			
			rhoz=rhoz,
			rhos=rhos,
			rhoa=rhoa,
			nze=nze,
			nzs=nzs,
			nza=nza,
			ns=ns,
			na=na,
			
			Vm=m/rhom,
			Vp=p/rhop,
			Vq=q/rhop,
			Vb=b/rhob,
			Vze=ze/rhoz,
			Vzs=zs/rhoz,
			Vza=za/rhoz,
			Vs=s/rhos,
			Va=a/rhoa,
			
			V=Vb+Vm+Vp+Vq+Vze+Vzs+Vza+Vs+Va,
			
			phi_m=Vm/V,
			phi_p=Vp/V,
			phi_q=Vq/V,
			phi_b=Vb/V,
			phi_ze=Vze/V,
			phi_zs=Vzs/V,
			phi_za=Vza/V,
			phi_s=Vs/V,
			phi_a=Va/V,
			
			Lorentz_Lorenz_RHS = phi_m*(nm*nm - 1)/(nm*nm + 2) + phi_p*(np*np - 1)/(np*np + 2) + phi_q*(nq*nq - 1)/(nq*nq + 2) + phi_b*(nb*nb - 1)/(nb*nb + 2) + phi_ze*(nze*nze - 1)/(nze*nze + 2) + phi_zs*(nzs*nzs - 1)/(nzs*nzs + 2) + phi_za*(nza*nza - 1)/(nza*nza + 2) + phi_s*(ns*ns - 1)/(ns*ns + 2) + phi_a*(na*na - 1)/(na*na + 2),
			
			n_after_exposure=sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
		
		) %>% select(-Lorentz_Lorenz_RHS,-Vb,-Vm,-Vp,-Vq,-Vze,-Vzs,-Vs,-V,-phi_b,-phi_m,-phi_p,-phi_q,-phi_ze,-phi_zs,-phi_s,-phi_a,-phi_za,-Va,-Vza) %>% arrange(time,x) -> df3
	
	#Refractive index modulation data frame
	n = df3 %>% arrange(time,x) %>% pull(n_after_exposure) %>% matrix(nrow=Nx)
	
	N0=rep(0, Nt)
  N1=rep(0, Nt)
  n.tilde=matrix(0, nrow=nrow(n), ncol=ncol(n))
  diff=matrix(0, nrow=nrow(n), ncol=ncol(n))
  n2=matrix(0, nrow=nrow(n), ncol=ncol(n))
  d2=rep(0, Nt)
  int.n2=rep(0,Nt)
  
  for(k in 1:Nt){# For each time point ...
    
    interior.points = 2:(Nx-1)
    times.4 = interior.points[interior.points %% 2 == 0]
    times.2 = interior.points[interior.points %% 2 != 0]
    
    n0=n[, k]
    n1=n[, k]*cos(2*pi*x)
    
    N0[k]=Delta_x/3*(n0[1] + n0[Nx] + sum(2*n0[times.2]) + sum(4*n0[times.4]))
    N1[k]=2*Delta_x/3*(n1[1] + n1[Nx] + sum(2*n1[times.2]) + sum(4*n1[times.4]))
    
    n.tilde[, k] = N0[k] + N1[k]*cos(2*pi*x)
    diff[, k] = (n[, k] - n.tilde[, k])**2
    d2[k] = Delta_x/3*(diff[1, k] + diff[Nx, k] + sum(2*diff[times.2, k]) + sum(4*diff[times.4, k]))
    n2[, k]=n[, k]**2
    int.n2[k]=Delta_x/3*(n2[1, k] + sum(2*n2[times.2, k]) + sum(4*n2[times.4, k]) + n2[Nx, k])
    
  }# End for loop
  
  df3 %>% arrange(time) %>% mutate(Delta_n=0, distortion=0, Mean_RI=0, theta_B=0) -> df4
  
  for(i in x){
  	df4[df4$x==i, "Delta_n"]=2*N1
  	df4[df4$x==i, "distortion"]=d2/int.n2
  	df4[df4$x==i, "Mean_RI"]=N0
  	df4[df4$x==i, "T2"]=T2
  }
  
  df4 %>% mutate(
  	theta_B=asin(lambda_probe/2/Mean_RI/Lambda), 
  	Moharam.Young=lambda_probe*lambda_probe/Mean_RI/Delta_n/Lambda/Lambda/cos(0), 
  	Klein.Cook=2*pi*lambda_probe*T2/Mean_RI/Lambda/Lambda/cos(0)
  ) -> df4
  
  if(df4 %>% subset(time==0) %>% pull(Klein.Cook) %>% unique < 10){
	
		df4 %>% mutate(J0=pi*Delta_n*T2/lambda_probe/2, J1=J0) -> df4
	
		for(l in 1:100){
			
			df4$J1 = df4$J1 + ((-1)**l)/factorial(l)/factorial(l+1)*df4$J0**(2*l + 1)
		
		}
			
		df4 %>% mutate(
			
			eta=J1*J1, 
			Geometry="Planar"
			
		) %>% select(-J1,-J0) -> df5
		
	} else {
			
			df4 %>% mutate(
				
				Delta_phi_r=0,
				nu=pi*Delta_n*T2/lambda_probe/cos(theta_B), 
				xi=pi*T2*Delta_phi_r/Lambda, 
				eta=sin(sqrt(xi*xi + nu*nu))**2/sqrt(1 + (xi*xi)/(nu*nu)), 
				Geometry="Volume"
				
			) -> df5
			
	}
	
	df5 %>% subset(time==0) %>% pull(eta) %>% unique -> pre_exposure_eta
	df5 %>% subset(time==0) %>% pull(theta_B) %>% unique -> pre_exposure_Bragg_angle
	
	df5 %>% mutate(
		
		time=as.numeric(time),
		x=as.numeric(x),
		pre_exposure_eta,
		pre_exposure_Bragg_angle,
		lambda_r_t=2*Mean_RI*Lambda*sin(theta_B),
		normalized_eta=eta/pre_exposure_eta,
		Delta_theta_B=theta_B-pre_exposure_Bragg_angle,
		nze,nzs,ns,z0,s0,a0,b0,tau_c_s,tau_e_s,tau_c_a,tau_e_a,
		tau_x_s=Lambda*Lambda/Ds,tau_x_a=Lambda*Lambda/Da,
		gamma_s,omega_s,rhoz,rhos,Ds,wt_pc,
		T1, Delta_T=T2-T1,
		lpmm,exposure_time,output_time_step,
		model="sensor_simulation_v11"
	) %>% arrange(time,x) -> df6

	return(df6)
	
}
