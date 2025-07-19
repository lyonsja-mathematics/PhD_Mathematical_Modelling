# sensor_2D_model_v2.R
# Author: Jack Lyons

# Fixed spatial domain, radiation BC

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
sensor_2D_model_v2 = function(
    
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
	Da=2.3e-5,# cm/2
	Ds=2.3e-5,# cm2/s
	lambda_probe=633e-7,# cm
	exposure_time=180,# seconds
	output_time_step=5# seconds
	
){
	
	if(is.null(sim_holo_grat)){
		return(
			warning(
				"Function needs a simulation of a holographic grating"
			)
		)
	}
	
	# if(!(unique(sim_holo_grat$Model)=="slanted_grating_simulation_v16")){
	# 	return(warning("Theoretically modelled nanocomposite needs to be simulated with slanted_grating_simulation_v16()."))
	# }

	if(unique(sim_holo_grat$slant.angle)==0){sim_holo_grat$slant.angle=1e-4}
	
	# 1.2 --- Define the holographic grating parameters
	sim_holo_grat %>% pull(Delta.x) %>% unique() -> Delta_x
	sim_holo_grat %>% pull(Delta.t) %>% unique() -> Delta_t
	r = Delta_t/Delta_x/Delta_x
	Nx=1/Delta_x + 1# Number of spatial points
	
	interior_points = 2:(Nx-1)
	times_4 = interior_points[interior_points %% 2 == 0]
	times_2 = interior_points[interior_points %% 2 != 0]
	
	df0=subset(sim_holo_grat, time==max(time))
	
	df0 %>% pull(lpmm) %>% unique() -> lpmm
	
	df0 %>% pull(b0) %>% unique() -> b0
	
	df0 %>% mutate(x_hat=Lambda.t/cos(phi.r.t)) %>% pull(x_hat) %>% unique() -> x_hat
	
	df0 %>% mutate(y_hat=Lambda.t/sin(phi.r.t)) %>% pull(y_hat) %>% unique() -> y_hat_1
	
	df0 %>% pull(T.t) %>% unique() -> T_1
	
	df0 %>% pull(Lambda.t) %>% unique() -> Lambda_1
	
	df0 %>% pull(phi.r.t) %>% unique() -> phi_r_1
	
	df0 %>% pull(theta_B) %>% unique() -> theta_B1
	
	df0 %>% mutate(lambda_r_1=2*Mean.RI*Lambda.t*sin(phi.r.t)) %>% pull(lambda_r_1) %>% unique() -> lambda_r_1
	
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
	if(tau_c_a==0){tau_c_a=1e-4}
	if(tau_e_a==0){tau_e_a=1e-4}
	
	gamma_a=1/z0/tau_c_a
	gamma_s=1/z0/tau_c_s
	omega_a=a0/z0/tau_e_a
	omega_s=s0/z0/tau_e_s
	
	df0 %>% arrange(x) %>% pull(n) -> n_before_exposure
	df0 %>% pull(nb) %>% unique -> nb
	df0 %>% pull(nm) %>% unique -> nm
	df0 %>% pull(np) %>% unique -> np
	df0 %>% pull(nq) %>% unique -> nq
	df0 %>% pull(nz) %>% unique -> nze
	free_surface=(Nx*Nx - Nx + 1):(Nx*Nx)
	fixed_surface=1:Nx
	
	rep(1, Nx*Nx) -> b1
	df0 %>% arrange(Y,x) %>% pull(m) -> m1
	df0 %>% arrange(Y,x) %>% pull(p) -> p1
	df0 %>% arrange(Y,x) %>% pull(q) -> q1
	df0 %>% arrange(Y,x) %>% mutate(z1=ifelse(z0==0, z, z/z0)) %>% pull(z1) -> ze1
	rep(0, Nx*Nx) -> s1; s1[free_surface]=1
	rep(0, Nx*Nx) -> a1; a1[free_surface]=1
	rep(0, Nx*Nx) -> zs1
	rep(0, Nx*Nx) -> za1
	
	matrix(a1, ncol=1) -> a
	matrix(s1, ncol=1) -> s
	matrix(b1, ncol=1) -> b
	matrix(m1, ncol=1) -> m
	matrix(p1, ncol=1) -> p
	matrix(q1, ncol=1) -> q
	matrix(ze1, ncol=1) -> ze
	matrix(zs1, ncol=1) -> zs
	matrix(za1, ncol=1) -> za
		
	m1 %>% matrix(nrow=Nx) -> matrix.m
	p1 %>% matrix(nrow=Nx) -> matrix.p
	q1 %>% matrix(nrow=Nx) -> matrix.q
	a1 %>% matrix(nrow=Nx) -> matrix.a
	s1 %>% matrix(nrow=Nx) -> matrix.s
	
	int_m_dx_dy=Delta_x*Delta_x/4*(matrix.m[1,1] + matrix.m[Nx,1] + matrix.m[1,Nx] + matrix.m[Nx,Nx] + sum(2*matrix.m[1,(2:(Nx-1))]) + sum(2*matrix.m[(2:(Nx-1)),1]) + sum(2*matrix.m[Nx,(2:(Nx-1))]) + sum(2*matrix.m[(2:(Nx-1)),Nx]) + sum(4*matrix.m[(2:(Nx-1)),(2:(Nx-1))]) )
	
	int_p_dx_dy=Delta_x*Delta_x/4*(matrix.p[1,1] + matrix.p[Nx,1] + matrix.p[1,Nx] + matrix.p[Nx,Nx] + sum(2*matrix.p[1,(2:(Nx-1))]) + sum(2*matrix.p[(2:(Nx-1)),1]) + sum(2*matrix.p[Nx,(2:(Nx-1))]) + sum(2*matrix.p[(2:(Nx-1)),Nx]) + sum(4*matrix.p[(2:(Nx-1)),(2:(Nx-1))]) )
	
	int_q_dx_dy=Delta_x*Delta_x/4*(matrix.q[1,1] + matrix.q[Nx,1] + matrix.q[1,Nx] + matrix.q[Nx,Nx] + sum(2*matrix.q[1,(2:(Nx-1))]) + sum(2*matrix.q[(2:(Nx-1)),1]) + sum(2*matrix.q[Nx,(2:(Nx-1))]) + sum(2*matrix.q[(2:(Nx-1)),Nx]) + sum(4*matrix.q[(2:(Nx-1)),(2:(Nx-1))]) )
	
	int_a_dx_dy=Delta_x*Delta_x/4*(matrix.a[1,1] + matrix.a[Nx,1] + matrix.a[1,Nx] + matrix.a[Nx,Nx] + sum(2*matrix.a[1,(2:(Nx-1))]) + sum(2*matrix.a[(2:(Nx-1)),1]) + sum(2*matrix.a[Nx,(2:(Nx-1))]) + sum(2*matrix.a[(2:(Nx-1)),Nx]) + sum(4*matrix.a[(2:(Nx-1)),(2:(Nx-1))]) )
	
	int_s_dx_dy=Delta_x*Delta_x/4*(matrix.s[1,1] + matrix.s[Nx,1] + matrix.s[1,Nx] + matrix.s[Nx,Nx] + sum(2*matrix.s[1,(2:(Nx-1))]) + sum(2*matrix.s[(2:(Nx-1)),1]) + sum(2*matrix.s[Nx,(2:(Nx-1))]) + sum(2*matrix.s[(2:(Nx-1)),Nx]) + sum(4*matrix.s[(2:(Nx-1)),(2:(Nx-1))]) )
	
	sim_holo_grat %>% pull(rhom) %>% unique() -> rhom
	sim_holo_grat %>% pull(rhop) %>% unique() -> rhop
	sim_holo_grat %>% pull(rhob) %>% unique() -> rhob
	sim_holo_grat %>% pull(rhoz) %>% unique() -> rhoz
	
	volume1=b0/rhob + int_m_dx_dy/rhom + int_p_dx_dy/rhop + int_q_dx_dy/rhop + z0/rhoz + a0*int_a_dx_dy/rhoa + s0*int_s_dx_dy/rhos
	
	# 3.1 --- Diffusion of target a
	alpha_s_x=Ds*t0/x_hat/x_hat
	alpha_a_x=Da*t0/x_hat/x_hat
	alpha_s_y=Ds*t0/T_1/T_1
	alpha_a_y=Da*t0/T_1/T_1
	gamma_ss=gamma_s*z0*t0
	omega_ss=ifelse(s0==0, 0, omega_s*z0*t0/s0)
	gamma_sz=gamma_s*s0*t0
	omega_sz=omega_s*t0
	gamma_aa=gamma_a*z0*t0
	omega_aa=ifelse(a0==0, 0, omega_a*z0*t0/a0)
	gamma_az=gamma_a*a0*t0
	omega_az=omega_a*t0
	
	n_iterations = exposure_time/Delta_t + 1# Total number of iterations

	seq(0, exposure_time, by=output_time_step) -> time_vals
	
	if(exposure_time < output_time_step){
		bind_iter="0"
	} else {
		time_vals/Delta_t + 1 -> bind_iter
		bind_iter[-1] -> bind_iter
	}
	
	# Simulation of holographic grating exposed to a loaded solvent
	for(j in 1:n_iterations){
		
		AA2 = diag(2 + 2*r*alpha_a_x + 2*r*alpha_a_y + gamma_aa*Delta_t*as.numeric(ze1), Nx*Nx)
	  AA1 = diag(2 - 2*r*alpha_a_x - 2*r*alpha_a_y - gamma_aa*Delta_t*as.numeric(ze1), Nx*Nx)
	  
	  SS2 = diag(2 + 2*r*alpha_s_x + 2*r*alpha_s_y + gamma_ss*Delta_t*as.numeric(ze1), Nx*Nx)
	  SS1 = diag(2 - 2*r*alpha_s_x - 2*r*alpha_s_y - gamma_ss*Delta_t*as.numeric(ze1), Nx*Nx)
	  
	  ZEZE2 = diag(2 + gamma_az*Delta_t*as.numeric(a1) + gamma_sz*Delta_t*as.numeric(s1), Nx*Nx)
	  ZEZE1 = diag(2 - gamma_az*Delta_t*as.numeric(a1) - gamma_sz*Delta_t*as.numeric(s1), Nx*Nx)
	  
	  ZAZA2 = diag(2 + omega_az*Delta_t, Nx*Nx)
	  ZAZA1 = diag(2 - omega_az*Delta_t, Nx*Nx)
	  
	  ZSZS2 = diag(2 + omega_sz*Delta_t, Nx*Nx)
	  ZSZS1 = diag(2 - omega_sz*Delta_t, Nx*Nx)
	  
	  BB2 = diag(2, Nx*Nx)
		BB1 = diag(2, Nx*Nx)
		
		MM2 = diag(2, Nx*Nx)
		MM1 = diag(2, Nx*Nx)
		
		PP2 = diag(2, Nx*Nx)
		PP1 = diag(2, Nx*Nx)
		
		QQ2 = diag(2, Nx*Nx)
		QQ1 = diag(2, Nx*Nx)
	  
	  for(i in 1:(Nx*Nx)){
	    
	    i_minus_1=ifelse((i + Nx)%%Nx == 1, i + Nx - 2, i - 1)
	    i_plus_1=ifelse((i + Nx)%%Nx == 0, i - Nx + 2, i + 1)
	    j_minus_1=ifelse(i <= Nx, i + Nx, i - Nx)
	    j_plus_1=ifelse(i > (Nx-1)*Nx, i - Nx, i + Nx)
	    # i_minus_2=ifelse((i + Nx)%%Nx == 1, i + Nx - 3, ifelse((i + Nx)%%Nx == 2, i + Nx - 2, i - 2))
	    # i_plus_2=ifelse((i + Nx)%%Nx == 0, i - Nx + 3, ifelse((i + Nx)%%Nx == Nx-1, i - Nx + 2, i + 2))
	    # j_minus_2=ifelse(i <= Nx, i, ifelse(i <= 2*Nx, i-Nx, i - 2*Nx))
	    # j_plus_2=ifelse(i > (Nx-1)*Nx, i - 2*Nx, ifelse(i > (Nx-2)*Nx, i, i + 2*Nx))
	    
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
				    	
	  }
		
		a2 = solve(AA2) %*% ((AA1 %*% matrix(a1, ncol=1)) + matrix(2*Delta_t*omega_aa*za1, ncol=1))
		s2 = solve(SS2) %*% ((SS1 %*% matrix(s1, ncol=1)) + matrix(2*Delta_t*omega_ss*zs1, ncol=1))
	  ze2 = solve(ZEZE2) %*% ((ZEZE1 %*% matrix(ze1, ncol=1)) + matrix(2*omega_az*Delta_t*za1, ncol=1) + matrix(2*omega_sz*Delta_t*zs1, ncol=1))
	  za2 = solve(ZAZA2) %*% ((ZAZA1 %*% matrix(za1, ncol=1)) + matrix(2*Delta_t*gamma_az*a1*ze1, ncol=1))
	  zs2 = solve(ZSZS2) %*% ((ZSZS1 %*% matrix(zs1, ncol=1)) + matrix(2*Delta_t*gamma_sz*s1*ze1, ncol=1))
	  b2 = solve(BB2) %*% BB1 %*% matrix(b1, ncol=1)
	  m2 = solve(MM2) %*% MM1 %*% matrix(m1, ncol=1)
	  p2 = solve(PP2) %*% PP1 %*% matrix(p1, ncol=1)
	  q2 = solve(QQ2) %*% QQ1 %*% matrix(q1, ncol=1)
	  
	  a2[free_surface]=1
	  s2[free_surface]=1
	  
	  a2 -> a1
	  s2 -> s1
	  ze2 -> ze1
	  za2 -> za1
	  zs2 -> zs1
	  b2 -> b1
	  m2 -> m1
	  p2 -> p1
	  q2 -> q1
	  
	  if(as.character(j) %in% as.character(bind_iter)){
	  	
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
	df_a %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, a=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""), time=as.numeric(time), a=a0*a) %>% arrange(time,Y,x) -> df_a
	
	# Solvent mass data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),s) -> df_s
	names(df_s)=c("x", "Y", paste0("t",time_vals[1:ncol(s)]))# Rename columns
	df_s %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, s=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=as.numeric(time), s=s0*s) %>% arrange(time,Y,x) -> df_s
	
	# Empty nanozeolite mass data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),ze) -> df_ze
	names(df_ze)=c("x", "Y", paste0("t",time_vals[1:ncol(ze)]))# Rename columns
	df_ze %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, ze=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=as.numeric(time), ze=z0*ze) %>% arrange(time,Y,x) -> df_ze
	
	# Analyte-filled nanozeolite mass data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),za) -> df_za
	names(df_za)=c("x", "Y", paste0("t",time_vals[1:ncol(za)]))# Rename columns
	df_za %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, za=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=as.numeric(time), za=z0*za) %>% arrange(time,Y,x) -> df_za
	
	# Solvent-filled nanozeolite mass data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),zs) -> df_zs
	names(df_zs)=c("x", "Y", paste0("t",time_vals[1:ncol(zs)]))# Rename columns
	df_zs %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, zs=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=as.numeric(time), zs=z0*zs) %>% arrange(time,Y,x) -> df_zs
	
	# Binder data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),b) -> df_b
	names(df_b)=c("x", "Y", paste0("t",time_vals[1:ncol(b)]))# Rename columns
	df_b %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, b=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=as.numeric(time),b=b0*b) %>% arrange(time,Y,x) -> df_b
	
	# Monomer data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),m) -> df_m
	names(df_m)=c("x", "Y", paste0("t",time_vals[1:ncol(m)]))# Rename columns
	df_m %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, m=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=as.numeric(time),m=m) %>% arrange(time,Y,x) -> df_m
	
	# Short polymer data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),p) -> df_p
	names(df_p)=c("x", "Y", paste0("t",time_vals[1:ncol(p)]))# Rename columns
	df_p %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, p=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=as.numeric(time),p=p) %>% arrange(time,Y,x) -> df_p
	
	# Cross-linked polymer data frame
	data.frame(x=rep(seq(0,1,Delta_x), Nx),Y=sort(rep(seq(0,1,Delta_x), Nx)),q) -> df_q
	names(df_q)=c("x", "Y", paste0("t",time_vals[1:ncol(q)]))# Rename columns
	df_q %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, q=value) %>% mutate(time=gsub(x=time,pattern="t", replacement=""),time=as.numeric(time),q=q) %>% arrange(time,Y,x) -> df_q
	
	# 3.1 --- Refractive index matrix
  Vb = b*b0/rhob # cm**3
  Vm = m*m0/rhom # cm**3
  Vp = p*m0/rhop # cm**3
  Vq = q*m0/rhop # cm**3
  Vze = ze*z0/rhoz # cm**3
  Vza = za*z0/rhoz # cm**3
  Vzs = zs*z0/rhoz # cm**3
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
  df_n %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, n_after_exposure=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) %>% arrange(time,Y,x) -> df_n# Additional columns containing info on input
  
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
  	
  	Mean_RI[k]=Delta_x*Delta_x/4*(n1[1,1] + n1[Nx,1] + n1[1,Nx] + n1[Nx,Nx] + sum(2*n1[1,(2:(Nx-1))]) + sum(2*n1[(2:(Nx-1)),1]) + sum(2*n1[Nx,(2:(Nx-1))]) + sum(2*n1[(2:(Nx-1)),Nx]) + sum(4*n1[(2:(Nx-1)),(2:(Nx-1))]) )
    
    for(i in 1:Nx){
      
      n1[, i] -> n1i
      N0[k, i]=Delta_x/3*(n1i[1] + sum(2*n1i[times_2]) + sum(4*n1i[times_4]) + n1i[Nx])
      
      n1[, i]*cos(2*pi*x) -> n1cos
      N1a[k, i]=2*Delta_x/3*(n1cos[1] + sum(2*n1cos[times_2]) + sum(4*n1cos[times_4]) + n1cos[Nx])
      
      n1[, i]*sin(2*pi*x) -> n1sin
      N1b[k, i]=2*Delta_x/3*(n1sin[1] + sum(2*n1sin[times_2]) + sum(4*n1sin[times_4]) + n1sin[Nx])
      
    }
    
    
    n.Fourier=matrix(0, nrow=Nx, ncol=Nx)
    
    for(i in 1:Nx){
      
      n.Fourier[,i] = N0[k, i] + N1a[k, i]*cos(2*pi*x) + N1b[k, i]*sin(2*pi*x)
      
    }
    
    diff = abs(n1 - n.Fourier)**2
    
    for(i in 1:Nx){
    	d2[k,i] = Delta_x/3*(diff[1, i] + sum(2*diff[times_2, i]) + sum(4*diff[times_4, i]) + diff[Nx, i])
    	int_n2[k,i] = Delta_x/3*(n2[1, i] + sum(2*n2[times_2, i]) + sum(4*n2[times_4, i]) + n2[Nx, i])
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
  	theta_B = asin(lambda_probe/2/Mean_RI/Lambda_1) - phi_r_1,
  	Delta_theta_B = theta_B1 - theta_B,
  	Delta_phi_r = 0
  ) -> df1
  
  if(df0 %>% subset(Y==0) %>% distinct(Klein.Cook) %>% pull(Klein.Cook) < 10){
	
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