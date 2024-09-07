# Slanted_Grating_Model_v21.R
# Author: Jack Lyons

# A function to run numerical simulation of mathematical model describing the formation of slanted holographic gratings.

# 0.1 --- Import libraries
library(dplyr)
library(reshape2)
library(ggplot2)

# 0.2 --- Increase memory limit
memory.limit(size=1e6)

# 0.3 --- Set working directory
setwd("~/PhD Documents/Inert Zeolite Nanoparticles Model/Slanted Gratings")

# 0.4 --- Plot aesthetics
load("~/PhD Documents/plt.theme1.RData")
load("~/PhD Documents/plt.theme2.RData")

# 1.1 --- Define function inputs and default values.
slanted_grating_simulation_v21 = function(
    
  start_exp=0,# Start of exposure
  end_exp=1e2, # End of exposure
  total.time=100,# Total simulation time
  lpmm=1e3,# Spatial frequency
  I0=5,# Intensity of recording beam
  slant.angle=1e-2,# Grating slant angle
  xi=0.3,# Scattering coefficient
  nm=1.55,# Monomer refractive index
  rhom=1.15,# Monomer density
  Dm=1.6e-7,# Monomer diffusion coefficient
	Dp=4e-3*1.6e-7,# Polymer diffusion coefficient
	rhop=1.3,# Polymer density
	rhoq=1.35,# Polymer density
	np=1.56,# Oligomer refractive index
  nq=1.64,# Polymer refractive index
  Gamma=1,# Rate of immobilization
  wt.pc=5e-2,# Doping %
  Dz=1e-8,# Nanodopant self-diffusion coefficient
  epsilon.mz=0,# Cross-diffusion ratio
  epsilon.pz=11,# Cross-diffusion ratio
	epsilon.qz=11,# Cross-diffusion ratio
	rhoz=1.74,# Nanodopant mass density
  nz=1.366,# Nanodopant refractive index
  b0=5.05,# Ratio of binder to monomer mass
  nb=1.5,# Binder refractive index
  rhob=1.19, # Binder mass density
  T0=50e-4,# Depth of photosensitive layer [cm]
	zeta=139,# absorption coefficient [cm**-1]
	lambda.probe=633e-7,# Wavelength of reconstruction beam
  Delta.t=1/100,# Numerical scheme time step
  Delta.x=1/20,# Numerical scheme spatial step
  output.time.step=1# Output time
  
){
  
  # 1.2 --- Define parameters
	Delta.Y=Delta.x
  iterations.per.second=1/Delta.t# Number of iterations each second
  Nx=1/Delta.x + 1# Number of spatial points
  Ny=1/Delta.Y + 1# Number of spatial points
  if(Nx%%2==0){return(warning("Number of x mesh points must be an odd number."))}
  if(Ny%%2==0){return(warning("Number of y mesh points must be an odd number."))}
  x=seq(0, 1, length.out=Nx)# Non-dimensional grating distance
  y=seq(0,1,length.out=Ny)# Non-dimensional depth
  Nt = total.time*iterations.per.second+1# Total number of iterations
  r=Delta.t/Delta.x/Delta.x# Ratio of finite time step to squared finite spatial step
  m0=1# Initial mass of monomer
  t0=1 #  Reference time [s]
  Lambda0=1/10/lpmm # Grating period [cm]
  Lambda1=Lambda0
  Lambda.t=c(Lambda0)
  j_start_exp=start_exp/Delta.t # Iteration of exposure start
  j_end_exp=end_exp/Delta.t # Iteration of exposure end
  z0 = wt.pc/(1 - wt.pc)*(m0 + b0)# Initial nanodopant to monomer
  
  # 1.3 --- Matrix initial conditions
  u1=1
  u.t=c(u1)
  du.dt=0
  rep(1, Nx*Nx) -> m1# m at j=0
  rep(0, Nx*Nx) -> p1# p at j=0
  rep(0, Nx*Nx) -> q1# q at j=0
  rep(1, Nx*Nx) -> z1# z at j=0
  rep(1, Nx*Nx) -> b1# b at j=0
  
  matrix(m1, ncol=1) -> m# monomer density matrix
  matrix(p1, ncol=1) -> p# oligomer density matrix
  matrix(q1, ncol=1) -> q# Immobile polymer matrix
  matrix(z1, ncol=1) -> z# Zeolite matrix
  matrix(b1, ncol=1) -> b# Binder density matrix
  
  Volume0=m0/rhom + b0/rhob + z0/rhoz
  Volume.t=c(Volume0)
  
  phi.m0=m0/rhom/Volume0
  phi.z0=z0/rhoz/Volume0
  phi.b0=b0/rhob/Volume0
  
  Lorentz.Lorenz.RHS = phi.m0*(nm*nm - 1)/(nm*nm + 2) + phi.b0*(nb*nb - 1)/(nb*nb + 2) + phi.z0*(nz*nz - 1)/(nz*nz + 2)
  Initial.RI=sqrt((2*Lorentz.Lorenz.RHS + 1)/(1 - Lorentz.Lorenz.RHS))
  slant.angle=ifelse(slant.angle==0,1e-5,slant.angle)
  phi.0=slant.angle/180*pi
  phi.r0=asin(sin(phi.0)/Initial.RI)
  phi.r1=phi.r0
  phi.r.t=c(phi.r0)
  theta_B0=asin(lambda.probe/2/Initial.RI/Lambda0) - phi.r0
  y.hat0=Lambda0/sin(phi.r0)
  y.hat1=y.hat0
  y.hat.t=c(y.hat0)
  x.hat=Lambda0/cos(phi.r0)
  
  # 1.4 --- time step iterations
  seq(0, total.time, by=output.time.step) -> time.vals
  time.vals/Delta.t + 1 -> bind.iter
  bind.iter = as.character(bind.iter[-1])
  
  # 1.5 --- Nondimensionalized parameters
  alpha.m.x=Dm*t0/x.hat/x.hat
  alpha.m.y=Dm*t0/T0/T0
  alpha.p.x=Dp*t0/x.hat/x.hat
  alpha.p.y=Dp*t0/T0/T0
  alpha.z.x=Dz*t0/x.hat/x.hat
  alpha.z.y=Dz*t0/T0/T0
  F0=0.1*I0**0.3
  beta=F0*t0
  zeta_star=zeta*T0
  T_star=T0/x.hat
  xi_star=xi*z0
  #beta=matrix(0, nrow=Nx, ncol=Nx)
  #for(k in 1:Nx){beta[,k]=F0*t0*exp(-zeta_star*Y[k])}
  #beta=melt(beta)$value
  #gamma=Gamma*t0
  
  interior.points = 2:(Nx-1)
  times.4 = interior.points[interior.points %% 2 == 0]
  times.2 = interior.points[interior.points %% 2 != 0]
  Y=seq(0,1,length.out=Ny)
  Y1=sort(rep(Y, Nx))
  int_M=c(1)
  int_P=c(0)
  int_Q=c(0)
  int_Z=c(z0)
  int_B=c(b0)
  Actual.Shrinkage=c(0)
  
  # 1.6 --- Calculate each time step via implicit finite difference method
  for(j in 1:Nt){
    
    if(j >= j_start_exp & j <= j_end_exp){Phi=1} else {Phi = 0}# Phi=1 if illumination is on, 0 otherwise
    
  	f = matrix(0, nrow=Nx, ncol=Nx)
  	#gamma=matrix(0, nrow=Nx, ncol=Nx)
    for(i in 1:Nx){
      matrix.z1=matrix(z1, nrow=Nx, ncol=Nx)
      z1.i=as.numeric(matrix.z1[,i])
      f[,i] = exp(-0.3*zeta_star*u1*(1-Y[i]))*(1 + exp(-xi_star*z1.i)*cos(2*pi*x - 2*pi*T_star*tan(phi.r1)*u1*Y[i]))
      #gamma = Gamma*t0*exp(-zeta_star*u1*Y[i])
    }
    f=melt(f)$value
    #gamma=melt(gamma)$value
    gamma=Gamma*m0*t0
    
    MM2 = diag(2 + Phi*Delta.t*beta*f, Nx*Nx)
    MM1 = diag(2 - Phi*Delta.t*beta*f, Nx*Nx)
    PP2 = diag(2 + Delta.t*gamma*as.numeric(p1), Nx*Nx)
    PP1 = diag(2 - Delta.t*gamma*as.numeric(p1), Nx*Nx)
    PM2 = diag(+Phi*Delta.t*beta*f, Nx*Nx)
    PM1 = diag(+Phi*Delta.t*beta*f, Nx*Nx)
    QQ2=diag(2, Nx*Nx)
    QQ1=diag(2, Nx*Nx)
    QP2=diag(+gamma*Delta.t*as.numeric(p1),Nx*Nx)
    QP1=diag(+gamma*Delta.t*as.numeric(p1),Nx*Nx)
    ZZ2 = diag(2, Nx*Nx)
    ZZ1 = diag(2, Nx*Nx)
    BB2 = diag(2, Nx*Nx)
    BB1 = diag(2, Nx*Nx)
    
    for(i in 1:(Nx*Nx)){
      
      i.minus.1=ifelse((i + Nx)%%Nx == 1, i + Nx - 2, i - 1)
      
      i.plus.1=ifelse((i + Nx)%%Nx == 0, i - Nx + 2, i + 1)
      
      j.minus.1=ifelse(i < Nx+1, i + Nx, i - Nx)
      
      j.plus.1=ifelse(i > (Nx-1)*Nx, i - Nx, i + Nx)
      
      # i.minus.2=ifelse((i + Nx)%%Nx == 1, i + Nx - 3, ifelse((i + Nx)%%Nx == 2, i + Nx - 2, i - 2))
      # 
      # i.plus.2=ifelse((i + Nx)%%Nx == 0, i - Nx + 3, ifelse((i + Nx)%%Nx == Nx-1, i - Nx + 2, i + 2))
      # 
      # j.minus.2=ifelse(i <= Nx, i + 2*Nx, ifelse(i <= 2*Nx, i, i - 2*Nx))
      # 
      # j.plus.2=ifelse(i > (Nx-1)*Nx, i - 2*Nx, ifelse(i > (Nx-2)*Nx, i, i + 2*Nx))
      
      MM2[i, i.minus.1] = MM2[i, i.minus.1] - r*alpha.m.x*(1 + epsilon.mz*z0*z1[i.minus.1])
	    MM2[i, j.minus.1] = MM2[i, j.minus.1] - r*alpha.m.y/u1/u1*(1 + epsilon.mz*z0*z1[j.minus.1])
	    MM2[i, i] = MM2[i, i] + 2*r*alpha.m.x*(1 + epsilon.mz*z0*z1[i])
	    MM2[i, i] = MM2[i, i] + 2*r*alpha.m.y/u1/u1*(1 + epsilon.mz*z0*z1[i])
	    MM2[i, i.plus.1] = MM2[i, i.plus.1] - r*alpha.m.x*(1 + epsilon.mz*z0*z1[i.plus.1])
	    MM2[i, j.plus.1] = MM2[i, j.plus.1] - r*alpha.m.y/u1/u1*(1 + epsilon.mz*z0*z1[j.plus.1])
	    
	    MM1[i, i.minus.1] = MM1[i, i.minus.1] + r*alpha.m.x*(1 + epsilon.mz*z0*z1[i.minus.1])
	    MM1[i, j.minus.1] = MM1[i, j.minus.1] + r*alpha.m.y/u1/u1*(1 + epsilon.mz*z0*z1[j.minus.1])
	    MM1[i, i] = MM1[i, i] - 2*r*alpha.m.x*(1 + epsilon.mz*z0*z1[i])
	    MM1[i, i] = MM1[i, i] - 2*r*alpha.m.y/u1/u1*(1 + epsilon.mz*z0*z1[i])
	    MM1[i, i.plus.1] = MM1[i, i.plus.1] + r*alpha.m.x*(1 + epsilon.mz*z0*z1[i.plus.1])
	    MM1[i, j.plus.1] = MM1[i, j.plus.1] + r*alpha.m.y/u1/u1*(1 + epsilon.mz*z0*z1[j.plus.1])
      
      PP2[i, i.minus.1] = PP2[i, i.minus.1] - r*alpha.p.x*(1 + epsilon.pz*z0*z1[i.minus.1])
	    PP2[i, j.minus.1] = PP2[i, j.minus.1] - r*alpha.p.y/u1/u1*(1 + epsilon.pz*z0*z1[j.minus.1])
	    PP2[i, i] = PP2[i, i] + 2*r*alpha.p.x*(1 + epsilon.pz*z0*z1[i])
	    PP2[i, i] = PP2[i, i] + 2*r*alpha.p.y/u1/u1*(1 + epsilon.pz*z0*z1[i])
	    PP2[i, i.plus.1] = PP2[i, i.plus.1] - r*alpha.p.x*(1 + epsilon.pz*z0*z1[i.plus.1])
	    PP2[i, j.plus.1] = PP2[i, j.plus.1] - r*alpha.p.y/u1/u1*(1 + epsilon.pz*z0*z1[j.plus.1])
	    
	    PP1[i, i.minus.1] = PP1[i, i.minus.1] + r*alpha.p.x*(1 + epsilon.pz*z0*z1[i.minus.1])
	    PP1[i, j.minus.1] = PP1[i, j.minus.1] + r*alpha.p.y/u1/u1*(1 + epsilon.pz*z0*z1[j.minus.1])
	    PP1[i, i] = PP1[i, i] - 2*r*alpha.p.x*(1 + epsilon.pz*z0*z1[i])
	    PP1[i, i] = PP1[i, i] - 2*r*alpha.p.y/u1/u1*(1 + epsilon.pz*z0*z1[i])
	    PP1[i, i.plus.1] = PP1[i, i.plus.1] + r*alpha.p.x*(1 + epsilon.pz*z0*z1[i.plus.1])
	    PP1[i, j.plus.1] = PP1[i, j.plus.1] + r*alpha.p.y/u1/u1*(1 + epsilon.pz*z0*z1[j.plus.1])
      
      ZZ2[i, i.minus.1] = ZZ2[i, i.minus.1] - r*alpha.z.x*(1 + epsilon.qz*q1[i.minus.1])
	    ZZ2[i, j.minus.1] = ZZ2[i, j.minus.1] - r*alpha.z.y/u1/u1*(1 + epsilon.qz*q1[j.minus.1])
	    ZZ2[i, i] = ZZ2[i, i] + 2*r*alpha.z.x*(1 + epsilon.qz*q1[i])
	    ZZ2[i, i] = ZZ2[i, i] + 2*r*alpha.z.y/u1/u1*(1 + epsilon.qz*q1[i])
	    ZZ2[i, i.plus.1] = ZZ2[i, i.plus.1] - r*alpha.z.x*(1 + epsilon.qz*q1[i.plus.1])
	    ZZ2[i, j.plus.1] = ZZ2[i, j.plus.1] - r*alpha.z.y/u1/u1*(1 + epsilon.qz*q1[j.plus.1])
	    
	    ZZ1[i, i.minus.1] = ZZ1[i, i.minus.1] + r*alpha.z.x*(1 + epsilon.qz*q1[i.minus.1])
	    ZZ1[i, j.minus.1] = ZZ1[i, j.minus.1] + r*alpha.z.y/u1/u1*(1 + epsilon.qz*q1[j.minus.1])
	    ZZ1[i, i] = ZZ1[i, i] - 2*r*alpha.z.x*(1 + epsilon.qz*q1[i])
	    ZZ1[i, i] = ZZ1[i, i] - 2*r*alpha.z.y/u1/u1*(1 + epsilon.qz*q1[i])
	    ZZ1[i, i.plus.1] = ZZ1[i, i.plus.1] + r*alpha.z.x*(1 + epsilon.qz*q1[i.plus.1])
	    ZZ1[i, j.plus.1] = ZZ1[i, j.plus.1] + r*alpha.z.y/u1/u1*(1 + epsilon.qz*q1[j.plus.1])

    	MM2[i, j.minus.1] = MM2[i, j.minus.1] + Y1[i]/u1/2/Delta.Y*du.dt
    	MM1[i, j.minus.1] = MM1[i, j.minus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	MM2[i, j.plus.1] = MM2[i, j.plus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	MM1[i, j.plus.1] = MM1[i, j.plus.1] + Y1[i]/u1/2/Delta.Y*du.dt

    	PP2[i, j.minus.1] = PP2[i, j.minus.1] + Y1[i]/u1/2/Delta.Y*du.dt
    	PP1[i, j.minus.1] = PP1[i, j.minus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	PP2[i, j.plus.1] = PP2[i, j.plus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	PP1[i, j.plus.1] = PP1[i, j.plus.1] + Y1[i]/u1/2/Delta.Y*du.dt
    	
    	QQ2[i, j.minus.1] = QQ2[i, j.minus.1] + Y1[i]/u1/2/Delta.Y*du.dt
    	QQ1[i, j.minus.1] = QQ1[i, j.minus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	QQ2[i, j.plus.1] = QQ2[i, j.plus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	QQ1[i, j.plus.1] = QQ1[i, j.plus.1] + Y1[i]/u1/2/Delta.Y*du.dt

    	ZZ2[i, j.minus.1] = ZZ2[i, j.minus.1] + Y1[i]/u1/2/Delta.Y*du.dt
    	ZZ1[i, j.minus.1] = ZZ1[i, j.minus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	ZZ2[i, j.plus.1] = ZZ2[i, j.plus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	ZZ1[i, j.plus.1] = ZZ1[i, j.plus.1] + Y1[i]/u1/2/Delta.Y*du.dt
    	
    	BB2[i, j.minus.1] = BB2[i, j.minus.1] + Y1[i]/u1/2/Delta.Y*du.dt
    	BB1[i, j.minus.1] = BB1[i, j.minus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	BB2[i, j.plus.1] = BB2[i, j.plus.1] - Y1[i]/u1/2/Delta.Y*du.dt
    	BB1[i, j.plus.1] = BB1[i, j.plus.1] + Y1[i]/u1/2/Delta.Y*du.dt

	    
    }
    
    m2 = solve(MM2) %*% MM1 %*% matrix(m1, ncol=1)
    
    p2 = solve(PP2) %*% ((PP1 %*% matrix(p1, ncol=1)) + (PM2 %*% matrix(m2, ncol=1)) + (PM1 %*% matrix(m1, ncol=1)))
    
    q2 = solve(QQ2) %*% ((QQ1 %*% matrix(q1, ncol=1)) + (QP2 %*% matrix(p2, ncol=1)) + (QP1 %*% matrix(p1, ncol=1)))
    
    if(z0==0){
    	z2 = matrix(z1, ncol=1)
    } else {
    	z2 = solve(ZZ2) %*% ZZ1 %*% matrix(z1, ncol=1)
    }
    
    b2 = solve(BB2) %*% BB1 %*% matrix(b1, ncol=1)
    
    
    dm.dt = matrix((m2-m1)/Delta.t, ncol=Nx)
    dp.dt = matrix((p2-p1)/Delta.t, ncol=Nx)
    dq.dt = matrix((q2-q1)/Delta.t, ncol=Nx)
    dz.dt = matrix(z0*(z2-z1)/Delta.t, ncol=Nx)
    
    int_dmdt=Delta.x*Delta.x/4*(dm.dt[1,1] + dm.dt[Nx,1] + dm.dt[1,Nx] + dm.dt[Nx,Nx] + sum(2*dm.dt[1,(2:(Nx-1))]) + sum(2*dm.dt[(2:(Nx-1)),1]) + sum(2*dm.dt[Nx,(2:(Nx-1))]) + sum(2*dm.dt[(2:(Nx-1)),Nx]) + sum(4*dm.dt[(2:(Nx-1)),(2:(Nx-1))]) )
    
    int_dpdt=Delta.x*Delta.x/4*(dp.dt[1,1] + dp.dt[Nx,1] + dp.dt[1,Nx] + dp.dt[Nx,Nx] + sum(2*dp.dt[1,(2:(Nx-1))]) + sum(2*dp.dt[(2:(Nx-1)),1]) + sum(2*dp.dt[Nx,(2:(Nx-1))]) + sum(2*dp.dt[(2:(Nx-1)),Nx]) + sum(4*dp.dt[(2:(Nx-1)),(2:(Nx-1))]) )
    
    int_dqdt=Delta.x*Delta.x/4*(dq.dt[1,1] + dq.dt[Nx,1] + dq.dt[1,Nx] + dq.dt[Nx,Nx] + sum(2*dq.dt[1,(2:(Nx-1))]) + sum(2*dq.dt[(2:(Nx-1)),1]) + sum(2*dq.dt[Nx,(2:(Nx-1))]) + sum(2*dq.dt[(2:(Nx-1)),Nx]) + sum(4*dq.dt[(2:(Nx-1)),(2:(Nx-1))]) )
    
    int_dzdt=Delta.x*Delta.x/4*(dz.dt[1,1] + dz.dt[Nx,1] + dz.dt[1,Nx] + dz.dt[Nx,Nx] + sum(2*dz.dt[1,(2:(Nx-1))]) + sum(2*dz.dt[(2:(Nx-1)),1]) + sum(2*dz.dt[Nx,(2:(Nx-1))]) + sum(2*dz.dt[(2:(Nx-1)),Nx]) + sum(4*dz.dt[(2:(Nx-1)),(2:(Nx-1))]) )
    
    m2 -> m1
    p2 -> p1
    q2 -> q1
    z2 -> z1
    b2 -> b1
    
    m.t=matrix(m0*m1, ncol=Nx)
    p.t=matrix(m0*p1, ncol=Nx)
    q.t=matrix(m0*q1, ncol=Nx)
    z.t=matrix(z0*z1, ncol=Nx)
    b.t=matrix(b0*b1, ncol=Nx)
    
    int_M1=Delta.x*Delta.x/4*(m.t[1,1] + m.t[Nx,1] + m.t[1,Nx] + m.t[Nx,Nx] + sum(2*m.t[1,(2:(Nx-1))]) + sum(2*m.t[(2:(Nx-1)),1]) + sum(2*m.t[Nx,(2:(Nx-1))]) + sum(2*m.t[(2:(Nx-1)),Nx]) + sum(4*m.t[(2:(Nx-1)),(2:(Nx-1))]) )
    
    int_P1=Delta.x*Delta.x/4*(p.t[1,1] + p.t[Nx,1] + p.t[1,Nx] + p.t[Nx,Nx] + sum(2*p.t[1,(2:(Nx-1))]) + sum(2*p.t[(2:(Nx-1)),1]) + sum(2*p.t[Nx,(2:(Nx-1))]) + sum(2*p.t[(2:(Nx-1)),Nx]) + sum(4*p.t[(2:(Nx-1)),(2:(Nx-1))]) )

    int_Q1=Delta.x*Delta.x/4*(q.t[1,1] + q.t[Nx,1] + q.t[1,Nx] + q.t[Nx,Nx] + sum(2*q.t[1,(2:(Nx-1))]) + sum(2*q.t[(2:(Nx-1)),1]) + sum(2*q.t[Nx,(2:(Nx-1))]) + sum(2*q.t[(2:(Nx-1)),Nx]) + sum(4*q.t[(2:(Nx-1)),(2:(Nx-1))]) )

    int_Z1=Delta.x*Delta.x/4*(z.t[1,1] + z.t[Nx,1] + z.t[1,Nx] + z.t[Nx,Nx] + sum(2*z.t[1,(2:(Nx-1))]) + sum(2*z.t[(2:(Nx-1)),1]) + sum(2*z.t[Nx,(2:(Nx-1))]) + sum(2*z.t[(2:(Nx-1)),Nx]) + sum(4*z.t[(2:(Nx-1)),(2:(Nx-1))]) )
    
    int_B1=Delta.x*Delta.x/4*(b.t[1,1] + b.t[Nx,1] + b.t[1,Nx] + b.t[Nx,Nx] + sum(2*b.t[1,(2:(Nx-1))]) + sum(2*b.t[(2:(Nx-1)),1]) + sum(2*b.t[Nx,(2:(Nx-1))]) + sum(2*b.t[(2:(Nx-1)),Nx]) + sum(4*b.t[(2:(Nx-1)),(2:(Nx-1))]) )
    
    du.dt = (int_dmdt/rhom + int_dpdt/rhop + int_dqdt/rhoq + int_dzdt/rhoz)/Volume0
    
    u2 = u1 + Delta.t*du.dt; u2 -> u1
    
    phi.r1=atan(tan(phi.r0)/u1)
    
    Lambda1=cos(phi.r1)/cos(phi.r0)*Lambda0
    
    y.hat1=Lambda1/sin(phi.r1)
    
    if(as.character(j) %in% bind.iter){
      
      matrix(cbind(m, m1), nrow=Nx*Nx) -> m
      matrix(cbind(p, p1), nrow=Nx*Nx) -> p
      matrix(cbind(q, q1), nrow=Nx*Nx) -> q
      matrix(cbind(z, z1), nrow=Nx*Nx) -> z
      matrix(cbind(b, b1), nrow=Nx*Nx) -> b
      Actual.Shrinkage=c(Actual.Shrinkage,1-u1)
      u.t=c(u.t,u1)
      Lambda.t=c(Lambda.t,Lambda1)
      phi.r.t=c(phi.r.t,phi.r1)
      y.hat.t=c(y.hat.t,y.hat1)
      int_M=c(int_M,int_M1)
      int_P=c(int_P,int_P1)
      int_Q=c(int_Q,int_Q1)
      int_Z=c(int_Z,int_Z1)
      int_B=c(int_B,int_B1)
    }
    
    
  }
  
  # 2.1 --- Monomer density data frame
  data.frame(x=rep(seq(0,1,Delta.x), Nx), Y=sort(rep(seq(0,1,Delta.x), Nx)), m) -> df.m# Turn the matrix into a data frame.
  names(df.m)=c("x","Y", paste0("t",time.vals[1:ncol(m)]))# Rename columns
  df.m %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, mass.m=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) -> melt.df.m# Additional columns containing info on input
  
  # 2.2 --- Oligomer density data frame
  data.frame(x=rep(seq(0,1,Delta.x), Nx), Y=sort(rep(seq(0,1,Delta.x), Nx)),p) -> df.p# Turn the matrix into a data frame.
  names(df.p)=c("x","Y", paste0("t",time.vals[1:ncol(p)]))# Rename columns
  df.p %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, mass.p=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) -> melt.df.p# Additional columns containing info on input
  
  # 2.3 --- Polymer density data frame
  data.frame(x=rep(seq(0,1,Delta.x), Nx), Y=sort(rep(seq(0,1,Delta.x), Nx)), q) -> df.q# Turn the matrix into a data frame.
  names(df.q)=c("x","Y", paste0("t",time.vals[1:ncol(q)]))# Rename columns
  df.q %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, mass.q=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) -> melt.df.q# Additional columns containing info on input
  
  # 2.4 --- Zeolite density data frame
  data.frame(x=rep(seq(0,1,Delta.x), Nx), Y=sort(rep(seq(0,1,Delta.x), Nx)), z) -> df.z# Turn the matrix into a data frame.
  names(df.z)=c("x","Y", paste0("t",time.vals[1:ncol(z)]))# Rename columns
  df.z %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, mass.z=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time),mass.z=z0*mass.z) -> melt.df.z# Additional columns containing info on input
  
  data.frame(x=rep(seq(0,1,Delta.x), Nx), Y=sort(rep(seq(0,1,Delta.x), Nx)), b) -> df.b# Turn the matrix into a data frame.
  names(df.b)=c("x","Y", paste0("t",time.vals[1:ncol(b)]))# Rename columns
  df.b %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, mass.b=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time),mass.b=b0*mass.b) -> melt.df.b# Additional columns containing info on input
  
  # 2.5 --- Element density data frame
  merge(
    melt.df.m %>% mutate(
      x=as.character(x),
      Y=as.character(Y),
      time=as.character(time)
    ),
    merge(
      melt.df.p %>% mutate(
        x=as.character(x),
        Y=as.character(Y),
        time=as.character(time)
      ),
      merge(
        melt.df.q %>% mutate(
          x=as.character(x),
          Y=as.character(Y),
          time=as.character(time)
        ),
        merge(
      	melt.df.z %>% mutate(
          x=as.character(x),
          Y=as.character(Y),
          time=as.character(time)
        ),
        	melt.df.b %>% mutate(
	          x=as.character(x),
	          Y=as.character(Y),
	          time=as.character(time)
	      	)
        )
      )
    )
  ) %>% mutate(
    x=as.numeric(x),
    Y=as.numeric(Y),
    time=as.numeric(time)
  ) -> df
  
  # 3.1 --- Refractive index matrix
  Vb = b*b0/rhob # cm**3
  Vm = m*m0/rhom # cm**3
  Vp = p*m0/rhop # cm**3
  Vq = q*m0/rhoq # cm**3
  Vz = z*z0/rhoz # cm**3
  Vtotal=Vm+Vp+Vq+Vz+Vb# cm**3
  
  phi.m=Vm/Vtotal
  phi.p=Vp/Vtotal
  phi.q=Vq/Vtotal
  phi.z=Vz/Vtotal
  phi.b=Vb/Vtotal
  
  Lorentz.Lorenz.RHS = phi.m*(nm*nm - 1)/(nm*nm + 2) + phi.b*(nb*nb - 1)/(nb*nb + 2) + phi.p*(np*np - 1)/(np*np + 2) + phi.q*(nq*nq - 1)/(nq*nq + 2) + phi.z*(nz*nz - 1)/(nz*nz + 2)
  n=sqrt((2*Lorentz.Lorenz.RHS + 1)/(1 - Lorentz.Lorenz.RHS))
  
  df.n=data.frame(x=rep(seq(0,1,Delta.x), Nx), Y=sort(rep(seq(0,1,Delta.x), Nx)), n)# Turn the matrix into a data frame.
  names(df.n)=c("x","Y", paste0("t",time.vals[1:ncol(n)]))# Rename columns
  df.n %>% melt(id.vars=c("x","Y")) %>% rename(time=variable, n=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) -> melt.df.n# Additional columns containing info on input
  
  df.2=merge(
    df %>% mutate(
      x=as.character(x),
      Y=as.character(Y),
      time=as.character(time)
    ),
    melt.df.n %>% mutate(
      x=as.character(x),
      Y=as.character(Y),
      time=as.character(time)
    ),
    by=c("x","Y","time")
  ) %>% mutate(
    x=as.numeric(x),
    Y=as.numeric(Y),
    time=as.numeric(time)
  ) %>% arrange(time,Y,x)
  
  # 4.1 --- Refractive index modulation data frame
  N = ncol(n)
  N0=matrix(0, nrow=N, ncol=Nx)
  N1a=matrix(0, nrow=N, ncol=Nx)
  N1b=matrix(0, nrow=N, ncol=Nx)
  d2=matrix(0, nrow=N, ncol=Nx)
  int_n2=matrix(0, nrow=N, ncol=Nx)
  
  for(k in 1:N){
    
    matrix(n[, k], nrow=Nx) -> n1
  	n2=n1**2
    
    for(i in 1:Nx){
      
      n1[, i] -> n1i
      N0[k, i]=Delta.x/3*(n1i[1] + sum(2*n1i[times.2]) + sum(4*n1i[times.4]) + n1i[Nx])
      
      #n1[, i]*cos(2*pi*x - 2*pi*T_star*tan(phi.r.t[k])*u.t[k]*Y[i]) -> n1cos
      n1[, i]*cos(2*pi*x) -> n1cos
      N1a[k, i]=2*Delta.x/3*(n1cos[1] + sum(2*n1cos[times.2]) + sum(4*n1cos[times.4]) + n1cos[Nx])
      
      #n1[, i]*sin(2*pi*x - 2*pi*T_star*tan(phi.r.t[k])*u.t[k]*Y[i]) -> n1sin
      n1[, i]*sin(2*pi*x) -> n1sin
      N1b[k, i]=2*Delta.x/3*(n1sin[1] + sum(2*n1sin[times.2]) + sum(4*n1sin[times.4]) + n1sin[Nx])
      
    }
    
    
    n.Fourier=matrix(0, nrow=Nx, ncol=Nx)
    
    for(i in 1:Nx){
      
      #n.Fourier[,i] = N0[k, i] + N1a[k, i]*cos(2*pi*x - 2*pi*T_star*tan(phi.r.t[k])*u.t[k]*Y[i]) + N1b[k, i]*sin(2*pi*x - 2*pi*T_star*tan(phi.r.t[k])*u.t[k]*Y[i])
    	
    	n.Fourier[,i] = N0[k, i] + N1a[k, i]*cos(2*pi*x) + N1b[k, i]*sin(2*pi*x)
      
    }
    
    diff = abs(n1 - n.Fourier)**2
    
    for(i in 1:Nx){
    	d2[k,i] = Delta.x/3*(diff[1, i] + sum(2*diff[times.2, i]) + sum(4*diff[times.4, i]) + diff[Nx, i])
    	int_n2[k,i] = Delta.x/3*(n2[1, i] + sum(2*n2[times.2, i]) + sum(4*n2[times.4, i]) + n2[Nx, i])
    }
    
  }# End for loop
  
  delta=d2/int_n2
  
  data.frame(
    time=rep(seq(0,total.time,output.time.step), Nx),
    Y=sort(rep(seq(0,1,Delta.x), N)),
    N0=melt(N0)$value,
  	N1a=melt(N1a)$value,
  	N1b=melt(N1b)$value,
    delta=melt(delta)$value
  ) -> df.3
  
  merge(
    df.2 %>% mutate(time=as.character(time), Y=as.character(Y)),
    df.3 %>% mutate(time=as.character(time), Y=as.character(Y)), 
    by=c("time","Y")
  ) -> df.4
  
  data.frame(
    time=seq(0,total.time,output.time.step), 
    T.t=T0*u.t,
  	Actual.Shrinkage,
    Lambda.t=Lambda.t,
  	phi.r.t=phi.r.t,
  	int_M=int_M,
  	int_P=int_P,
  	int_Q=int_Q,
  	int_Z=int_Z,
  	int_B=int_B
  ) -> df.5
  
  merge(df.4 %>% mutate(time=as.character(time)),
  			df.5 %>% mutate(time=as.character(time)), 
    		by="time") -> df.6
  
  # Refractive index harmonics as function of y and t
  mean_RI=rep(0,N)
  
  for(k in 1:N){
  	
  	matrix(n[, k], nrow=Nx) -> n1
  	
  	mean_RI[k]=Delta.x*Delta.x/4*(n1[1,1] + n1[Nx,1] + n1[1,Nx] + n1[Nx,Nx] + sum(2*n1[1,(2:(Nx-1))]) + sum(2*n1[(2:(Nx-1)),1]) + sum(2*n1[Nx,(2:(Nx-1))]) + sum(2*n1[(2:(Nx-1)),Nx]) + sum(4*n1[(2:(Nx-1)),(2:(Nx-1))]) )
  	
  }
  
   merge(
    df.6 %>% mutate(time=as.character(time)),
    data.frame(time=as.character(seq(0,total.time,output.time.step)),mean_RI),
    by="time"
  ) %>% mutate(
  	theta_B=asin(lambda.probe/2/mean_RI/Lambda.t) - phi.r.t,
  	Delta_theta_B=theta_B0 - theta_B,
  	Delta_phi_r=phi.r0 - phi.r.t,
  	Apparent.Shrinkage=1 - tan(phi.r0)/tan(phi.r0 + Delta_theta_B)
  	) -> df.8
  
  df.8 %>%
    mutate(
      time=as.numeric(time),
      Y=as.numeric(Y),
    	y=T0*u.t*Y,
    	start_exp=start_exp,
      end_exp=end_exp,
      total.time=total.time,
      lpmm=lpmm,
      I0=I0,
      slant.angle=ifelse(abs(slant.angle) < 1e-3, 0, slant.angle),
      xi=xi,
      Dm=Dm,
    	Dp=Dp,
      epsilon.mz=epsilon.mz,
    	epsilon.pz=epsilon.pz,
    	epsilon.qz=epsilon.qz,
      nm=nm,
      rhom=rhom,
      rhop=rhop,
    	rhoq=rhoq,
      np=np,
      nq=nq,
      Gamma=Gamma,
      wt.pc=wt.pc,
      Dz=Dz,
      rhoz=rhoz,
      nz=nz,
      nb=nb,
      rhob=rhob,
      b0=b0,
      z0=z0,
      T0=T0,
      zeta=zeta,
    	lambda.probe=lambda.probe,
    	Delta.n=2*sqrt(N1a*N1a + N1b*N1b),
      Moharam.Young=lambda.probe*lambda.probe/mean_RI/Delta.n/Lambda.t/Lambda.t/cos(phi.r.t),
    	Klein.Cook=2*pi*lambda.probe*T.t/mean_RI/Lambda.t/Lambda.t/cos(phi.r.t),
      Delta.t=Delta.t,
      Delta.x=Delta.x,
      output.time.step=output.time.step,
      Model="slanted_grating_simulation_v21"
    ) %>% 
    arrange(time,Y,x) -> df.9
  
  if(
    nrow(
      subset(
        df.9, 
        mass.z < 0
      )
    ) > 0
  ){
    return(warning("Fail: z(x,t) < 0"))
  } else {
    return(df.9)
  }
  
}
