# coupled_cross_diffusion_v7.R
# Author: Jack Lyons

# A function to run a numerical simulation of a mathematical model to describe the formation of unslanted holographic gratings.

# 0.1 --- Import libraries
library(dplyr)
library(reshape2)
library(ggplot2)

# 0.2 --- Increase memory limit
memory.limit(size=1e6)

# 0.3 --- Set working directory
setwd("~/PhD Documents/Inert Zeolite Nanoparticles Model/Unslanted Gratings/Coupled Cross-Diffusion/")

# 0.4 --- Import plot aesthetics
load("~/PhD Documents/plt.theme1.RData")
load("~/PhD Documents/plt.theme2.RData")

# 1.1 --- Define function inputs and default values.
coupled_cross_diffusion_v7 = function(
    
  end_exp=1e2, # End of exposure
  total.time=100,# Total simulation time
  lpmm=1e3,# Spatial frequency
	T0=50e-4,# cm
  F0=0.1*5**0.3,# Intensity of first recording beam
  xi=0.3,# Scattering coefficient
  Dm=1.6e-7,# Monomer diffusion coefficient
  nm=1.55,# Monomer refractive index
  rhom=1.15,# Monomer mass density
  Dp=6.35e-10,# Oligomer diffusivity ratio
  rhop=1.3,# Fractional van der Waals space loss
  np=1.56,# Oligomer refractive index
  nq=1.64,# Polymer refractive index
  Gamma=1,# Rate of immobilization
  Dz=1e-10,# Nanoparticle cross-diffusion ratio
	epsilon.mz=0,# Cross-diffusion
	epsilon.pz=13,# Cross-diffusion
	epsilon.qz=13,# Cross-diffusion
  wt.pc=5e-2,# Doping %
  rhoz=1.74,# Nanoparticle mass density
  nz=1.366,# Nanoparticle refractive index
  b0=5.05,# Binder mass
	nb=1.5,# Binder refractive inedx
  rhob=1.19, # Binder mass density
  Delta.t=1/100,# Numerical scheme time step
  Delta.x=1/20,# Numerical scheme spatial step
  output.time.step=1# Output time
  
){
  
  # 1.2 --- Define parameters
	iterations.per.second=1/Delta.t# Number of iterations each second
  Nx=1/Delta.x + 1# Number of spatial points
  x=seq(0, 1, length.out=Nx)# Non-dimensional grating distance
  Nt = total.time*iterations.per.second +1# Total number of iterations
  r=Delta.t/Delta.x/Delta.x# Ratio of finite time step to squared finite spatial step
  m0=1# Initial grams of monomer
  t0=1 #  Reference time [s]
  Lambda=1/10/lpmm # Grating period [cm]
  j_start_exp=0/Delta.t # Iteration of exposure start
  j_end_exp=end_exp/Delta.t # Iteration of exposure end
  z0 = wt.pc/(1 - wt.pc)*(m0 + b0)# Initial nanoparticle mass
  rep(1, Nx) -> m1# m at j=0
  rep(0, Nx) -> p1# p at j=0
  rep(0, Nx) -> q1# q at j=0
  rep(1, Nx) -> z1# z at j=0
  
  # 1.3 --- Non-dimensional parameters
  alpha.m=Dm*t0/Lambda/Lambda# Monomer diffusion
  beta=F0*t0# Monomer comsumption
  alpha.p=Dp*t0/Lambda/Lambda# Oligomer diffusion
  alpha.z=Dz*t0/Lambda/Lambda# Nondimensional self-nanoparticle diffusion
  gamma = m0*Gamma*t0# Immobilization
  alpha.mz=ifelse(wt.pc==0, 0, z0*epsilon.mz*alpha.m)
  alpha.pz=ifelse(wt.pc==0, 0, z0*epsilon.pz*alpha.p)
  alpha.zm=ifelse(wt.pc==0, 0, epsilon.mz*alpha.z)
  alpha.zp=ifelse(wt.pc==0, 0, epsilon.pz*alpha.z)
  alpha.zq=ifelse(wt.pc==0, 0, epsilon.qz*alpha.z)
  
  matrix(m1, ncol=1) -> m# monomer density matrix
  matrix(p1, ncol=1) -> p# oligomer density matrix
  matrix(q1, ncol=1) -> q# Immobile polymer matrix
  matrix(z1, ncol=1) -> z# Zeolite matrix
  
  seq(0, total.time, by=output.time.step) -> time.vals
  time.vals/Delta.t + 1 -> bind.iter
  bind.iter = as.character(bind.iter[-1])
  
  # 1.3 --- Calculate each time step via implicit finite difference method
  for(j in 1:Nt){
    
    if(j >= j_start_exp & j <= j_end_exp){Phi=1} else {Phi = 0}# Phi=1 if illumination is on, 0 otherwise
    
    f = 1 + exp(-xi*z0*z1)*cos(2*pi*x)# Illumination pattern
    f = as.numeric(f)
    
    MM2 = diag(2 + Phi*Delta.t*beta*f, Nx)
    MM1 = diag(2 - Phi*Delta.t*beta*f, Nx)
    PP2 = diag(2 + Phi*Delta.t*gamma*as.numeric(p1), Nx)
    PP1 = diag(2 - Phi*Delta.t*gamma*as.numeric(p1), Nx)
    PM2 = diag(+Phi*Delta.t*beta*f, Nx)
    PM1 = diag(+Phi*Delta.t*beta*f, Nx)
    QQ2=diag(2, Nx)
    QQ1=diag(2, Nx)
    QP2=diag(+Phi*gamma*Delta.t*as.numeric(p1),Nx)
    QP1=diag(+Phi*gamma*Delta.t*as.numeric(p1),Nx)
    ZZ2 = diag(2, Nx)
    ZZ1 = diag(2, Nx)
    
    for(i in 1:(Nx)){
      
      i.minus.1=ifelse(i == 1, i + 1, i - 1)
      i.plus.1=ifelse(i == Nx, i - 1, i + 1)
      i.minus.2=ifelse(i == 1, i + 2,ifelse(i == 2, i, i - 2))
      i.plus.2=ifelse(i == Nx, i - 2,ifelse(i == Nx-1, i, i + 2))
      
      MM2[i, i.minus.1] = MM2[i, i.minus.1] - r*alpha.m - r*alpha.mz*z1[i.minus.1]
	    MM2[i, i] = MM2[i, i] + 2*r*alpha.m + 2*r*alpha.mz*z1[i]
	    MM2[i, i.plus.1] = MM2[i, i.plus.1] - r*alpha.m - r*alpha.mz*z1[i.plus.1]
	    
	    MM1[i, i.minus.1] = MM1[i, i.minus.1] + r*alpha.m + r*alpha.mz*z1[i.minus.1]
	    MM1[i, i] = MM1[i, i] - 2*r*alpha.m - 2*r*alpha.mz*z1[i]
	    MM1[i, i.plus.1] = MM1[i, i.plus.1] + r*alpha.m + r*alpha.mz*z1[i.plus.1]
	    
      PP2[i, i.minus.1] = PP2[i, i.minus.1] - r*alpha.p -  r*alpha.pz*z1[i.minus.1]
	    PP2[i, i] = PP2[i, i] + 2*r*alpha.p + 2*r*alpha.pz*z1[i]
	    PP2[i, i.plus.1] = PP2[i, i.plus.1] - r*alpha.p - r*alpha.pz*z1[i.plus.1]
	    
	    PP1[i, i.minus.1] = PP1[i, i.minus.1] + r*alpha.p + r*alpha.pz*z1[i.minus.1]
	    PP1[i, i] = PP1[i, i] - 2*r*alpha.p - 2*r*alpha.pz*z1[i]
	    PP1[i, i.plus.1] = PP1[i, i.plus.1] + r*alpha.p + r*alpha.pz*z1[i.plus.1]
	    
      ZZ2[i, i.minus.1] = ZZ2[i, i.minus.1] - r*alpha.z - r*alpha.zq*q1[i.minus.1] - r*alpha.zp*p1[i.minus.1]  - r*alpha.zm*m1[i.minus.1]
	    ZZ2[i, i] = ZZ2[i, i] + 2*r*alpha.z + 2*r*alpha.zq*q1[i] + 2*r*alpha.zp*p1[i] + 2*r*alpha.zm*m1[i]
	    ZZ2[i, i.plus.1] = ZZ2[i, i.plus.1] - r*alpha.z - r*alpha.zq*q1[i.plus.1]  - r*alpha.zp*p1[i.plus.1]  - r*alpha.zm*m1[i.plus.1]
	    
	    ZZ1[i, i.minus.1] = ZZ1[i, i.minus.1] + r*alpha.z + r*alpha.zq*q1[i.minus.1] + r*alpha.zp*p1[i.minus.1] + r*alpha.zm*m1[i.minus.1]
	    ZZ1[i, i] = ZZ1[i, i] - 2*r*alpha.z - 2*r*alpha.zq*q1[i]  - 2*r*alpha.zp*p1[i]  - 2*r*alpha.zm*m1[i]
	    ZZ1[i, i.plus.1] = ZZ1[i, i.plus.1] + r*alpha.z + r*alpha.zq*q1[i.plus.1] + r*alpha.zp*p1[i.plus.1] + r*alpha.zm*m1[i.plus.1]
	    
    }
    
    m2 = solve(MM2) %*% MM1 %*% matrix(m1, ncol=1)
    
    p2 = (solve(PP2) %*% PP1 %*% matrix(p1, ncol=1)) + (solve(PP2) %*% PM2 %*% matrix(m2, ncol=1)) + (solve(PP2) %*% PM1 %*% matrix(m1, ncol=1))
    
    q2 = (solve(QQ2) %*% QQ1 %*% matrix(q1, ncol=1)) + (solve(QQ2) %*% QP2 %*% matrix(p2, ncol=1)) + (solve(QQ2) %*% QP1 %*% matrix(p1, ncol=1))
    
    z2 = solve(ZZ2) %*% ZZ1 %*% matrix(z1, ncol=1)
    
    if(as.character(j) %in% bind.iter){
      
      matrix(cbind(m, m2), nrow=Nx) -> m
      matrix(cbind(p, p2), nrow=Nx) -> p
      matrix(cbind(q, q2), nrow=Nx) -> q
      matrix(cbind(z, z2), nrow=Nx) -> z
      
    }
    
    m2 -> m1
    p2 -> p1
    q2 -> q1
    z2 -> z1
    
  }
  
  # 2.1 --- Monomer density data frame
  data.frame(x=seq(0,1,Delta.x),m) -> df.m# Turn the matrix into a data frame.
  names(df.m)=c("x", paste0("t",time.vals[1:ncol(m)]))# Rename columns
  df.m %>% melt(id.vars=c("x")) %>% rename(time=variable, mass.m=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) -> melt.df.m# Additional columns containing info on input
  
  # 2.2 --- Oligomer density data frame
  data.frame(x=seq(0,1,Delta.x), p) -> df.p# Turn the matrix into a data frame.
  names(df.p)=c("x", paste0("t",time.vals[1:ncol(p)]))# Rename columns
  df.p %>% melt(id.vars=c("x")) %>% rename(time=variable, mass.p=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) -> melt.df.p# Additional columns containing info on input
  
  # 2.3 --- Polymer density data frame
  data.frame(x=seq(0,1,Delta.x), q) -> df.q# Turn the matrix into a data frame.
  names(df.q)=c("x", paste0("t",time.vals[1:ncol(q)]))# Rename columns
  df.q %>% melt(id.vars=c("x")) %>% rename(time=variable, mass.q=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) -> melt.df.q# Additional columns containing info on input
  
  # 2.4 --- Zeolite density data frame
  data.frame(x=seq(0,1,Delta.x), z) -> df.z# Turn the matrix into a data frame.
  names(df.z)=c("x", paste0("t",time.vals[1:ncol(z)]))# Rename columns
  df.z %>% melt(id.vars=c("x")) %>% rename(time=variable, mass.z=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time),mass.z=z0*mass.z) -> melt.df.z# Additional columns containing info on input
  
  # 2.5 --- Element density data frame
  merge(melt.df.m, melt.df.p) -> df
  merge(df, melt.df.q) -> df
  merge(df, melt.df.z) -> df
  subset(df, time %in% time.vals) -> df
  
  # 3.1 --- Refractive index matrix
  matrix(1, nrow=nrow(m), ncol=ncol(m)) -> b# Binder density matrix
  
  Vb = b*b0/rhob # cm**3
  Vm = m*m0/rhom # cm**3
  Vp = p*m0/rhop # cm**3
  Vq = q*m0/rhop # cm**3
  Vz = z*z0/rhoz # cm**3
  Vtotal=Vb+Vm+Vp+Vq+Vz # cm**3
  
  phi.m = Vm/Vtotal
  phi.b = Vb/Vtotal
  phi.p = Vp/Vtotal
  phi.q = Vq/Vtotal
  phi.z = Vz/Vtotal
  
  Lorentz.Lorenz.RHS = phi.m*(nm*nm - 1)/(nm*nm + 2) + phi.b*(nb*nb - 1)/(nb*nb + 2) + phi.p*(np*np - 1)/(np*np + 2) + phi.q*(nq*nq - 1)/(nq*nq + 2) + phi.z*(nz*nz - 1)/(nz*nz + 2)
  sqrt((2*Lorentz.Lorenz.RHS + 1)/(1 - Lorentz.Lorenz.RHS)) -> n
  
  # 3.2 --- Refractive index data frame
  data.frame(x=seq(0,1,Delta.x), n) -> df.n# Turn the matrix into a data frame.
  names(df.n)=c("x", paste0("t",time.vals[1:ncol(n)]))# Rename columns
  
  df.n %>% melt(id.vars=c("x")) %>% rename(time=variable, n=value) %>% mutate(time=gsub(x=time, pattern="t", replacement=""), time=as.numeric(time)) -> df.n
  df=merge(df, df.n, by=c("x","time"))
  
  # 4.1 --- Refractive index modulation data frame
  N = ncol(n)
  N0=rep(0, N)
  N1=rep(0, N)
  n.tilde=matrix(0, nrow=nrow(n), ncol=ncol(n))
  diff=matrix(0, nrow=nrow(n), ncol=ncol(n))
  n2=matrix(0, nrow=nrow(n), ncol=ncol(n))
  d2=rep(0, N)
  int.n2=rep(0,N)
  h=matrix(0, nrow=Nx, ncol=N)
  h0=m0/rhom + z0/rhoz + b0/rhob
  shrinkage=rep(0,N)
  
  for(k in 1:N){# For each time point ...
    
    interior.points = 2:(Nx-1)
    times.4 = interior.points[interior.points %% 2 == 0]
    times.2 = interior.points[interior.points %% 2 != 0]
    
    n0=n[, k]
    n1=n[, k]*cos(2*pi*x)
    
    N0[k]=Delta.x/3*(n0[1] + n0[Nx] + sum(2*n0[times.2]) + sum(4*n0[times.4]))
    N1[k]=2*Delta.x/3*(n1[1] + n1[Nx] + sum(2*n1[times.2]) + sum(4*n1[times.4]))
    
    n.tilde[, k] = N0[k] + N1[k]*cos(2*pi*x)
    diff[, k] = (n[, k] - n.tilde[, k])**2
    d2[k] = Delta.x/3*(diff[1, k] + diff[Nx, k] + sum(2*diff[times.2, k]) + sum(4*diff[times.4, k]))
    n2[, k]=n[, k]**2
    int.n2[k]=Delta.x/3*(n2[1, k] + sum(2*n2[times.2, k]) + sum(4*n2[times.4, k]) + n2[Nx, k])
    h[,k]=(m[,k]/rhom + p[,k]/rhop + q[,k]/rhop + z[,k]*z0/rhoz + b0/rhob)/h0
    int.h.dx=Delta.x/3*(h[1,k] + sum(4*h[times.4,k]) + sum(2*h[times.2,k]) + h[Nx,k])
    shrinkage[k] = 1 - int.h.dx
    
  }# End for loop
  
  Delta.n=2*abs(N1)
  distortion=d2/int.n2
  time=time.vals[1:ncol(n)]# Time vector
  redistributed=c()
  for(i in time){
    df %>% subset(time==i) %>% arrange(x) %>% mutate(v1=mass.z/z0) %>% pull(v1) -> zi
    ifelse(zi > 1, zi - 1, 0) -> zi.gt.1
    Delta.x/3*(zi.gt.1[1] + sum(2*zi.gt.1[times.2]) + sum(4*zi.gt.1[times.4]) + zi.gt.1[Nx]) -> v
    c(redistributed,v) -> redistributed
  }
  df.Dn=data.frame(time, shrinkage, distortion, Delta.n, redistributed, N0)
  df=merge(df, df.Dn, by=c("time"))
  
  df %>% mutate(end_exp=end_exp,total.time,lpmm,F0,xi,Dm,nm,rhom,rhop,Dp,np,nq,Gamma,wt.pc,Dz,epsilon.mz,epsilon.pz,epsilon.qz,rhoz,nz,nb,rhob,b0,z0,T0,T.t=T0*(1-shrinkage),alpha.m,beta,Delta.t,Delta.x,output.time.step,Model="coupled_cross_diffusion_v7") -> df
  
  if(
    nrow(
      subset(
        df, 
        mass.z < 0
      )
    ) > 0
  ){
    return(cat("Fail: z(x,t) < 0"))
  } else {
    return(df)
  }
  
}
