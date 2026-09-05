from warnings import filterwarnings
filterwarnings("ignore")
import numpy as np
import pandas as pd
from math import pi, factorial
from numpy.linalg import inv
from time import time as gettime
from numerical_integration import simpsons_rule_1D, trapezoidal_rule_integration
from jl_errors import non_negative_args

class UnslantedBraggGrating:
    
    @non_negative_args
    def __init__(
            self, 
            start_exp=0,# Start of exposure
            end_exp=1e2,# End of exposure
            total_time=1e2,# Total simulation time
            lpmm=1e3,# Spatial frequency
            I0=5,# Intensity of recording beam
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
            Dz=1e-10,# Nanoparticle self-diffusion coefficient
            epsilon_pz=13,# Cross-diffusion ratio
            epsilon_qz=13,# Cross-diffusion ratio
            rhoz=1.74,# Nanoparticle mass density
            n_z=1.366,# Nanoparticle refractive index
            b0=5.05,# Ratio of binder to monomer mass
            n_b=1.5,# Binder refractive index
            rhob=1.19,# Binder mass density
            lambda_probe=633e-7,# Wavelength of reconstruction beam
            Delta_t=1/100,# Numerical scheme time step
            Delta_x=1/20,# Numerical scheme spatial step
            output_time_step=1# Seconds
        ):
        
        self.__total_time = total_time
        self.__end_exp = end_exp
        self.__lpmm = lpmm
        self.__I0 = I0
        self.__xi = xi
        self.__Dm = Dm
        self.__n_m = n_m
        self.__rhom = rhom
        self.__Dp = Dp
        self.__rhop = rhop
        self.__n_p = n_p
        self.__n_q = n_q
        self.__Gamma = Gamma
        self.__Dz = Dz
        self.__epsilon_pz = epsilon_pz
        self.__epsilon_qz = epsilon_qz
        self.__wt_pc = wt_pc
        self.__rhoz = rhoz
        self.__n_z = n_z
        self.__b0 = b0
        self.__n_b = n_b
        self.__rhob = rhob
        self.__lambda_probe = lambda_probe
        self.__Delta_x = Delta_x
        try:
            self.__Nx = int(1/Delta_x) + 1
        except ZeroDivisionError:
            print("Error: Delta_x > 0")
            return None
        self.__Delta_t = Delta_t
        self.__output_time_step = output_time_step
        
    def __str__(self):
        
        table_width = 107
        left_column_width = int((table_width-7)*.8)
        right_column_width = int((table_width-7)*.2)
        hl = "-"*table_width + "\n"
        table_header = "| Holographic Grating".ljust(table_width-1) + "|" + "\n"
        
        table_subheader_1 = "| Properties of the Host Photopolymer".ljust(table_width-1) + "|" + "\n"
        host_properties = ['Diffusion Coefficient (Monomer) [cm2/s]','Refractive Index (Monomer)', 
                           "Refractive Index (Polymer)",'Diffusion Coefficient (Polymer) [cm2/s]', "Immobilization Constant [/s]",
                          'Density (Monomer) [g/cm3]', 'Density (Polymer) [g/cm3]','Binder to Monomer Ratio','Refractive Index (Binder)','Density (Binder) [g/cm3]']
        host_values = [self.__Dm, self.__n_m, self.__n_q, self.__Dp,self.__Gamma,self.__rhom,self.__rhop,self.__b0,self.__n_b,self.__rhob]
        table_body_photopolymer = ""
        for i, j in zip(host_properties, host_values):
            table_body_photopolymer += "| " + i.ljust(left_column_width) + " | " + str(j).ljust(right_column_width) + " |" + "\n"
        
        table_subheader_2 = "| Recording Conditions".ljust(table_width-1) + "|" + "\n"
        recording_conditions = ['Recording Duration [s]','Spatial Frequency [lines/mm]','Recording Intensity [mW/cm2]','Wavelength of the Reconstruction Beam [cm]']
        recording_values = [self.__total_time,self.__lpmm,self.__I0,self.__lambda_probe]
        table_body_recording = ""
        for i, j in zip(recording_conditions, recording_values):
            table_body_recording += "| " + i.ljust(left_column_width) + " | " + str(j).ljust(right_column_width) + " |" + "\n"
            
        table_subheader_3 = "| Nanoparticle Properties".ljust(table_width-1) + "|" + "\n"
        nanoparticle_properties = ['Refractive Index (Nanoparticles)', 'Density (Nanoparticles) [g/cm3]','Diffusion Coefficient (Nanoparticles) [cm2/s]',
                                   'Oligomer-Nanoparticle Cross-Diffusion Constant','Polymer-Nanoparticle Cross-Diffusion Constant','Scattering']
        nanoparticle_values = [self.__n_z,self.__rhoz,self.__Dz,self.__epsilon_pz,self.__epsilon_qz,self.__xi]
        table_body_nanoparticle = ""
        for i, j in zip(nanoparticle_properties, nanoparticle_values):
            table_body_nanoparticle += "| " + i.ljust(left_column_width) + " | " + str(j).ljust(right_column_width) + " |" + "\n"
            
        return hl + table_header + hl + table_subheader_1 + hl+ table_body_photopolymer + hl +table_subheader_2 + hl +table_body_recording + hl + table_subheader_3 + hl + table_body_nanoparticle + hl
        
    def get_host_photopolymer_properties(self):
        
        return {'Dm':self.__Dm, 'n_m':self.__n_m, 'n_q':self.__n_q, 'Gamma':self.__Gamma, 'T0':self.__T0, 
                'xi':self.__xi,'rhom':self.__rhom,'rhop':self.__rhop,'b0':self.__b0,'n_b':self.__n_b}
    
    def get_recording_conditions(self):
        
        return {'total_time':self.__total_time,'lpmm':self.__lpmm,'I0':self.__I0,
                'slant_angle':self.__slant_angle,'lambda_probe':self.__lambda_probe}
    
    def get_nanoparticle_properties(self):
        
        return {'n_z':self.__n_z,'rhoz':self.__rhoz,'Dz':self.__Dz,'epsilon_pz':self.__epsilon_pz,'epsilon_qz':self.__epsilon_qz}
    
    def holographic_recording(self):
        
        start_computation = gettime()
        # 1.2 --- Define parameters
        Nx=int(1/self.__Delta_x) + 1# Number of spatial points
        if Nx%2==0:
            raise Exception("Number of x mesh points must be an odd number.")
            return None
        
        x=np.linspace(0,1,Nx)# Non-dimensional grating distance
        try:
            n_iterations = int(self.__total_time/self.__Delta_t)+1# Total number of iterations
        except ZeroDivisionError:
            print("Error: Delta_t > 0")
            return None
        
        r=self.__Delta_t/self.__Delta_x/self.__Delta_x# Ratio of finite time step to squared finite spatial step
        m0=1# Initial mass of monomer
        t0=1 #  Reference time [s]
        try:
            Lambda=1/10/self.__lpmm # Grating period [cm]
        except ZeroDivisionError:
            print("Error: lpmm > 0")
            return None
        j_end_exp = self.__end_exp/self.__Delta_t # Iteration of exposure end
        try:
            z0 = self.__wt_pc/(1 - self.__wt_pc)*(m0 + self.__b0)# Initial nanoparticle to monomer
        except ZeroDivisionError:
            print("Error: 0 <= wt_pc < 1")
            return None
        
        # 1.3 --- Matrix initial conditions
        m1=np.ones(Nx); m2=m1
        p1=np.zeros(Nx); p2=p1
        q1=np.zeros(Nx); q2=q1
        z1=np.ones(Nx); z2=z1
        b1=np.ones(Nx)
        
        Vb = self.__b0*b1/self.__rhob # cm**3
        Vm = m1*m0/self.__rhom # cm**3
        Vp = p1*m0/self.__rhop # cm**3
        Vq = q1*m0/self.__rhop # cm**3
        Vz = z1*z0/self.__rhoz # cm**3
        Vtotal=Vb+Vm+Vp+Vq+Vz # cm**3
        
        phi_m = Vm/Vtotal
        phi_b = Vb/Vtotal
        phi_p = Vp/Vtotal
        phi_q = Vq/Vtotal
        phi_z = Vz/Vtotal
        
        Lorentz_Lorenz_RHS = phi_m*(self.__n_m*self.__n_m - 1)/(self.__n_m*self.__n_m + 2) + phi_b*(self.__n_b*self.__n_b - 1)/(self.__n_b*self.__n_b + 2) + phi_p*(self.__n_p*self.__n_p - 1)/(self.__n_p*self.__n_p + 2) + phi_q*(self.__n_q*self.__n_q - 1)/(self.__n_q*self.__n_q + 2) + phi_z*(self.__n_z*self.__n_z - 1)/(self.__n_z*self.__n_z + 2)
        
        n1=np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
        
        interior_points = list(range(1, Nx-1))
        
        times_4 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2 != 0]
        
        times_2 = [interior_points[i] for i in range(len(interior_points)) if interior_points[i]%2==0]
        
        N0=[self.__Delta_x/3*(n1[0] + sum(2*n1[times_2]) + sum(4*n1[times_4]) + n1[Nx-1])]
        
        n1_cos=n1*np.cos(2*pi*x)
        
        N1=[2*self.__Delta_x/3*(n1_cos[0] + sum(2*n1_cos[times_2]) + sum(4*n1_cos[times_4]) + n1_cos[Nx-1])]
            
        sq_diff=(n1 - N0[0] - N1[0]*np.cos(2*pi*x))**2
        
        d2 = [self.__Delta_x/3*(sq_diff[0] + sum(2*sq_diff[times_2]) + sum(4*sq_diff[times_4]) + sq_diff[Nx-1])]
        
        # 1.2 --- Non-dimensional parameters
        alpha_m=self.__Dm*t0/Lambda/Lambda# Monomer diffusion
        F0=0.1*self.__I0**0.3
        beta=F0*t0# Monomer comsumption
        alpha_p=self.__Dp*t0/Lambda/Lambda# Oligomer diffusion
        alpha_z=self.__Dz*t0/Lambda/Lambda# Nondimensional self-nanoparticle diffusion
        gamma = m0*self.__Gamma*t0# Immobilization
        if self.__wt_pc==0:
            alpha_pz=0
            alpha_zp=0
            alpha_zq=0
        else:
            alpha_pz=z0*self.__epsilon_pz*alpha_p
            alpha_zp=self.__epsilon_pz*alpha_z
            alpha_zq=self.__epsilon_qz*alpha_z
        
        spatial_profile_DF=pd.DataFrame({"x": x, 
                         "monomer": m1,
                         "short_polymer": p1,
                         "immobile_polymer": q1,
                         'nanoparticles': z0*z1,
                         'binder': self.__b0*b1,
                         'refractive_index': n1,
                         "time": np.zeros(Nx)})
        
        time_vals=np.arange(0, self.__total_time+1, self.__output_time_step)
        
        # 1.3 --- Calculate each time step via implicit finite difference method
        for j in range(1,n_iterations):
            # Phi=1 if illumination is on, 0 otherwise
            if ((j >= 0) and (j <= j_end_exp)):
                Phi=1
            else:
                Phi = 0
                
            # Illumination pattern
            f = 1 + np.exp(-self.__xi*z0*z1)*np.cos(2*pi*x)
            MM2 = (2 + Phi*self.__Delta_t*beta*f)*np.identity(Nx)
            MM1 = (2 - Phi*self.__Delta_t*beta*f)*np.identity(Nx)
            PP2 = (2 + Phi*self.__Delta_t*gamma*p1)*np.identity(Nx)
            PP1 = (2 - Phi*self.__Delta_t*gamma*p1)*np.identity(Nx)
            PM2 = Phi*self.__Delta_t*beta*f*np.identity(Nx)
            PM1 = Phi*self.__Delta_t*beta*f*np.identity(Nx)
            QQ2=2*np.identity(Nx)
            QQ1=2*np.identity(Nx)
            QP2=(Phi*gamma*self.__Delta_t*p1)*np.identity(Nx)
            QP1=(Phi*gamma*self.__Delta_t*p1)*np.identity(Nx)
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
                
                MM2[i, i_minus_1] = MM2[i, i_minus_1] - r*alpha_m
                MM2[i, i] = MM2[i, i] + 2*r*alpha_m
                MM2[i, i_plus_1] = MM2[i, i_plus_1] - r*alpha_m
        
                MM1[i, i_minus_1] = MM1[i, i_minus_1] + r*alpha_m
                MM1[i, i] = MM1[i, i] - 2*r*alpha_m
                MM1[i, i_plus_1] = MM1[i, i_plus_1] + r*alpha_m
        
                PP2[i, i_minus_1] = PP2[i, i_minus_1] - r*alpha_p -  r*alpha_pz*z1[i_minus_1]
                PP2[i, i] = PP2[i, i] + 2*r*alpha_p + 2*r*alpha_pz*z1[i]
                PP2[i, i_plus_1] = PP2[i, i_plus_1] - r*alpha_p - r*alpha_pz*z1[i_plus_1]
                
                PP1[i, i_minus_1] = PP1[i, i_minus_1] + r*alpha_p + r*alpha_pz*z1[i_minus_1]
                PP1[i, i] = PP1[i, i] - 2*r*alpha_p - 2*r*alpha_pz*z1[i]
                PP1[i, i_plus_1] = PP1[i, i_plus_1] + r*alpha_p + r*alpha_pz*z1[i_plus_1]
        
                ZZ2[i, i_minus_1] = ZZ2[i, i_minus_1] - r*alpha_z - r*alpha_zq*q1[i_minus_1] - r*alpha_zp*p1[i_minus_1]
                ZZ2[i, i] = ZZ2[i, i] + 2*r*alpha_z + 2*r*alpha_zq*q1[i] + 2*r*alpha_zp*p1[i]
                ZZ2[i, i_plus_1] = ZZ2[i, i_plus_1] - r*alpha_z - r*alpha_zq*q1[i_plus_1]  - r*alpha_zp*p1[i_plus_1]
        
                ZZ1[i, i_minus_1] = ZZ1[i, i_minus_1] + r*alpha_z + r*alpha_zq*q1[i_minus_1] + r*alpha_zp*p1[i_minus_1]
                ZZ1[i, i] = ZZ1[i, i] - 2*r*alpha_z - 2*r*alpha_zq*q1[i]  - 2*r*alpha_zp*p1[i]
                ZZ1[i, i_plus_1] = ZZ1[i, i_plus_1] + r*alpha_z + r*alpha_zq*q1[i_plus_1] + r*alpha_zp*p1[i_plus_1]
                
            m2 = np.matmul(inv(MM2), np.matmul(MM1,m1))
            
            p2 = np.matmul(inv(PP2), (np.matmul(PP1,p1) + np.matmul(PM2,m2) + np.matmul(PM1,m1)))
            
            q2 = np.matmul(inv(QQ2), (np.matmul(QQ1,q1) + np.matmul(QP2,p2) + np.matmul(QP1,p1)))
            
            z2 = np.matmul(inv(ZZ2), np.matmul(ZZ1,z1))
            
            m1=m2; p1=p2; q1=q2; z1=z2
            
            Vb = self.__b0*b1/self.__rhob # cm**3
            Vm = m1*m0/self.__rhom # cm**3
            Vp = p1*m0/self.__rhop # cm**3
            Vq = q1*m0/self.__rhop # cm**3
            Vz = z1*z0/self.__rhoz # cm**3
            Vtotal=Vb+Vm+Vp+Vq+Vz # cm**3
            phi_m = Vm/Vtotal
            phi_b = Vb/Vtotal
            phi_p = Vp/Vtotal
            phi_q = Vq/Vtotal
            phi_z = Vz/Vtotal
        
            Lorentz_Lorenz_RHS = phi_m*(self.__n_m*self.__n_m - 1)/(self.__n_m*self.__n_m + 2) + phi_b*(self.__n_b*self.__n_b - 1)/(self.__n_b*self.__n_b + 2) + phi_p*(self.__n_p*self.__n_p - 1)/(self.__n_p*self.__n_p + 2) + phi_q*(self.__n_q*self.__n_q - 1)/(self.__n_q*self.__n_q + 2) + phi_z*(self.__n_z*self.__n_z - 1)/(self.__n_z*self.__n_z + 2)
        
            n1=np.sqrt((2*Lorentz_Lorenz_RHS + 1)/(1 - Lorentz_Lorenz_RHS))
        
            N0_new=self.__Delta_x/3*(n1[0] + sum(2*n1[times_2]) + sum(4*n1[times_4]) + n1[Nx-1])
            
            n1_cos=n1*np.cos(2*pi*x)
        
            N1_new=2*self.__Delta_x/3*(n1_cos[0] + sum(2*n1_cos[times_2]) + sum(4*n1_cos[times_4]) + n1_cos[Nx-1])
            
            sq_diff=(n1 - N0_new - N1_new*np.cos(2*pi*x))**2
        
            d2_new = [self.__Delta_x/3*(sq_diff[0] + sum(2*sq_diff[times_2]) + sum(4*sq_diff[times_4]) + sq_diff[Nx-1])]
            
            if self.__Delta_t*j in time_vals:
                
                new_df=pd.DataFrame({"x": x, 
                                 "monomer": m1,
                                 "short_polymer": p1,
                                 "immobile_polymer": q1,
                                 'nanoparticles': z0*z1,
                                 'binder': self.__b0*b1,
                                 'refractive_index': n1,
                                 "time": np.zeros(Nx) + self.__Delta_t*j})
                
                spatial_profile_DF=pd.concat([spatial_profile_DF,new_df])
                
                N0.append(N0_new)
                
                N1.append(N1_new)
                
                d2.append(d2_new)
        
        optical_properties_DF=pd.DataFrame({
            'time': time_vals, 
            'N0': N0, 
            'Delta_n': [2*i for i in N1]})         
    
        self.spatial_profile_DF = spatial_profile_DF
        self.optical_properties_DF = optical_properties_DF

        end_computation = gettime()
        self.__computation_time = end_computation-start_computation
        
        table_width=107
        print("-"*table_width)
        print(" Numerical Simulation Complete ".center(table_width, '-'))
        print("-"*table_width)
        print("| Computation Time: {:.1f} s".format(self.__computation_time).ljust(table_width-1)+'|')        
        print(UnslantedBraggGrating.__str__(self))
        
    def get_simulation_properties(self):
        
        return {'Delta_t':self.__Delta_t, 'Delta_x':self.__Delta_x, 'Nx':self.__Nx,'computation_time':self.__computation_time}

if __name__ == "__main__":
    a = UnslantedBraggGrating(total_time=0)
    a.holographic_recording()
    assert len(a.optical_properties_DF) == 1
    assert len(a.spatial_profile_DF) == a.get_simulation_properties()['Nx']
    del a
    print("All tests passed\a.")

