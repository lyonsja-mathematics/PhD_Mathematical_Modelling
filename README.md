This repository contains the final version of the programming scripts I wrote for my PhD thesis on *Mathematical Modelling of Hybrid Photonic Structures for Holographic Sensors*. My research was funded by the Irish Research Council (GOIPG/2021/214).

![The formation of holographic diffraction gratings in hybrid photopolymer media](https://github.com/lyonsja-mathematics/PhD_Mathematical_Modelling/blob/main/Holographic%20Photopolymerization%20v3.png)

# Motivation & Research Questions

The motivation for this research project was as follows:

- Optical properties of holographic gratings are sensitive to external stimulus and can be exploited for the purpose of environmental sensing.
- Hybrid photopolymers are a strong candidate as recording media for holographic gratings because they offer a wide dynamic range and good selectivity.
- A mathematical framework modelling the formation and operation of a holographic sensor is needed to optimize the design.

The successful completion of this research was characterized by the ability to answer the following research questions.

- How can we use the ...
1. host photopolymer material properties (monomer, binder matrix, dye, etc.)
2. recording conditions (recording intensity, spatial frequency, etc.)
3. nanoparticle properties (initial concentration, refractive index, etc.)

in order to control the final grating in a hybrid polymer system and hence optimize the functionality of a holographic sensor? 
		
- Can the theoretical model predict IEO experimental results:
1. Significantly increased dynamic range over conventional photopolymer media.
2. Non-linear response of refractive index modulation to increased doping.
3. Photopolymerization-induced shrinkage is significantly reduced by the addition of zeolite nanoparticles.
4. Increased shrinkage at high spatial frequencies.

# Mathematical Model

![Optical setup and geometry for a slanted holographic diffraction grating](https://github.com/lyonsja-mathematics/PhD_Mathematical_Modelling/blob/main/holographic_grating_geometry_setup.png)

The hybrid photopolymer system consists of a binder matrix ($b$) such as polyvynyl alcohol (PVA), a monomer ($m$) such as acrylamide (AA), an electron donor, a dye, a crosslinking monomer and inorganic nanoparticles ($z$) such as zeolites. The two-way diffusion model treats the monomer, dye, electron donor and crosslinker as a single constituent. The rates of change of each constituent in the photopolymer system is modelled by the following coupled system of partial differential equations,

```math
\begin{align}
\frac{\partial b}{\partial t} &= 0,\\
\frac{\partial m}{\partial t} &+ \nabla \cdot \vec{J_m} =  - \Phi(t) F(x,y,t) m,\\
\frac{\partial p}{\partial t} &+ \nabla \cdot \vec{J_p} = \Phi(t) F(x,y,t) m - \Phi(t) \Gamma p^2,\\
\frac{\partial q}{\partial t} &= \Phi(t) \Gamma p^2,\\
\frac{\partial z}{\partial t} &+ \nabla \cdot \vec{J_z} =0,
\end{align}
```

whereby $F$ describes the non-uniform photochemical reaction, $\Phi(t)=1$ for all $t$ when illumination is switched on and 0 otherwise, $\Gamma$ is the rate of crosslinking 

```math
\begin{equation}
F(x,y,t) = k_p \left[ I_0 e^{-\zeta (T-y) } \right]^a \left\{ 1 +  e^{- \xi z} \cos \left[ \frac{2 \pi \cos \phi_r(t)}{\Lambda(t)} x - \frac{2 \pi \sin \phi_r(t)}{\Lambda(t)} y \right] \right\}.
\end{equation}
```

Here $k_p$ and $a$ are constants related to the rate of polymerization, $\zeta$ is the rate of absorption within the photosensitive layer, $xi$ is the scattering constant, $\Lambda$ and $T$ are the grating period and thickness respectively. The slant angle $\phi_r$ illustrated below from the assumptions of the fringe-plane rotation model.

![The geometry for a slanted holographic diffraction grating.](https://github.com/lyonsja-mathematics/PhD_Mathematical_Modelling/blob/main/fringe_plane_rotation_v3.png)

The diffusion flux vectors are

```math
\begin{align}
\vec{J_m} &=  -D_m \frac{\partial m}{\partial x}\vec{i} - D_m \frac{\partial m}{\partial y}\vec{j},\\
\vec{J_p} &=  -D_p \left\{ \left[ \frac{\partial p}{\partial x}\vec{i} + \frac{\partial p}{\partial y}\vec{j} \right] + \epsilon_{z} \left[ \frac{\partial (pz)}{\partial x}\vec{i} + \frac{\partial (pz)}{\partial y}\vec{j} \right] \right\},\\
\vec{J_z} &=  -D_z \left\{ \left[ \frac{\partial z}{\partial x}\vec{i} + \frac{\partial z}{\partial y}\vec{j} \right] + \epsilon_{z} \left[ \frac{\partial (pz)}{\partial x}\vec{i} + \frac{\partial (pz)}{\partial y}\vec{j} \right] + \epsilon_{z} \left[ \frac{\partial (qz)}{\partial x}\vec{i} + \frac{\partial (qz)}{\partial y}\vec{j} \right] \right\}.
\end{align}
```
The variable domain is
```math
\begin{align}
    0 &\leq x \leq \frac{\Lambda(t)}{ \cos \phi_r(t)}, & 0 &\leq y \leq T(t), & t &\geq 0.
\end{align}
```
# Boundary Immobilization

```math
\begin{align}
	\frac{\partial B}{\partial t} &= \frac{Y}{u}\frac{du}{dt}\frac{\partial B}{\partial Y},\\
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	\frac{\partial M}{\partial t} &= \frac{Y}{u}\frac{du}{dt}\frac{\partial M}{\partial Y} + \alpha^{(x)}_m \frac{\partial^2 m}{\partial x^2} +\alpha^{(y)}_m\frac{1}{u^2} \frac{\partial^2 M}{\partial Y^2} - \Phi(t) \beta F^*(x,Y, t) M ,\\
	%%%%%%%%%%%%%%%%%%%%%
	\frac{\partial P}{\partial t} &= \frac{Y}{u}\frac{du}{dt}\frac{\partial P}{\partial Y} + \alpha^{(x)}_p \frac{\partial^2 P}{\partial x^2} + \alpha^{(y)}_p\frac{1}{u^2} \frac{\partial^2 P}{\partial Y^2}+ \alpha^{(x)}_{pz}\frac{\partial^2 (PZ)}{\partial x^2} + \nonumber\\
	&\qquad \alpha^{(y)}_{pz}\frac{1}{u^2} \frac{\partial^2 (PZ)}{\partial Y^2} + \Phi \beta F^*(x, Y, t) M - \Phi(t) \gamma P^2,\\
	%%%%%%%%%%%%%%%%%%%%%%%
	\frac{\partial Q}{\partial t} &= \frac{Y}{u}\frac{du}{dt}\frac{\partial Q}{\partial Y} + \Phi(t) \gamma P^2,\\
	%%%%%%%%%%%%%%%%%%%%%%%%
	\frac{\partial Z}{\partial t} &= \frac{Y}{u}\frac{du}{dt}\frac{\partial Z}{\partial Y} +\alpha^{(x)}_z \frac{\partial^2 Z}{\partial x^2} + \alpha^{(y)}_z \frac{1}{u^2}\frac{\partial^2 Z}{\partial Y^2} + \alpha^{(x)}_{pz}\frac{\partial^2 (PZ)}{\partial x^2} + \nonumber\\
	&\qquad \frac{1}{u^2}\alpha^{(y)}_{pz} \frac{\partial^2 (PZ)}{\partial Y^2} + \alpha^{(x)}_{qz}\frac{\partial^2 (QZ)}{\partial x^2} + \alpha^{(y)}_{qz}\frac{1}{u^2} \frac{\partial^2 (QZ)}{\partial Y^2},
\end{align}
```

```math
\begin{equation}
    F^*(x, Y, t)= e^{-a\zeta^* u(1- Y)} \left\{1 + e^{- \xi^* Z} \cos \left[ 2 \pi\left(x - \frac{T_0}{\hat{x}}\tan \phi_r u Y\right)\right]\right\}
\end{equation}
```

# Initial & Boundary Conditions

```math
\begin{align}
    M(x,Y,0)&=1, & P(x,Y,0) &= 0, & Q(x,Y,0) &= 0, & Z(x,Y,0)&=1,\nonumber \\
    B(x,Y,0) &= 1, & u(0) &= 1, & u'(0) &= 0.
\end{align}
```

The variables $m$, $p$, $q$ and $z$ are all periodic on the $x$ domain, i.e. $m(x+\hat{x},y,t)=m(x,y,t)$, etc. A periodic boundary can be modelled by the condition

```math
\begin{align}
        \frac{\partial^n M}{\partial x^n}(0,Y,t) &=\frac{\partial^n M}{\partial x^n}(1,Y,t) & n &=\{0,1,2,...\},\\
        \frac{\partial^n P}{\partial x^n}(0,Y,t) &=\frac{\partial^n P}{\partial x^n}(1,Y,t) & n &=\{0,1,2,...\},\\
        \frac{\partial^n Z}{\partial x^n}(0,Y,t) &=\frac{\partial^n Z}{\partial x^n}(1,Y,t) & n &=\{0,1,2,...\}.
\end{align}
```

```math
\begin{align}
\frac{\partial M}{\partial Y}(x,0,t)&=\frac{\partial P}{\partial Y}(x,0,t)=\frac{\partial Q}{\partial Y}(x,0,t)=\frac{\partial Z}{\partial Y}(x,0,t)=\frac{\partial B}{\partial Y}(x,0,t)=0,\\
\frac{\partial M}{\partial Y}(x,1,t)&=\frac{\partial P}{\partial Y}(x,1,t)=\frac{\partial Q}{\partial Y}(x,1,t)=\frac{\partial Z}{\partial Y}(x,1,t)=\frac{\partial B}{\partial Y}(x,1,t)=0.
\end{align}
```

# Numerical Scheme

Numerical simulation can be done using the Crank-Nicolson implicit finite-difference scheme.

```math
\begin{align}
	\frac{M_{i,j}^{k+1} - M_{i,j}^k}{\Delta t} = &\frac{Y_j}{u_k}\frac{u_{k}-u_{k-1}}{\Delta t}\left(\frac{M_{i,j+1}^{k+1} - M_{i,j-1}^{k+1}}{4\Delta Y}+\frac{M_{i,j+1}^{k} - M_{i,j-1}^{k}}{4\Delta Y}\right) +\\
&\frac{\alpha_{mm}}{2}\left[ \frac{M_{i-1,j}^{k+1} - 2M_{i,j}^{k+1} + M_{i+1,j}^{k+1}}{\Delta x^2} + \frac{M_{i-1,j}^{k} - 2M_{i,j}^{k} + M_{i+1,j}^{k}}{\Delta x^2} \right] +\nonumber\\ 
&\frac{\alpha_{mm}}{2u_k^2}\left[ \frac{M_{i,j-1}^{k+1} - 2M_{i,j}^{k+1} + M_{i,j+1}^{k+1}}{\Delta Y^2} + \frac{M_{i,j-1}^{k} - 2M_{i,j}^{k} + M_{i,j+1}^{k}}{\Delta Y^2} \right] +\nonumber\\ 
&-\Phi^k \beta F_{i,j}^k \left(\frac{M_{i,j}^{k+1} + M_{i,j}^k}{2} \right),\nonumber
\end{align}
```
```math
\begin{align}
	\frac{P_{i,j}^{k+1} - P_{i,j}^k}{\Delta t} = &\frac{Y_j}{u_k}\frac{u_{k}-u_{k-1}}{\Delta t}\left(\frac{P_{i,j+1}^{k+1} - M_{i,j-1}^{k+1}}{4\Delta Y}+\frac{P_{i,j+1}^{k} - M_{i,j-1}^{k}}{4\Delta Y}\right) +\\
&\frac{\alpha_{pp}}{2}\left[ \frac{P_{i-1,j}^{k+1} - 2P_{i,j}^{k+1} + P_{i+1,j}^{k+1}}{\Delta x^2} + \frac{P_{i-1,j}^{k} - 2P_{i,j}^{k} + P_{i+1,j}^{k}}{\Delta x^2} \right] +\nonumber\\ 
&\frac{\alpha_{pp}}{2u_k^2}\left[ \frac{P_{i,j-1}^{k+1} - 2P_{i,j}^{k+1} + P_{i,j+1}^{k+1}}{\Delta Y^2} + \frac{P_{i,j-1}^{k} - 2P_{i,j}^{k} + P_{i,j+1}^{k}}{\Delta Y^2} \right] +\nonumber\\ 
&\frac{\alpha_{pz}}{2}\left[ \frac{Z^k_{i-1,j}P_{i-1,j}^{k+1} - 2Z^k_{i,j}P_{i,j}^{k+1} + Z^k_{i+1,j}P_{i+1,j}^{k+1}}{\Delta x^2} + \right.\\
&\left.\frac{Z^k_{i-1,j}P_{i-1,j}^{k} - 2Z^k_{i,j}P_{i,j}^{k} + Z^k_{i+1,j}P_{i+1,j}^{k}}{\Delta x^2} \right] + ... (\text{cut due to markdown constraints})
\end{align}
```

```math
\begin{align}
        \frac{Q_{i,j}^{k+1} - Q_{i,j}^k}{\Delta t} =  &\frac{Y_j}{u_k}\frac{u_{k}-u_{k-1}}{\Delta t}\left(\frac{Q_{i,j+1}^{k+1} - Q_{i,j-1}^{k+1}}{4\Delta Y}+\frac{Q_{i,j+1}^{k} - Q_{i,j-1}^{k}}{4\Delta Y}\right) + \Phi^k \gamma P^k_{i,j}\left(\frac{P^k_{i,j} + P^{k+1}_{i,j}}{2}\right),
\end{align}
```

```math
\begin{align}
	\frac{Z_{i,j}^{k+1} - Z_{i,j}^k}{\Delta t} = &\frac{Y_j}{u_k}\frac{u_{k}-u_{k-1}}{\Delta t}\left(\frac{Z_{i,j+1}^{k+1} - Z_{i,j-1}^{k+1}}{4\Delta Y}+\frac{Z_{i,j+1}^{k} - Z_{i,j-1}^{k}}{4\Delta Y}\right) +\\
&\frac{\alpha_{zz}}{2}\left[ \frac{Z_{i-1,j}^{k+1} - 2Z_{i,j}^{k+1} + Z_{i+1,j}^{k+1}}{\Delta x^2} + \frac{Z_{i-1,j}^{k} - 2Z_{i,j}^{k} + Z_{i+1,j}^{k}}{\Delta x^2} \right] +\nonumber\\ 
&\frac{\alpha_{zz}}{2u_k^2}\left[ \frac{Z_{i,j-1}^{k+1} - 2Z_{i,j}^{k+1} + Z_{i,j+1}^{k+1}}{\Delta Y^2} + \frac{Z_{i,j-1}^{k} - 2Z_{i,j}^{k} + Z_{i,j+1}^{k}}{\Delta Y^2} \right] +\nonumber\\ 
&\frac{\alpha_{zp}}{2}\left[ \frac{P^k_{i-1,j}Z_{i-1,j}^{k+1} - 2P^k_{i,j}Z_{i,j}^{k+1} + P^k_{i+1,j}Z_{i+1,j}^{k+1}}{\Delta x^2}\right] + \nonumber\\
&\frac{\alpha_{zp}}{2}\left[\frac{P^k_{i-1,j}Z_{i-1,j}^{k} - 2P^k_{i,j}Z_{i,j}^{k} + P^k_{i+1,j}Z_{i+1,j}^{k}}{\Delta x^2} \right] + ... (\text{cut due to markdown constraints})
\end{align}
```

```math
\begin{align}
    F_{i,j}^k &= \exp\left[-a \zeta^* u_k (1-Y_j)\right] \left[1 + \exp\left(-\xi^* Z_{i,j}^k\right) \cos\left(2 \pi x - 2 \pi \tan \phi_r^{k} u_k Y_j\right)\right].
\end{align}
```

```math
\begin{align}
    \int_0^1\int_0^1 M\text{ }dx\text{ }dY \approx M^*_k = &\frac{\Delta x^2}{4}\left[M_{0,0}^k + M_{0,J}^k + M_{J,0}^k + M_{J,J}^k + 2M_{0,1}^k + ... +2M^k_{0,J-1} + \right.\\
    &\left.2M^k_{1,0} + ... + 2M^k_{J-1,0} + 2M^k_{1,J} + ... + 2M^k_{J-1,J} + 2M^k_{J,1} + ... \right.\nonumber\\
    &\left.+ 2M^k_{J,J-1} + 4M^k_{2,2} + ... + 4M^k_{J-2,J-2} \right],\nonumber
\end{align}
```

```math
\begin{align}
    \int_0^1\int_0^1 P\text{ }dx\text{ }dY \approx P^*_k = &\frac{\Delta x^2}{4}\left[P_{0,0}^k + P_{0,J}^k + P_{J,0}^k + P_{J,J}^k + 2P_{0,1}^k + ... +2P^k_{0,J-1} + \right.\\
    &\left.2P^k_{1,0} + ... + 2P^k_{J-1,0} + 2P^k_{1,J} + ... + 2P^k_{J-1,J} + 2P^k_{J,1} + ... \right.\nonumber\\
    &\left.+ 2P^k_{J,J-1} + 4P^k_{2,2} + ... + 4P^k_{J-2,J-2} \right],\nonumber
\end{align}
```

```math
\begin{align}
    \int_0^1\int_0^1 Q\text{ }dx\text{ }dY \approx Q^*_k = &\frac{\Delta x^2}{4}\left[Q_{0,0}^k + Q_{0,J}^k + Q_{J,0}^k + Q_{J,J}^k + 2Q_{0,1}^k + ... +2Q^k_{0,J-1} + \right.\\
    &\left.2Q^k_{1,0} + ... + 2Q^k_{J-1,0} + 2Q^k_{1,J} + ... + 2Q^k_{J-1,J} + 2Q^k_{J,1} + ... \right.\nonumber\\
    &\left.+ 2Q^k_{J,J-1} + 4Q^k_{2,2} + ... + 4Q^k_{J-2,J-2} \right],\nonumber
\end{align}
```
```math
\begin{align}
    u_k &= \left[\frac{b_0}{\rho_b}+ \frac{1}{\rho_m} + \frac{z_0}{\rho_z}\right]^{-1}\left[\frac{b_0}{\rho_b} + \frac{M_k^*}{\rho_m} + \frac{P_k^*}{\rho_p} + \frac{Q^*_k}{\rho_p} + \frac{z_0}{\rho_z} \right],\\
     \phi_r^k &= \tan^{-1}\left(\frac{\tan \phi_r^0}{u_k}\right),\\
     \Lambda_k &= \Lambda_0 \frac{\cos \phi_r^k}{\cos \phi_r^0}.
\end{align}
```

# Refractive Index Modulation

```math
\begin{equation}
	\frac{n^2 - 1}{n^2 + 2} = \phi_m \frac{n_m^2 - 1}{n_m^2 + 2} + \phi_p \frac{n_p^2 - 1}{n_p^2 + 2} + \phi_q \frac{n_q^2 - 1}{n_q^2 + 2} + \phi_z \frac{n_z^2 - 1}{n_z^2 + 2} + \phi_b \frac{n_b^2 - 1}{n_b^2 + 2}.
\end{equation}
```

Solving the Lorentz-Lorenz equation will give the RI of the nanocomposite as a function of $x$ and $t$. The nanocomposite RI, $n(x,t)$, can be represented by a Fourier expansion series

```math
\begin{equation}
    n(x,y,t) \approx \sum_{i=0} A_i(y,t) \cos \left(\frac{2 \pi}{\Lambda} i x\right) + B_i(y,t) \sin \left(\frac{2 \pi}{\Lambda} i x\right),
\end{equation}
```

```math
\begin{align}
    A_0(y,t) &= \frac{1}{\Lambda}\int_0^{\Lambda} n(x,y,t)\text{ }dx,\\
    A_1(y,t) &= \frac{2}{\Lambda}\int_0^{\Lambda} n(x,y,t) \cos \left(\frac{2 \pi}{\Lambda} x\right)\text{ }dx,\\
    B_1(y,t) &= \frac{2}{\Lambda}\int_0^{\Lambda} n(x,y,t) \sin \left(\frac{2 \pi}{\Lambda} x\right)\text{ }dx.
\end{align}
```

RI modulation can be modelled as

```math
\begin{equation}
    \Delta n(y, t) = 2\sqrt{A_1^2 + B_1^2}.
\end{equation}
```

# Shrinkage Modelling

We can calculate the volume at time $t$ if we have expressions for the total volume of monomer, short polymer, cross-linked polymer and nanoparticles inside the grating

```math
\begin{align}
    v(t)=&\frac{1}{\rho_m}\left[\frac{1}{\hat{x} T(t)}\int_0^{T(t)}\int_0^{\hat{x}}m\text{ }dx\text{ }dy\right] + \frac{1}{\rho_p}\left[\frac{1}{\hat{x} T(t)}\int_0^{T(t)}\int_0^{\hat{x}}p\text{ }dx\text{ }dy\right] +\nonumber\\
    &\frac{1}{\rho_p}\left[\frac{1}{\hat{x} T(t)}\int_0^{T(t)}\int_0^{\hat{x}}q\text{ }dx\text{ }dy\right] + \frac{1}{\rho_z}\left[\frac{1}{\hat{x} T(t)}\int_0^{T(t)}\int_0^{\hat{x}}z\text{ }dx\text{ }dy\right] +\nonumber\\
    &\frac{1}{\rho_b}\left[\frac{1}{\hat{x} T(t)}\int_0^{T(t)}\int_0^{\hat{x}}b\text{ }dx\text{ }dy\right].
\end{align}
```

An important assumptions of the fringe-plane rotation model is that all loss of volume due to polymerization takes place in the thickness of the recording medium

```math
\begin{equation}
    u(t)=\frac{T(t)}{T_0}= \frac{v(t)}{v(0)}.
\end{equation}
```
```math
\begin{equation}
    u(t)=\left[\frac{1}{\rho_b}\frac{b_0}{m_0} + \frac{1}{\rho_m} + \frac{1}{\rho_z}\frac{z_0}{m_0}\right]^{-1}\left[\int_0^{1}\int_0^{1}\frac{M}{\rho_m} + \frac{P}{\rho_p} + \frac{Q}{\rho_p} + \frac{z_0/m_0 Z}{\rho_z} + \frac{b_0/m_0}{\rho_b}\text{ }dx\text{ }dY\right],
\end{equation}
```

```math
\begin{equation}
    \text{Actual Shrinkage} = \frac{u(0) - u(t)}{u(0)} = 1 - u(t),
\end{equation}
```

```math
\begin{align}
        \phi_r(t) &= \tan^{-1}\left[\frac{\tan \phi_r(0)}{u(t)}\right],\\
        \Lambda(t) &= \hat{x} \cos \phi_r(t),
\end{align}
```

```math
\begin{equation}
    \overline{n}(t)= \frac{1}{\hat{x} T} \int_0^T \int_0^{\hat{x}} n(x,y,t)\text{ }dx\text{ }dy = \int_0^{1} \int_0^{1} n(x,Y,t)\text{ }dx\text{ }dY,
\end{equation}
```

```math
\begin{equation}
    \theta_B(t) = \sin^{-1}\left(\frac{\lambda_r}{2 \overline{n}(t) \Lambda(t)}\right) - \phi_r(t),
\end{equation}
```

```math
\begin{equation}
    \text{Apparent Shrinkage} = 1 - \frac{\tan  \phi_r(0)}{\tan \left[\phi_r(0) + \Delta \theta_B\right]}.
\end{equation}
```

# Results

Numerical simulation of these equations can be executed via the Python script [holographic_grating.py](https://github.com/lyonsja-mathematics/PhD_Mathematical_Modelling/blob/main/holographic_grating.py)

![Observed and predicted time evolution of RI modulation of an unslanted holographic grating recorded in AA/PVA photopolymer with increased doping of BEA nanozeolites.](https://github.com/lyonsja-mathematics/PhD_Mathematical_Modelling/blob/main/actual_and_predicted_cody.png)

![The predicted actual and apparent shrinkage for slant angles up to 10 degrees and spatial frequencies ranging from 250 - 2000 lines/mm](https://github.com/lyonsja-mathematics/PhD_Mathematical_Modelling/blob/main/shrinkage_angle_lpmm.png)
