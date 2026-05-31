This repository contains the final version of the programming scripts I wrote for my PhD thesis on *Mathematical Modelling of Hybrid Photonic Structures for Holographic Sensors*. My research was funded by the Irish Research Council (GOIPG/2021/214).

![The formation of holographic diffraction gratings in hybrid photopolymer media](https://github.com/lyonsja-mathematics/PhD_Mathematical_Modelling/blob/main/Holographic%20Photopolymerization%20v3.png)

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
