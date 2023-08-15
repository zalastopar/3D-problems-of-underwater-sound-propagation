# 3D problems of underwater sound propagation

Project for seminar Mathematische Modellierung in Anwendungen at Bergische Universität Wuppertal. 


Human presence on the sea is mainly on the surface, but we also affect the life
under it. Noise caused by propellers on the ships loses its intensity, when it
crosses the air-sea border, and the underwater sound we cause is rarely taken
into consideration. 

## Equations
Acoustic pressure $P(x, y, z)$ is described by a boundary value problem with the Helmholtz equation, which characterizes the sound propagation in shallow waters: 
```math
\begin{equation}
    P_{xx}+ P_{yy} + P_{zz} + \frac{\omega^2}{c^2} P = - \delta(x, y, z - z_s)\,,
\end{equation}
```
$\omega = 2\pi f$ angular frequency. $P|_{z=0} = 0$; at the sea bottom, described by function $h(x, y)$, $P_z|_{z=h} = 0 \,. Sveshnikov partial
radiation boundary condition at radius to infinity.


## Summary of the algorithm
\text{i}
* Calculate modes: $K_j = K_j(x, y)\,,$ $\phi_j = \phi_j(z, x, y)\,.$
* Apply Crank-Nicholson algorithm to PDE: $2\text{i} K_j^0 A_x + A_{yy} + (K_j^2 - (K_j^0)^2)A = 0\,.$
* Solve tridiagonal system: $Qa_{n+1} = Pa_{n}\,.$
* Repeat for all $j$.
* $$P(x, y, z) = \sum_{j=1}^{M} \text{e}^{\text{i} xK_j^0} A_j(x, y) \phi_j(z)\,.$$
* Apply Perfectly matched layer (PML).
