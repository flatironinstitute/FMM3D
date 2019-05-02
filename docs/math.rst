Mathematical definitions of transforms
======================================

Vectorized Laplace FMM
**************
Let $y_{i} \in \mathbb{R}^{3}$, $i=1,2,\ldots N$, 
denote a collection of ``source'' locations, $c_{i,j} \in \mathbb{R}$ 
$i=1,2,\ldots n_{d}$, $j=1,2,\ldots N$, 
denote a collection of $n_{d}$ charge densities, $v_{i,j} \in \mathbb{R}^{3}$,
$i=1,2,\ldots n_{d}$, $j=1,2,\ldots N$, 
denote a collection of $n_{d}$ dipole densities.
Note that the soruce locations are the same for the different 
charge and dipole densities. 
Let $x_{i} \in \mathbb{R}^{3}$, $i=1,2,\ldots M$, denote a collection 
of ``target'' locations. 
Let $G_{0}(x): \mathbb{R}^{3} \to \mathbb{R}$ denote the scaled Green's function for Laplace's equation
given by
.. math::
   G_{0}(x) = \frac{1}{|x|} \bI_{|x|>\varepsilon}\, ,

where $\bI$ is the indicator function of a set, and if the sources
and targets are contained in a cube of side length $L$, then
$\varepsilon = \varepsilon_{\text{mach}} L$, with 
$\varepsilon_{\text{mach}}$ is machine precision ($2^{-52}$). 

The Laplace $N-$body interaction is the computation of the 
the potential $u_{i}(x)$ and the it's gradient $\nabla u_{i}(x)$,
$i=1,2,\ldots n_{d}$
at the source and target locations where $u(x)$ is defined 
by the formula,

.. math::
   :label: lap_nbody
   u_{i}(x) = \sum_{j=1}^{N} c_{i,j} G(x-y_{j}) 
   - v_{i,j} \cdot \nabla G(x-y_{j}) \, .

