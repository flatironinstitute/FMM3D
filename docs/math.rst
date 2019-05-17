Definitions 
===========
Let $\mathbf{I}_{A}$ denote the indicator function of a set $A$,
and let $\varepsilon_{\textrm{mach}}$ denote machine precision ($2^{-52}$).
Let $y_{i} \in \mathbb{R}^{3}$, $i=1,2,\ldots N$, denote a collection
of source locations and let $x_{i} \in \mathbb{R}^{3}$ denote a collection
of target locations. Suppose that the source and target locations
are contained in a cube of side length $L$.


Laplace FMM
***********
Let $c_{i} \in \mathbb{R}$ 
$i=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{i} \in \mathbb{R}^{3}$,
$i=1,2,\ldots N$, 
denote a collection of dipole strengths.
Let $G_{0}(x): \mathbb{R}^{3} \to \mathbb{R}$ denote 
the scaled Green's function for Laplace's equation given by

.. math::

   G_{0}(x) = \frac{1}{|x|} \mathbf{I}_{|x|>\varepsilon}\, ,

where $\varepsilon = \varepsilon_{\text{mach}} L$.

The Laplace FMM computes 
the potential $u(x)$ and the it's gradient $\nabla u(x)$,
at the source and target locations where $u(x)$ is defined 
by the formula,

.. math::
    :label: lap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} G_{0}(x-y_{j}) - v_{j} \cdot \nabla G_{0}(x-y_{j}) \, .


Helmholtz FMM
*************
Let $c_{i} \in \mathbb{C}$ 
$i=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{i} \in \mathbb{C}^{3}$,
$i=1,2,\ldots N$, 
denote a collection of dipole strengths.
Let $k\in\mathbb{C}$ denote the wave number or the Helmholtz 
parameter. 
Let $G_{k}(x): \mathbb{R}^{3} \to \mathbb{C}$ denote 
the scaled Green's function for Helmholtz's equation given by

.. math::

    G_{k}(x) = \frac{e^{ik |x|}}{|x|} \mathbf{I}_{|x|>\varepsilon}\, ,

where $\varepsilon = \varepsilon_{\textrm{mach}} \lvert \omega \rvert  L$. 

The Helmholtz FMM computes 
the potential $u(x)$ and the it's gradient $\nabla u(x)$,
at the source and target locations where $u(x)$ is defined 
by the formula,

.. math::
   :label: helm_nbody

    u(x) = \sum_{j=1}^{N} c_{j} G_{k}(x-y_{j}) - v_{j} \cdot \nabla G_{k}(x-y_{j}) \, .

Vectorized versions   
*******************
The vectorized versions of the Laplace and Helmholtz FMM, compute the $n_{d}$ collection
of potentials corresponding to $n_{d}$ charge or dipole densities, located
at the same set of source and target locations. 
For example, for the Laplace FMM, let $c_{\ell,j}\in\mathbb{R}$, $j=1,2,\ldots N$, $\ell=1,2,\ldots n_{d}$
denote a collection of $n_{d}$ charge densities, and
let $v_{\ell,j} \in \mathbb{R}^{3}$ denote a collection of $n_{d}$ dipole densities, 
then the vectorized Laplace FMM computes the potentials $u_{\ell}(x)$ 
and it's gradients $\nabla u_{\ell}(x)$ defined by the formula

.. math::
    :label: lap_nbody_vec

    u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} G_{0}(x-y_{j}) - v_{\ell,j} \cdot \nabla G_{0}(x-y_{j}) \, , \, \quad \ell=1,2,\ldots n_{d}.

