Definitions 
===========
Let $x_{j} \in \mathbb{R}^{3}$, $i=1,2,\ldots N$, denote a collection
of source locations and let $t_{i} \in \mathbb{R}^{3}$ denote a collection
of target locations. 


Laplace FMM
***********
Let $c_{j} \in \mathbb{R}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{j} \in \mathbb{R}^{3}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths.

The Laplace FMM computes 
the potential $u(x)$ and the its gradient $\nabla u(x)$
given by

.. math::
    :label: lap_nbody

    u(x) = \sum_{j=1}^{N} \frac{c_{j}}{\|x-x_{j}\|} - v_{j} \cdot \nabla \left( \frac{1}{\|x-x_{j}\|}\right)  \, , 

at the source and target locations. When $x=x_{j}$, the term
corresponding to $x_{j}$ is dropped from the sum.

Helmholtz FMM
*************
Let $c_{j} \in \mathbb{C}$, 
$j=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{j} \in \mathbb{C}^{3}$,
$j=1,2,\ldots N$, 
denote a collection of dipole strengths.
Let $k\in\mathbb{C}$ denote the wave number or the Helmholtz 
parameter. 

The Helmholtz FMM computes 
the potential $u(x)$ and the its gradient $\nabla u(x)$
given by

.. math::
   :label: helm_nbody

    u(x) = \sum_{j=1}^{N} c_{j} \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} - v_{j} \cdot \nabla \left( \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)  \, , 

at the source and target locations. When $x=x_{j}$, the term
corresponding to $x_{j}$ is dropped from the sum.

Vectorized versions   
*******************
The vectorized versions of the Laplace and Helmholtz FMM, 
computes repeated FMMs for new charge and dipole strengths
located at the same source locations, where the potential and its
gradient are evaluated at the same set of target locations.

For example, for the vectorized Laplace FMM, let $c_{\ell,j}\in\mathbb{R}$, 
$j=1,2,\ldots N$, $\ell=1,2,\ldots n_{d}$
denote a collection of $n_{d}$ charge strengths, and
let $v_{\ell,j} \in \mathbb{R}^{3}$ denote a collection of $n_{d}$ dipole strengths. 
Then the vectorized Laplace FMM computes the potentials $u_{\ell}(x)$ 
and its gradients $\nabla u_{\ell}(x)$ defined by the formula

.. math::
    :label: lap_nbody_vec

    u_{\ell}(x) = \sum_{j=1}^{N} \frac{c_{\ell,j}}{\|x-x_{j}\|} - v_{\ell,j} \cdot \nabla \left( \frac{1}{\|x-x_{j}\|}\right)  \, , \quad \ell=1,2,\ldots n_{d}\,

at the source and target locations. 

Similarly, for the vectorized Helmholtz FMM, let $c_{\ell,j}\in\mathbb{C}$, 
$j=1,2,\ldots N$, $\ell=1,2,\ldots n_{d}$
denote a collection of $n_{d}$ charge strengths, and
let $v_{\ell,j} \in \mathbb{C}^{3}$ denote a collection of $n_{d}$ dipole strengths. 
Then the vectorized Helmholtz FMM computes the potentials $u_{\ell}(x)$ 
and its gradients $\nabla u_{\ell}(x)$ defined by the formula

.. math::
    :label: helm_nbody_vec

    u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} - v_{\ell,j} \cdot \nabla \left( \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)  \, ,\quad \ell =1,2,\ldots n_{d}  

at the source and target locations. 

.. note::

   In double precision arithmetic, two numbers which are
   within machine precision of each other cannot be
   distinguished. In order to account for this, suppose that the sources
   and targets are contained in a cube with side length $L$, then
   for all $x$ such that $\| x-x_{j} \| \leq L \varepsilon_{\textrm{mach}}$,
   the term corresponding to $x_{j}$ is dropped from the sum.
   Here $\varepsilon_{\textrm{mach}} = 2^{-52}$ is machine precision.

