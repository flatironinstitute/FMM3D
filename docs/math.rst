Mathematical definitions of transforms
======================================

Laplace FMM
***********
Let $y_{i} \in \mathbb{R}^{3}$, $i=1,2,\ldots N$, 
denote a collection of ``source'' locations, $c_{i} \in \mathbb{R}$ 
$i=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{i} \in \mathbb{R}^{3}$,
$i=1,2,\ldots N$, 
denote a collection of dipole strengths.
Let $x_{i} \in \mathbb{R}^{3}$, $i=1,2,\ldots M$, denote a collection 
of ``target'' locations. 
Let $G_{0}(x): \mathbb{R}^{3} \to \mathbb{R}$ denote 
the scaled Green's function for Laplace's equation given by

.. math::

   G_{0}(x) = \frac{1}{|x|} \mathbf{I}_{|x|>\varepsilon}\, ,

where $\mathbf{I}$ is the indicator function of a set, and if the sources
and targets are contained in a cube of side length $L$, then
$\varepsilon = \varepsilon_{\text{mach}} L$, with 
$\varepsilon_{\text{mach}}$ is machine precision ($2^{-52}$). 

The Laplace $N-$body interaction computes 
the potential $u(x)$ and the it's gradient $\nabla u(x)$,
at the source and target locations where $u(x)$ is defined 
by the formula,

.. math::
    :label: lap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} G_{0}(x-y_{j}) - v_{j} \cdot \nabla G_{0}(x-y_{j}) \, .


Helmholtz FMM
*************
Let $y_{i} \in \mathbb{R}^{3}$, $i=1,2,\ldots N$, 
denote a collection of ``source'' locations, $c_{i} \in \mathbb{C}$ 
$i=1,2,\ldots N$, 
denote a collection of charge strengths, $v_{i} \in \mathbb{C}^{3}$,
$i=1,2,\ldots N$, 
denote a collection of dipole strengths.
Let $x_{i} \in \mathbb{R}^{3}$, $i=1,2,\ldots M$, denote a collection 
of ``target'' locations.
Let $\omega\in\mathbb{C}$ denote the wave number or the Helmholtz 
parameter. 
Let $G_{\omega}(x): \mathbb{R}^{3} \to \mathbb{C}$ denote 
the scaled Green's function for Helmholtz's equation given by

.. math::

    G_{\omega}(x) = \frac{e^{i\omega |x|}}{|x|} \mathbf{I}_{|x|>\varepsilon}\, ,

where $\mathbf{I}$ is the indicator function of a set, and if the sources
and targets are contained in a cube of side length $L$, then $\varepsilon = \varepsilon_{\textrm{mach}} \lvert \omega \rvert  L$, 
with $\varepsilon_{\textrm{mach}}$ is machine precision ($2^{-52}$). 

The Helmholtz $N-$body interaction computes 
the potential $u(x)$ and the it's gradient $\nabla u(x)$,
at the source and target locations where $u(x)$ is defined 
by the formula,

.. math::
   :label: lap_nbody

    u(x) = \sum_{j=1}^{N} c_{j} G_{\omega}(x-y_{j}) - v_{j} \cdot \nabla G_{\omega}(x-y_{j}) \, .

Vectorized versions   
*******************
