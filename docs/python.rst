.. _pyt:

Python
=======

The Python interface has four callable subroutines:

*  `Laplace wrappers <python.html#lap-pyt>`__: Fast multipole implementation (lfmm3d) and direct evaluation (l3ddir) for Laplace N-body interactions
*  `Helmholtz wrappers <python.html#helm-pyt>`__: Fast multipole implementation (hfmm3d) and direct evaluation (h3ddir) for Helmholtz N-body interactions
*  `Stokes wrappers <python.html#stok-pyt>`__: Fast multipole implementation (stfmm3d) and direct evaluation (st3ddir) for Stokes N-body interactions
*  `Maxwell wrappers <python.html#em-pyt>`__: Fast multipole implementation (emfmm3d) and direct evaluation (em3ddir) for Maxwell N-body interactions


.. _lap-pyt:

Laplace wrappers
*******************


This subroutine computes the N-body Laplace
interactions and its gradients in three dimensions where 
the interaction kernel is given by $\frac{1}/{4\pi r}$
 
.. math::

    u(x) = \sum_{j=1}^{N} \frac{c_{j}}{4\pi\|x-x_{j}\|} - v_{j} \cdot \nabla \left( \frac{1}{4\pi\|x-x_{j}\|}\right)   

where $c_{j}$ are the charge densities
$v_{j}$ are the dipole orientation vectors, and
$x_{j}$ are the source locations.
When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
from the sum.

.. code:: python
   
   def lfmm3d(*,eps,sources,charges=None,dipvec=None,targets=None,pg=0,pgt=0,nd=1)

Wrapper for fast multipole implementation for Laplace N-body
interactions.

Args:

-  eps: double   
      precision requested
-  sources: double(3,n)    
     source locations, $x_{j}$
-  charges: double(n,) or double(nd,n) 
     charge densities, $c_{j}$ 
-  dipvec: double(3,n) or double(nd,3,n)
     dipole orientation vectors, $v_{j}$ 
-  nd: integer
     number of charge/dipole vectors 
-  pg: integer
      | source eval flag
      | potential at sources evaluated if pg = 1
      | potenial and gradient at sources evaluated if pg=2
-  targets: double(3,nt)
      target locations ($t_{i}$) (optional)
-  pgt: integer
      | target eval flag (optional)
      | potential at targets evaluated if pgt = 1
      | potenial and gradient at targets evaluated if pgt=2  

Returns:
The subroutine returns an object out of type Output with the following
variables

-  out.pot: potential at source locations, if requested, $u(x_{j})$
-  out.grad: gradient at source locations, if requested, $\nabla u(x_{j})$
-  out.pottarg: potential at target locations, if requested, $u(t_{i})$
-  out.gradtarg: gradient at target locations, if requested, $\nabla u(t_{i})$

------------------------------------------------------------------

Wrapper for direct evaluation of Laplace N-body interactions.
Note that this wrapper only returns potentials and gradients at the
target locations.
              
.. code:: python
   
   def l3ddir(*,sources,charges=None,dipvec=None,targets=None,pgt=0,nd=1)

------------------------------------------------------------------

Example:

-  see ``lfmm3d_example.py``

.. container:: rttext

  `Back to top <python.html#pyt>`__


.. _helm-pyt:

Helmholtz wrappers
*******************


This subroutine computes the N-body Helmholtz
interactions and its gradients in three dimensions where 
the interaction kernel is given by $\frac{e^{ikr}}{4\pi r}$
 
.. math::

    u(x) = \sum_{j=1}^{N} c_{j} \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} - v_{j} \cdot \nabla \left( \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|}\right)   

where $c_{j}$ are the charge densities
$v_{j}$ are the dipole orientation vectors, and
$x_{j}$ are the source locations.
When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
from the sum.

.. code:: python
   
   def hfmm3d(*,eps,zk,sources,charges=None,dipvec=None,targets=None,pg=0,pgt=0,nd=1)

Wrapper for fast multipole implementation for Helmholtz N-body
interactions.

Args:

-  eps: double   
      precision requested
-  zk: complex
      Helmholtz parameter, k
-  sources: double(3,n)    
     source locations, $x_{j}$
-  charges: complex(n,) or complex(nd,n) 
     charge densities, $c_{j}$
-  dipvec: complex(3,n) or complex(nd,3,n)
     dipole orientation vectors, $v_{j}$ 
-  nd: integer
     number of charge/dipole vectors 
-  pg: integer
      | source eval flag
      | potential at sources evaluated if pg = 1
      | potenial and gradient at sources evaluated if pg=2
-  targets: double(3,nt)
      target locations, $t_{i}$ (optional)
-  pgt: integer
      | target eval flag (optional)
      | potential at targets evaluated if pgt = 1
      | potenial and gradient at targets evaluated if pgt=2  

Returns:
The subroutine returns an object out of type Output with the following
variables

-  out.pot: potential at source locations, if requested, $u(x_{j})$
-  out.grad: gradient at source locations, if requested, $\nabla u(x_{j})$
-  out.pottarg: potential at target locations, if requested, $u(t_{i})$
-  out.gradtarg: gradient at target locations, if requested, $\nabla u(t_{i})$

------------------------------------------------------------------

Wrapper for direct evaluation of Helmholtz N-body interactions.
Note that this wrapper only returns potentials and gradients at the
target locations.
              
.. code:: python
   
   def h3ddir(*,zk,sources,charges=None,dipvec=None,targets=None,pgt=0,nd=1)

------------------------------------------------------------------

Example:

-  see ``hfmm3d_example.py``

.. container:: rttext

  `Back to top <python.html#pyt>`__


.. _stok-pyt:

Stokes wrappers
*******************


Let $\mathcal{G}^{\textrm{stok}}(x,y)$ 
denote the Stokeslet given by


.. math::
   \mathcal{G}^{\textrm{stok}}(x,y)=\frac{1}{8\pi \|x-y\|^3}
   \begin{bmatrix}
   (x_{1}-y_{1})^2 + \|x-y \|^2 & (x_{1}-y_{1})(x_{2}-y_{2}) &
   (x_{1}-y_{1})(x_{3}-y_{3}) \\ 
   (x_{2}-y_{2})(x_{1}-y_{1}) & (x_{2}-y_{2})^2 + \|x-y \|^2 & 
   (x_{2}-y_{2})(x_{3}-y_{3}) \\ 
   (x_{3}-y_{3})(x_{1}-y_{1})  & (x_{3}-y_{3})(x_{2}-y_{2}) & 
   (x_{3}-y_{3})^2 + \|x-y \|^2 
   \end{bmatrix} \, ,

let $\mathcal{T}^{\textrm{stok}}(x,y)$ denote the Stresslet whose action on
a vector $v$ is given by

.. math::
   v\cdot \mathcal{T}^{\textrm{stok}}(x,y)  = 
   \frac{3 v \cdot (x-y)}{4\pi\|x-y \|^5}
   \begin{bmatrix}
   (x_{1}-y_{1})^2 & (x_{1}-y_{1})(x_{2}-y_{2}) &
   (x_{1}-y_{1})(x_{3}-y_{3}) \\ 
   (x_{2}-y_{2})(x_{1}-y_{1}) & (x_{2}-y_{2})^2 & 
   (x_{2}-y_{2})(x_{3}-y_{3}) \\ 
   (x_{3}-y_{3})(x_{1}-y_{1})  & (x_{3}-y_{3})(x_{2}-y_{2}) & 
   (x_{3}-y_{3})^2  
   \end{bmatrix} \, ,

let $\mathcal{R}^{\textrm{stok}}(x,y)$ denote the Rotlet whose action on
a vector $v$ is given by

.. math::
   v\cdot \mathcal{R}^{\textrm{stok}}(x,y)  = 
   \frac{v \cdot (x-y)}{4\pi\|x-y \|^3}
   \begin{bmatrix}
   1 & 0 & 0 \\ 
   0 & 1 & 0 \\ 
   0 & 0 & 1  
   \end{bmatrix} \, ,

and $\mathcal{D}^{\textrm{stok}}(x,y)$ denote the symmetric part of Doublet whose action on
a vector $v$ is given by

.. math::
   v\cdot \mathcal{D}^{\textrm{stok}}(x,y)  = 
   \frac{3 v \cdot (x-y)}{4\pi\|x-y \|^5}
   \begin{bmatrix}
   (x_{1}-y_{1})^2 & (x_{1}-y_{1})(x_{2}-y_{2}) &
   (x_{1}-y_{1})(x_{3}-y_{3}) \\ 
   (x_{2}-y_{2})(x_{1}-y_{1}) & (x_{2}-y_{2})^2 & 
   (x_{2}-y_{2})(x_{3}-y_{3}) \\ 
   (x_{3}-y_{3})(x_{1}-y_{1})  & (x_{3}-y_{3})(x_{2}-y_{2}) & 
   (x_{3}-y_{3})^2  
   \end{bmatrix} - \\
   \frac{1}{4\pi\|x-y \|^3}
   \begin{bmatrix}
   v_1(x_{1}-y_{1}) & v_2(x_{1}-y_{1}) & v_3(x_{1}-y_{1}) \\
   v_2(x_{2}-y_{2}) & v_2(x_{2}-y_{2}) & v_3(x_{2}-y_{2}) \\
   v_3(x_{3}-y_{3}) & v_3(x_{3}-y_{3}) & v_3(x_{3}-y_{3})
   \end{bmatrix} \, . 

This subroutine computes the N-body Stokes
interactions, its gradients and the corresponding pressure 
in three dimensions given by 


.. math::

    u(x) = \sum_{m=1}^{N} \mathcal{G}^{\textrm{stok}}(x,x_{j}) \sigma_{j}  + \nu_{j} \cdot \mathcal{T}^{\textrm{stok}}(x,x_{j}) \cdot \mu_{j} +
                          \nu^{r}_{j} \cdot \mathcal{R}^{\textrm{stok}}(x,x_{j}) \cdot \mu^{r}_{j} - \mu^{r}_{j} \cdot \mathcal{R}^{\textrm{stok}}(x,x_{j}) \cdot \nu^{r}_{j} + \\
                          \nu^{d}_{j} \cdot \mathcal{R}^{\textrm{stok}}(x,x_{j}) \cdot \mu^{d}_{j} - \mu^{d}_{j} \cdot \mathcal{R}^{\textrm{stok}}(x,x_{j}) \cdot \nu^{d}_{j} +
                          \nu^{d}_{j} \cdot \mathcal{D}^{\textrm{stok}}(x,x_{j}) \cdot \mu^{d}_{j} \, .

where $\sigma_{j}$ are the Stokeslet densities,
$\nu_{j}$ are the stresslet orientation vectors, $\mu_{j}$ 
are the stresslet densities,
$\nu^{r}_{j}$ are the rotlet orientation vectors, $\mu^{r}_{j}$ 
are the rotlet densities,
$\nu^{d}_{j}$ are the doublet orientation vectors, $\mu^{d}_{j}$ 
are the doublet densities, and
$x_{j}$ are the source locations.
When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
from the sum.

.. code:: python
   
   def stfmm3d(*,eps,sources,stoklet=None,strslet=None,strsvec=None,rotlet=None,rotvec=None,doublet=None,doubvec=None,targets=None,ifppreg=0,ifppregtarg=0,nd=1)

Wrapper for fast multipole implementation for Stokes N-body
interactions.

Args:

-  eps: double   
      precision requested
-  sources: float(3,n)   
      source locations
-  stoklet: float(nd,3,n) or float(3,n)
      Stokeslet charge strengths ($\sigma_{j}$ above)
-  strslet: float(nd,3,n) or float(3,n)
      stresslet strengths ($\mu_{j}$ above)
-  strsvec: float(nd,3,n) or float(3,n)
      stresslet orientations ($\nu_{j}$ above)
-  rotlet: float(nd,3,n) or float(3,n)
      rotlet strengths ($\mu^{r}_{j}$ above)
-  rotvec: float(nd,3,n) or float(3,n)
      rotlet orientations ($\nu^{r}_{j}$ above)
-  doublet: float(nd,3,n) or float(3,n)
      doublet strengths ($\mu^{d}_{j}$ above)
-  doubvec: float(nd,3,n) or float(3,n)
      doublet orientations ($\nu^{d}_{j}$ above)
-  targets: float(3,nt)
      target locations (x)
-  ifppreg: integer
      | flag for evaluating potential, gradient, and pressure at sources
      | potential at sources evaluated if ifppreg = 1
      | potential and pressure at sources evaluated if ifppreg=2
      | potential, pressure and gradient at sources evaluated if ifppreg=3
-  ifppregtarg: integer
      | flag for evaluating potential, gradient, and pressure at targets
      | potential at targets evaluated if ifppregtarg = 1
      | potential and pressure at targets evaluated if ifppregtarg = 2 
      | potential, pressure and gradient at targets evaluated if ifppregtarg = 3

Returns:

-  out.pot: velocity at source locations if requested
-  out.pre: pressure at source locations if requested
-  out.grad: gradient of velocity at source locations if requested
-  out.pottarg: velocity at target locations if requested
-  out.pretarg: pressure at target locations if requested
-  out.gradtarg: gradient of velocity at target locations if requested

------------------------------------------------------------------

Wrapper for direct evaluation of Stokes N-body interactions. 
Note that this wrapper only returns potentials and gradients at the
target locations.
              
.. code:: python
   
   def st3ddir(*,eps,sources,stoklet=None,strslet=None,strsvec=None,rotlet=None,rotvec=None,doublet=None,doubvec=None,targets=None,ifppreg=0,ifppregtarg=0,nd=1,thresh=1e-16):

------------------------------------------------------------------

Example:

-  see ``stfmm3d_example.py``

.. container:: rttext

  `Back to top <python.html#pyt>`__



.. _em-pyt:

Maxwell wrappers
*******************


This subroutine computes the N-body Maxwell
interactions, its curl and its divergence in three dimensions
given by
 
.. math::

    E(x) = \sum_{j=1}^{N} \nabla \times \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} M_{j} + \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} J_{j} +  \nabla \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} \rho_{j}       

where $M_{j}$ are the magnetic current densities,
$J_{j}$ are the electric current densities, 
$\rho_{j}$ are the electric charge densities, and
$x_{j}$ are the source locations.
When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
from the sum.

.. code:: python
   
   def emfmm3d(*,eps,zk,sources,h_current=None,e_current=None,e_charge=None,targets=None,ifE=0,ifcurlE=0,ifdivE=0,nd=1):

Wrapper for fast multipole implementation for Maxwell N-body
interactions.
Note that this wrapper only returns fields, divergences, and curls at the
target locations.

Args:

-  eps: double   
      precision requested
-  zk: complex
      Wavenumber, k
-  sources: float(3,n)   
      source locations
-  h_current: complex(3,n) or complex(nd,3,n)
      Magnetic currents, $M_{j}$
-  e_current: complex(3,n) or complex(nd,3,n)
      Electric currents, $J_{j}$
-  e_charge: complex(n,) or complex(nd,n)
      Electric charges, $\rho_{j}$
-  targets: float(3,nt)
      target locations, $t_{i}$ 
-  ifE: integer
      E is returned at the target locations if ifE = 1
-  ifcurlE: integer
      curl E is returned at the target locations if ifcurlE = 1
-  ifdivE: integer
      div E is returned at the target locations if ifdivE = 1

Returns:

-  out.E: E field defined above at target locations if requested $(E(t_{j}))$
-  out.curlE: curl of E field at target locations if requested $(\nabla \times E(t_{j}))$
-  out.divE: divergence of E at target locations if requested $(\nabla \cdot E(t_{j}))$

------------------------------------------------------------------

Wrapper for direct evaluation of Maxwell N-body interactions.
Note that this wrapper only returns fields, divergences, and curls at the
target locations.
              
.. code:: python
   
   def em3ddir(*,eps,zk,sources,h_current=None,e_current=None,e_charge=None,targets=None,ifE=0,ifcurlE=0,ifdivE=0,nd=1,thresh=1e-16):

------------------------------------------------------------------

Example:

-  see ``emfmm3d_example.py``

.. container:: rttext

  `Back to top <python.html#pyt>`__

