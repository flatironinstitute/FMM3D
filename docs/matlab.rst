.. _mat:

MATLAB
=======

The MATLAB interface has four callable subroutines:

*  `Laplace wrappers <matlab.html#lap-mat>`__: Fast multipole implementation (lfmm3d) and direct evaluation (l3ddir) for Laplace N-body interactions
*  `Helmholtz wrappers <matlab.html#helm-mat>`__: Fast multipole implementation (hfmm3d) and direct evaluation (h3ddir) for Helmholtz N-body interactions
*  `Stokes wrappers <matlab.html#stok-mat>`__: Fast multipole implementation (stfmm3d) and direct evaluation (st3ddir) for Stokes N-body interactions
*  `Maxwell wrappers <matlab.html#em-mat>`__: Fast multipole implementation (emfmm3d) and direct evaluation (em3ddir) for Maxwell N-body interactions


.. _lap-mat:

Laplace wrappers
*******************


This subroutine computes the N-body Laplace
interactions and its gradients in three dimensions where 
the interaction kernel is given by $1/r$
 
.. math::

    u(x) = \sum_{j=1}^{N} \frac{c_{j}}{4\pi\|x-x_{j}\|} - v_{j} \cdot \nabla \left( \frac{1}{4\pi\|x-x_{j}\|}\right)   

where $c_{j}$ are the charge densities
$v_{j}$ are the dipole orientation vectors, and
$x_{j}$ are the source locations.
When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
from the sum.

.. code:: matlab
   
   function [U] = lfmm3d(eps,srcinfo,pg,targ,pgt)

Wrapper for fast multipole implementation for Laplace N-body
interactions.

Args:

-  eps: double   
      precision requested
-  srcinfo: structure
      structure containing sourceinfo
   
   *  srcinfo.sources: double(3,n)    
         source locations, $x_{j}$
   *  srcinfo.nd: integer
         number of charge/dipole vectors (optional, 
         default - nd = 1)
   *  srcinfo.charges: double(nd,n) 
         charge densities, $c_{j}$ (optional, 
         default - term corresponding to charges dropped)
   *  srcinfo.dipoles: double(nd,3,n) 
         dipole orientation vectors, $v_{j}$ (optional
         default - term corresponding to dipoles dropped) 

-  pg: integer
      | source eval flag
      | potential at sources evaluated if pg = 1
      | potenial and gradient at sources evaluated if pg=2
-  targ: double(3,nt)
      target locations, $t_{i}$ (optional)
-  pgt: integer
      | target eval flag (optional)
      | potential at targets evaluated if pgt = 1
      | potenial and gradient at targets evaluated if pgt=2  

Returns:

-  U.pot: potential at source locations, if requested, $u(x_{j})$
-  U.grad: gradient at source locations, if requested, $\nabla u(x_{j})$
-  U.pottarg: potential at target locations, if requested, $u(t_{i})$
-  U.gradtarg: gradient at target locations, if requested, $\nabla u(t_{i})$

------------------------------------------------------------------

Wrapper for direct evaluation of Laplace N-body interactions. 
Note that this wrapper only returns potentials and gradients at the
target locations.
              
.. code:: matlab
   
   function [U] = l3ddir(srcinfo,targ,pgt)

------------------------------------------------------------------

Example:

-  see ``lfmm3d_example.m``

.. container:: rttext

  `Back to top <matlab.html#mat>`__



.. _helm-mat:

Helmholtz wrappers
*******************


This subroutine computes the N-body Helmholtz
interactions and its gradients in three dimensions where 
the interaction kernel is given by $e^{ikr}/r$
 
.. math::

    u(x) = \sum_{j=1}^{N} c_{j} \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} - v_{j} \cdot \nabla \left( \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|}\right)   

where $c_{j}$ are the charge densities
$v_{j}$ are the dipole orientation vectors, and
$x_{j}$ are the source locations.
When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
from the sum.

.. code:: matlab
   
   function [U] = hfmm3d(eps,zk,srcinfo,pg,targ,pgt)

Wrapper for fast multipole implementation for Helmholtz N-body
interactions.

Args:

-  eps: double   
      precision requested
-  zk: complex
      Helmholtz parameter, k
-  srcinfo: structure
      structure containing sourceinfo
   
   *  srcinfo.sources: double(3,n)    
         source locations, $x_{j}$
   *  srcinfo.nd: integer
         number of charge/dipole vectors (optional, 
         default - nd = 1)
   *  srcinfo.charges: complex(nd,n) 
         charge densities, $c_{j}$ (optional, 
         default - term corresponding to charges dropped)
   *  srcinfo.dipoles: complex(nd,3,n) 
         dipole orientation vectors, $v_{j}$ (optional
         default - term corresponding to dipoles dropped) 

-  pg: integer
      | source eval flag
      | potential at sources evaluated if pg = 1
      | potenial and gradient at sources evaluated if pg=2
-  targ: double(3,nt)
      target locations, $t_{i}$ (optional)
-  pgt: integer
      | target eval flag (optional)
      | potential at targets evaluated if pgt = 1
      | potenial and gradient at targets evaluated if pgt=2  

Returns:

-  U.pot: potential at source locations, if requested, $u(x_{j})$
-  U.grad: gradient at source locations, if requested, $\nabla u(x_{j})$
-  U.pottarg: potential at target locations, if requested, $u(t_{i})$
-  U.gradtarg: gradient at target locations, if requested, $\nabla u(t_{i})$

------------------------------------------------------------------

Wrapper for direct evaluation of Helmholtz N-body interactions.
Note that this wrapper only returns potentials and gradients at the
target locations.
              
.. code:: matlab
   
   function [U] = h3ddir(zk,srcinfo,targ,pgt)

------------------------------------------------------------------

Example:

-  see ``hfmm3d_example.m``

.. container:: rttext

  `Back to top <matlab.html#mat>`__


.. _stok-mat:

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

.. code:: matlab
   
   function [U] = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg)

Wrapper for fast multipole implementation for Stokes N-body
interactions.

Args:

-  eps: double   
      precision requested
-  srcinfo: structure
      structure containing sourceinfo
   
   *  srcinfo.sources: double(3,n)    
         source locations, $x_{j}$
   *  srcinfo.nd: integer
         number of charge/dipole vectors (optional, 
         default - nd = 1)
   *  srcinfo.stoklet: double(nd,3,n) 
         Stokeslet densities, $\sigma_{j}$ (optional, 
         default - term corresponding to Stokeslet dropped)
   *  srcinfo.strslet: double(nd,3,n) 
         Stresslet densities, $\mu_{j}$ (optional
         default - term corresponding to stresslet dropped) 
   *  srcinfo.strsvec: double(nd,3,n) 
         Stresslet orientiation vectors, $\nu_{j}$ (optional
         default - term corresponding to stresslet dropped) 
   *  srcinfo.rotlet: double(nd,3,n) 
         Rotlet densities, $\mu^{r}_{j}$ (optional
         default - term corresponding to rotlet dropped) 
   *  srcinfo.rotvec: double(nd,3,n) 
         Rotlet orientiation vectors, $\nu^{r}_{j}$ (optional
         default - term corresponding to rotlet dropped) 
   *  srcinfo.doublet: double(nd,3,n) 
         Doublet densities, $\mu^{d}_{j}$ (optional
         default - term corresponding to doublet dropped) 
   *  srcinfo.doubvec: double(nd,3,n) 
         Doublet orientiation vectors, $\nu^{d}_{j}$ (optional
         default - term corresponding to doublet dropped) 


-  ifppreg: integer
      | source eval flag
      | potential at sources evaluated if ifppreg = 1
      | potential and pressure at sources evaluated if ifppreg=2
      | potential, pressure and gradient at sources evaluated if ifppreg=3
-  targ: double(3,nt)
      target locations, $t_{i}$ (optional)
-  ifppregtarg: integer
      | target eval flag (optional)
      | potential at targets evaluated if ifppregtarg = 1
      | potential and pressure at targets evaluated if ifppregtarg = 2 
      | potential, pressure and gradient at targets evaluated if ifppregtarg = 3

Returns:

-  U.pot: velocity at source locations if requested
-  U.pre: pressure at source locations if requested
-  U.grad: gradient of velocity at source locations if requested
-  U.pottarg: velocity at target locations if requested
-  U.pretarg: pressure at target locations if requested
-  U.gradtarg: gradient of velocity at target locations if requested

------------------------------------------------------------------

Wrapper for direct evaluation of Stokes N-body interactions. 
Note that this wrapper only returns potentials and gradients at the
target locations.
              
.. code:: matlab
   
   function [U] = st3ddir(srcinfo,targ,ifppregtarg)

------------------------------------------------------------------

Example:

-  see ``stfmm3d_example.m``

.. container:: rttext

  `Back to top <matlab.html#mat>`__



.. _em-mat:

Maxwell wrappers
*******************


This subroutine computes the N-body Maxwell
interactions, its curl and its divergence in three dimensions
given by
 
.. math::

    E(x) = \sum_{j=1}^{N} \nabla \times \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} M_{j} + \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} J_{j} +  \nabla \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} \rho_{j}       

where $M_{j}$ are the magnetic current densities,
$J_{j}$ are the electric current densities, 
$\rho_{j}$ are the electric charge densities, and
$x_{j}$ are the source locations.
When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
from the sum.

.. code:: matlab
   
   function [U] = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE)

Wrapper for fast multipole implementation for Maxwell N-body
interactions.
Note that this wrapper only returns fields, divergences, and curls at the
target locations.

Args:

-  eps: double   
      precision requested
-  zk: complex
      Wavenumber, k
-  srcinfo: structure
      structure containing sourceinfo
   
   *  srcinfo.sources: double(3,n)    
         source locations, $x_{j}$
   *  srcinfo.nd: integer
         number of charge/dipole vectors (optional, 
         default - nd = 1)
   *  srcinfo.h_current: complex(nd,3,n) 
         Magnetic current densities, $M_{j}$ (optional,
         default - term corresponding to magnetic current dropped) 
   *  srcinfo.e_current: complex(nd,3,n) 
         Electric current densities, $J_{j}$ (optional,
         default - term corresponding to electric current dropped) 
   *  srcinfo.e_charge: complex(nd,n) 
         Electric charge densities, $\rho_{j}$ (optional, 
         default - term corresponding to electric charge dropped)

-  targ: double(3,nt)
      target locations, $t_{i}$ 
-  ifE: integer
      E is returned at the target locations if ifE = 1
-  ifcurlE: integer
      curl E is returned at the target locations if ifcurlE = 1
-  ifdivE: integer
      div E is returned at the target locations if ifdivE = 1

Returns:

-  U.E: E field defined above at target locations if requested $(E(t_{j}))$
-  U.curlE: curl of E field at target locations if requested $(\nabla \times E(t_{j}))$
-  U.divE: divergence of E at target locations if requested $(\nabla \cdot E(t_{j}))$

------------------------------------------------------------------

Wrapper for direct evaluation of Maxwell N-body interactions.
Note that this wrapper only returns fields, divergences, and curls at the
target locations.
              
.. code:: matlab
   
   function [U] = em3ddir(zk,srcinfo,targ,ifE,ifcurlE,ifdivE)

------------------------------------------------------------------

Example:

-  see ``emfmm3d_example.m``

.. container:: rttext

  `Back to top <matlab.html#mat>`__

