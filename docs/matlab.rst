.. _mat:

MATLAB
=======

The MATLAB interface has four callable subroutines:

*  `Helmholtz wrappers <matlab.html#helm-mat>`__: Fast multipole implementation (hfmm3d) and direct evaluation (h3ddir) for Helmholtz N-body interactions
*  `Laplace wrappers <matlab.html#lap-mat>`__: Fast multipole implementation (lfmm3d) and direct evaluation (l3ddir) for Laplace N-body interactions


.. _helm-mat:

Helmholtz wrappers
*******************


This subroutine computes the N-body Helmholtz
interactions and its gradients in three dimensions where 
the interaction kernel is given by $e^{ikr}/r$
 
.. math::

    u(x) = \sum_{j=1}^{N} c_{j} \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} - v_{j} \cdot \nabla \left( \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)   

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

-  see ``hfmmexample.m``

.. container:: rttext

  `Back to top <matlab.html#mat>`__


.. _lap-mat:

Laplace wrappers
*******************


This subroutine computes the N-body Laplace
interactions and its gradients in three dimensions where 
the interaction kernel is given by $1/r$
 
.. math::

    u(x) = \sum_{j=1}^{N} \frac{c_{j}}{\|x-x_{j}\|} - v_{j} \cdot \nabla \left( \frac{1}{\|x-x_{j}\|}\right)   

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

-  see ``lfmmexample.m``

.. container:: rttext

  `Back to top <matlab.html#mat>`__


