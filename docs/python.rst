.. _pyt:

Python
=======

The Python interface has four callable subroutines:

*  `Helmholtz wrappers <python.html#helm-pyt>`__: Fast multipole implementation (hfmm3d) and direct evaluation (h3ddir) for Helmholtz N-body interactions
*  `Laplace wrappers <python.html#lap-pyt>`__: Fast multipole implementation (lfmm3d) and direct evaluation (l3ddir) for Laplace N-body interactions


.. _helm-pyt:

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

-  see ``hfmmexample.py``

.. container:: rttext

  `Back to top <python.html#pyt>`__


.. _lap-pyt:

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

-  see ``lfmmexample.py``

.. container:: rttext

  `Back to top <python.html#pyt>`__

