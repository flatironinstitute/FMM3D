.. _fcexmp:

Fortran and C interfaces
========================

-  :ref:`lap`
-  :ref:`helm`
-  :ref:`stokes`
-  :ref:`maxwell`
-  :ref:`cinter`


.. _lap:

Laplace FMM
------------

The Laplace FMM evaluates the following potential and its
gradient

.. math::
 
    u(x) = \sum_{j=1}^{N} \frac{c_{j}}{4\pi\|x-x_{j}\|}   - v_{j} \cdot \nabla  \left( \frac{1}{4\pi\|x-x_{j}\|} \right)  \, .
     
Here $x_{j}$ are the source locations, $c_{j}$ are the 
charge strengths and $v_{j}$ are the dipole strengths, 
and the collection of $x$ at which the potential
and its gradient are evaluated are referred to as the
evalution points.
     
There are 18 different Fortran wrappers for the Laplace FMM 
to account for collection of evaluation points (sources only, 
targets only, sources+targets), interaction kernel (charges only,
dipoles only, charges + dipoles), output request (potential,
potential+gradient).

For example, the subroutine to evaluate the potential and gradient, at a collection
of targets $t_{i}$ due to a collection of charges is::

   lfmm3d_t_c_g


In general, the subroutine names take the following form::

   lfmm3d_<eval-pts>_<int-ker>_<out>



- <eval-pts>: evaluation points. Collection of `x` where $u$ and its gradient is to be evaluated

    - s: Evaluate $u$ and its gradient at the source locations $x_{i}$ 
    - t: Evaluate $u$ and its gradient at $t_{i}$, a collection of target locations specified  by the user.
    - st: Evaluate $u$ and its gradient at both source and target locations $x_{i}$ and $t_{i}$.

- <int-ker>: kernel of interaction (charges/dipoles/both). The charge interactions are given by $c_{j}/4\pi\|x-x_{j}\| $, and the dipole interactions are given by $-v_{j} \cdot \nabla (1/4\pi\|x-x_{j}\|)$

    - c: charges
    - d: dipoles
    - cd: charges + dipoles
 
- <out>: Flag for evaluating potential or potential + gradient

    - p: on output only $u$ is evaluated
    - g: on output both $u$ and its gradient are evaluated
    - h: on output $u$, its gradient and its hessian are evaluated

These are all the single density routines. To get a vectorized version 
of any of the routines use::

   <subroutine name>_vec

.. note::

   For the vectorized subroutines, the charge strengths, dipole
   strengths, potentials, and gradients are interleaved as opposed to
   provided in a sequential manner. For example for three sets of charge
   strengths, they should be stored as $c_{1,1}, c_{2,1}, c_{3,1},
   c_{1,2}, c_{2,2},c_{3,2} \ldots c_{1,N}, c_{2,N}, c_{3,N}$. 


Example drivers:

-   ``examples/lfmm3d_example.f``. The corresponding makefile is
    ``examples/lfmm3d_example.make``
-   ``examples/lfmm3d_vec_example.f``. The corresponding makefile is
    ``examples/lfmm3d_vec_example.make``

   

.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__



List of interfaces
******************
 
- Evaluation points: Sources  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`lscp`)          
    - Gradient  (:ref:`lscg`)          
    - Hessian   (:ref:`lsch`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`lsdp`)          
    - Gradient  (:ref:`lsdg`)          
    - Hessian   (:ref:`lsdh`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`lscdp`)         
    - Gradient  (:ref:`lscdg`)         
    - Hessian   (:ref:`lscdh`)          


- Evaluation points: Targets  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`ltcp`)          
    - Gradient  (:ref:`ltcg`)          
    - Hessian   (:ref:`ltch`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`ltdp`)          
    - Gradient  (:ref:`ltdg`)          
    - Hessian   (:ref:`ltdh`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`ltcdp`)         
    - Gradient  (:ref:`ltcdg`)         
    - Hessian   (:ref:`ltcdh`)          

- Evaluation points: Sources + Targets  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`lstcp`)          
    - Gradient  (:ref:`lstcg`)          
    - Hessian   (:ref:`lstch`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`lstdp`)          
    - Gradient  (:ref:`lstdg`)          
    - Hessian   (:ref:`lstdh`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`lstcdp`)         
    - Gradient  (:ref:`lstcdg`)         
    - Hessian   (:ref:`lstcdh`)          

.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__


.. include:: fortrandocs_lap.raw


.. _helm:

Helmholtz FMM
--------------

The Helmholtz FMM evaluates the following potential and its
gradient 

.. math::
 
    u(x) = \sum_{j=1}^{N} \frac{c_{j} e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|}   - v_{j} \cdot \nabla  \left( \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} \right) \, .
     
Here $x_{j}$ are the source locations, $c_{j}$ are the 
charge strengths and $v_{j}$ are the dipole strengths, 
and the collection of $x$ at which the potential
and its gradient are evaluated are referred to as the
evalution points.
     

There are 18 different Fortran wrappers for the Helmholtz FMM 
to account for collection of evaluation points (sources only, 
targets only, sources+targets), interaction kernel (charges only,
dipoles only, charges + dipoles), output request (potential,
potential+gradient).

For example, the subroutine to evaluate the potential and gradient, at a collection
of targets $t_{i}$ due to a collection of charges is::

   hfmm3d_t_c_g

In general, the subroutine names take the following form::

   hfmm3d_<eval-pts>_<int-ker>_<out>



- <eval-pts>: evaluation points. Collection of `x` where $u$ and its gradient is to be evaluated

    - s: Evaluate $u$ and its gradient at the source locations $x_{i}$ 
    - t: Evaluate $u$ and its gradient at $t_{i}$, a collection of target locations specified  by the user.
    - st: Evaluate $u$ and its gradient at both source and target locations $x_{i}$ and $t_{i}$.

- <int-ker>: kernel of interaction (charges/dipoles/both). The charge interactions are given by $c_{j}/4\pi\|x-x_{j}\| $, and the dipole interactions are given by $-v_{j} \cdot \nabla (1/4\pi\|x-x_{j}\|)$

    - c: charges
    - d: dipoles
    - cd: charges + dipoles
 
- <out>: Flag for evaluating potential or potential + gradient

    - p: on output only $u$ is evaluated
    - g: on output both $u$ and its gradient are evaluated

These are all the single density routines. To get a vectorized version 
of any of the routines use::

   <subroutine name>_vec

.. note::

   For the vectorized subroutines, the charge strengths, dipole
   strengths, potentials, and gradients are interleaved as opposed to
   provided in a sequential manner. For example for three sets of charge
   strengths, they should be stored as $c_{1,1}, c_{2,1}, c_{3,1},
   c_{1,2}, c_{2,2},c_{3,2} \ldots c_{1,N}, c_{2,N}, c_{3,N}$. 


Example drivers:

-   ``examples/hfmm3d_example.f``. The corresponding makefile is
    ``examples/hfmm3d_example.make``
-   ``examples/hfmm3d_vec_example.f``. The corresponding makefile is
    ``examples/hfmm3d_vec_example.make``

   
.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__


List of interfaces
******************
 
- Evaluation points: Sources  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`hscp`)          
    - Gradient  (:ref:`hscg`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`hsdp`)          
    - Gradient  (:ref:`hsdg`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`hscdp`)         
    - Gradient  (:ref:`hscdg`)         


- Evaluation points: Targets  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`htcp`)          
    - Gradient  (:ref:`htcg`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`htdp`)          
    - Gradient  (:ref:`htdg`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`htcdp`)         
    - Gradient  (:ref:`htcdg`)         

- Evaluation points: Sources + Targets  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`hstcp`)          
    - Gradient  (:ref:`hstcg`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`hstdp`)          
    - Gradient  (:ref:`hstdg`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`hstcdp`)         
    - Gradient  (:ref:`hstcdg`)         


.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__


.. include:: fortrandocs_helm.raw


.. _stokes:

Stokes FMM
------------


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

The Stokes FMM evaluates the following velocity, its gradient
and the associated pressure

.. math::

    u(x) = \sum_{m=1}^{N} \mathcal{G}^{\textrm{stok}}(x,x_{j}) \sigma_{j}  + \nu_{j} \cdot \mathcal{T}^{\textrm{stok}}(x,x_{j}) \cdot \mu_{j} +
                          \nu^{r}_{j} \cdot \mathcal{R}^{\textrm{stok}}(x,x_{j}) \cdot \mu^{r}_{j} - \mu^{r}_{j} \cdot \mathcal{R}^{\textrm{stok}}(x,x_{j}) \cdot \nu^{r}_{j} + \\
                          \nu^{d}_{j} \cdot \mathcal{R}^{\textrm{stok}}(x,x_{j}) \cdot \mu^{d}_{j} - \mu^{d}_{j} \cdot \mathcal{R}^{\textrm{stok}}(x,x_{j}) \cdot \nu^{d}_{j} +
                          \nu^{d}_{j} \cdot \mathcal{D}^{\textrm{stok}}(x,x_{j}) \cdot \mu^{d}_{j} \, .

Here $x_{j}$ are the source locations,
$\sigma_{j}$ are the Stokeslet densities,
$\nu_{j}$ are the stresslet orientation vectors, $\mu_{j}$ 
are the stresslet densities,
$\nu^{r}_{j}$ are the rotlet orientation vectors, $\mu^{r}_{j}$ 
are the rotlet densities,
$\nu^{d}_{j}$ are the doublet orientation vectors, $\mu^{d}_{j}$ 
are the doublet densities,
and the locations $x$ 
at which the velocity and its gradient are evaluated are referred to
as the evaluation points.

Unlike the Laplace and Helmholtz FMM, currently we  
have only the guru interface for the Stokes FMM  (for both the single density
and the vectorized density cases)
with appropriate flags for including or excluding the 
stokeslet/stresslet term in the interaction, and flags
for computing velocity/velocity and pressure/velocity, pressure, and
gradients at the evaluation points.


.. code::

   subroutine stfmm3d(nd,eps,nsource,source,ifstoklet,stoklet,ifstrslet,strslet,strsvec,ifrotlet,rotlet,rotvec,ifdoublet,doublet,doubvec,ifppreg,pot,pre,grad,ntarg,targ,ifppregtarg,pottarg,pretarg,gradtarg,ier) 

Input arguments:

  -    nd: integer
          Number of densities
  -    eps: double precision
          Precision requested
  -    nsource: integer
          Number of sources
  -    source: double precision(3,nsource)
          Source locations, $x_{j}$
  -    ifstoklet: integer
          Flag for including Stokeslet ($\sigma_{j}$) term in interaction kernel
          Stokeslet term will be included if ifstoklet
          = 1
  -    stoklet: double precision(nd,3,nsource)
          Stokeslet strengths, $\sigma_{j}$
  -    ifstrslet: integer
          Flag for including Stresslet ($\mu_{j},\nu_{j}$) term in interaction kernel
          Stresslet term will be included if ifstrslet
          = 1
  -    strslet: double precision(nd,3,nsource)
          Stresslet strengths, $\mu_{j}$
  -    strsvec: double precision(nd,3,nsource)
          Stresslet orientation vectors, $\nu_{j}$
  -    ifrotlet: integer
          Flag for including Rotlet ($\mu^{r}_{j},\nu^{r}_{j}$) term in interaction kernel
          Rotlet term will be included if ifrotlet
          = 1
  -    rotlet: double precision(nd,3,nsource)
          Rotlet strengths, $\mu^{r}_{j}$
  -    rotvec: double precision(nd,3,nsource)
          Rotlet orientation vectors, $\nu^{r}_{j}$
  -    ifdoublet: integer
          Flag for including Doublet ($\mu^{d}_{j},\nu^{d}_{j}$) term in interaction kernel
          Doublet term will be included if ifdoublet
          = 1
  -    doublet: double precision(nd,3,nsource)
          Doublet strengths, $\mu^{d}_{j}$
  -    doubvec: double precision(nd,3,nsource)
          Doublet orientation vectors, $\nu^{d}_{j}$
  -    ifppreg: integer
          | Flag for computing velocity, pressure and/or gradients at source locations
          | ifppreg = 1, compute velocity
          | ifppreg = 2, compute velocity and pressure
          | ifppreg = 3, compute veloicty, pressure and gradient
  -    ntarg: integer
          Number of targets
  -    targets: double precision (3,ntarg)
          Target locations $x$
  -    ifppregtarg: integer
          | Flag for computing velocity, pressure and/or gradients at target locations
          | ifppregtarg = 1, compute velocity
          | ifppregtarg = 2, compute velocity and pressure
          | ifppregtarg = 3, compute veloicty, pressure and gradient

Output arguments:

  -    pot: double precision (nd,3,nsource)
          Velocity at source locations if requested
  -    pre: double precision (nd,nsource)
          Pressure at source locations if requested
  -    grad: double precision (nd,3,3,nsource)
          Gradient at source locations if requested
  -    pottarg: double precision (nd,3,ntarg)
          Velocity at target locations if requested
  -    pretarg: double precision (nd,ntarg)
          Pressure at target locations if requested
  -    gradtarg: double precision (nd,3,3,ntarg)
          Gradient at target locations if requested
  -    ier: integer
          Error flag; ier=0 implies successful execution, and ier=4/8
          implies insufficient memory

Example drivers:

-   ``examples/stfmm3d_example.f``. The corresponding makefile is
    ``examples/stfmm3d_example.make``
   

.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__


.. _maxwell:

Maxwell FMM
--------------

The Maxwell FMM evaluates the following field, its curl, and its
divergence

 
.. math::

    E(x) = \sum_{j=1}^{N} \nabla \times \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} M_{j} + \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} J_{j} +  \nabla \frac{e^{ik\|x-x_{j}\|}}{4\pi\|x-x_{j}\|} \rho_{j}  \, .     

Here $x_{j}$ are the source locations,
$M_{j}$ are the magnetic current densities,
$J_{j}$ are the electric current densities, 
$\rho_{j}$ are the electric charge densities, and
the collection of $x$ at which the field, its curl and its
divergence are evalauated are referred to as the evaluation points.

Unlike the Laplace and Helmholtz FMM, currently we  
have only the guru interface for the Maxwell FMM  (for both the single density
and the vectorized density cases)
with appropriate flags for including or excluding the 
magnetic current/electric current/electric charge term in the interaction, and flags
for computing field/curl/divergence at the evaluation points.

.. code:: fortran

   subroutine emfmm3d(nd,eps,zk,ns,source,ifh_current,h_current,ife_current,e_current,ife_charge,e_charge,nt,targets,ifE,E,ifcurlE,curlE,ifdivE,divE,ier) 

Input arguments:

  -    nd: integer
          Number of densities
  -    eps: double precision
          Precision requested
  -    zk: double complex
          Wave number $k$
  -    ns: integer
          Number of sources
  -    source: double precision(3,ns)
          Source locations, $x_{j}$
  -    ifh_current: integer
          Flag for including magnetic current ($M_{j}$) term in
          interaction kernel.
          Magnetic current term will be included if ifh_current
          = 1
  -    h_current: double complex(nd,3,ns)
          Magnetic currents, $M_{j}$
  -    ife_current: integer
          Flag for including electric current ($J_{j}$) term in interaction kernel.
          Electric current term will be included if ife_current
          = 1
  -    e_current: double complex(nd,3,ns)
          Electric currents, $J_{j}$
  -    ife_charge: integer
          Flag for including electric charge ($\rho_{j}$) term in interaction kernel.
          Electric charge term will be included if ife_charge
          = 1
  -    e_charge: double complex(nd,ns)
          Electric charges, $\rho_{j}$
  -    nt: integer
          Number of targets
  -    targets: double precision (3,nt)
          Target locations $x$
  -    ifE: integer
          Flag for computing field. The field $E$ will be returned if
          ifE = 1
  -    ifcurlE: integer
          Flag for computing curl of field. $\nabla \times E$ will be returned if
          ifcurlE = 1
  -    ifdivE: integer
          Flag for computing divergence of field. $\nabla \cdot E$ will be returned if
          ifdivE = 1

Output arguments:

  -    E: double complex (nd,3,nt)
          Field at the evaluation points if requested
  -    curlE: double complex (nd,3,nt)
          Curl of field at the evaluation points if requested
  -    divE: double complex (nd,nt)
          Divergence of field at the evaluation points if requested
  -    ier: integer
          Error flag; ier=0 implies successful execution, and ier=4/8
          implies insufficient memory


Example drivers:

-   ``examples/emfmm3d_example.f``. The corresponding makefile is
    ``examples/emfmm3d_example.make``

   
.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__


List of interfaces
******************
 

.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__


.. _cinter:

C interfaces
------------

All of the above fortran routines can be called from c by including the
header ``utils.h`` and ``lfmm3d_c.h`` for Laplace FMMs or ``hfmm3d_c.h`` for
Helmholtz FMMs. 

For example, the subroutine to evaluate the potential and gradient, at a collection
of targets $t_{i}$ due to a collection of Helmholtz charges is::

   hfmm3d_t_c_g_


In general, to call a fortran subroutine from ``c`` use::

   "<fortran subroutine name>"_("<calling sequence>") 


.. note:: 
   All the variables in the calling sequence must be passed
   as pointers from ``c``. 

.. note::

   For the vectorized subroutines, the charge strengths, dipole
   strengths, potentials, and gradients are interleaved as opposed to
   provided in a sequential manner. For example for three sets of charge
   strengths, they should be stored as $c_{1,1}, c_{2,1}, c_{3,1},
   c_{1,2}, c_{2,2},c_{3,2} \ldots c_{1,N}, c_{2,N}, c_{3,N}$. 

Example drivers:

    - Laplace:

        -   ``c/lfmm3d_example.c``. The corresponding makefile is
            ``c/lfmm3d_example.make``
        -   ``c/lfmm3d_vec_example.c``. The corresponding makefile is
            ``c/lfmm3d_vec_example.make``

    - Helmholtz:

        -   ``c/hfmm3d_example.c``. The corresponding makefile is
            ``c/hfmm3d_example.make``
        -   ``c/hfmm3d_vec_example.c``. The corresponding makefile is
            ``c/hfmm3d_vec_example.make``

The Maxwell and Stokes interfaces are currently unavailable in C, and
will be made available soon.

.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__

