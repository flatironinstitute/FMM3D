.. _fcexmp:

Fortran and C interfaces
========================

-  :ref:`lap`
-  :ref:`helm`
-  :ref:`cinter`


.. _lap:

Laplace FMM
------------

The Laplace FMM evaluates the following potential and its
gradient

.. math::
 
    u(x) = \sum_{j=1}^{N} \frac{c_{j}}{\|x-x_{j}\|}   - v_{j} \cdot \nabla  \left( \frac{1}{\|x-x_{j}\|} \right)  \, .
     
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

- <int-ker>: kernel of interaction (charges/dipoles/both). The charge interactions are given by $c_{j}/\|x-x_{j}\| $, and the dipole interactions are given by $-v_{j} \cdot \nabla (1/\|x-x_{j}\|)$

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
 
    u(x) = \sum_{j=1}^{N} \frac{c_{j} e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}   - v_{j} \cdot \nabla  \left( \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} \right) \, .
     
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

- <int-ker>: kernel of interaction (charges/dipoles/both). The charge interactions are given by $c_{j}/\|x-x_{j}\| $, and the dipole interactions are given by $-v_{j} \cdot \nabla (1/\|x-x_{j}\|)$

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

   
.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__

