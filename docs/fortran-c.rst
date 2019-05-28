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
     

The subroutine names take the following form:

.. highlights::  
   
   "lfmm3d_a_b_c"

    - <1>: evaluation points. Collection of `x' where $u$ and its gradient

        - stos: Evaluate $u$ and its gradient at the source locations $x_{i}$ 
        - stot: Evaluate $u$ and its gradient at $t_{i}$, a collection of target locations specified  by the user.
        - stost: Evaluate $u$ and its gradient at both source and target locations $x_{i}$ and $t_{i}$.

    - <2>: kind of interaction (charges/dipoles/both). The charge interactions are given by $c_{j}/\|x-x_{j}\| $, and the dipole interactions are given by $-v_{j} \cdot \nabla (1/\|x-x_{j}\|)$

        - c: charges
        - d: dipoles
        - cd: charges + dipoles
 
    - <3>: Flag for evaluating potential or potential + gradient

        - p: on output only $u$ is evaluated
        - g: on output both $u$ and its gradient are evaluated

These are all the single density routines. To get a vectorized version 
of any of the routines use

.. epigraph::  
   
   "<subroutine name>_vec"

.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__



List of interfaces
******************
 
- Evaluation points: Sources  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`lscp`)          
    - Gradient  (:ref:`lscg`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`lsdp`)          
    - Gradient  (:ref:`lsdg`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`lscdp`)         
    - Gradient  (:ref:`lscdg`)         


- Evaluation points: Targets  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`ltcp`)          
    - Gradient  (:ref:`ltcg`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`ltdp`)          
    - Gradient  (:ref:`ltdg`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`ltcdp`)         
    - Gradient  (:ref:`ltcdg`)         

- Evaluation points: Sources + Targets  
                                        
  - Interaction Type: Charges           
                                        
    - Potential (:ref:`lstcp`)          
    - Gradient  (:ref:`lstcg`)          
                                       
  - Interaction Type: Dipoles           
                                       
    - Potential (:ref:`lstdp`)          
    - Gradient  (:ref:`lstdg`)          
                                       
  - Interaction Type: Charges + Dipoles 
                                       
    - Potential (:ref:`lstcdp`)         
    - Gradient  (:ref:`lstcdg`)         

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
     

The subroutine names take the following form:

.. highlights::  
   
   hfmm3dpart<1><2><3>

    - <1>: evaluation points. Collection of `x' where $u$ and its gradient

        - stos: Evaluate $u$ and its gradient at the source locations $x_{i}$ 
        - stot: Evaluate $u$ and its gradient at $t_{i}$, a collection of target locations specified  by the user.
        - stost: Evaluate $u$ and its gradient at both source and target locations $x_{i}$ and $t_{i}$.

    - <2>: kind of interaction (charges/dipoles/both). The charge interactions are given by $c_{j} e^{ik\|x-x_{j}\|}/\|x-x_{j}\|$, and the dipole interactions are given by $-v_{j} \cdot \nabla (e^{ik\|x-x_{j}\|}/\|x-x_{j}\|)$

        - c: charges
        - d: dipoles
        - cd: charges + dipoles
 
    - <3>: Flag for evaluating potential or potential + gradient

        - p: on output only $u$ is evaluated
        - g: on output both $u$ and its gradient are evaluated

These are all the single density routines. To get a vectorized version 
of any of the routines use

.. epigraph::  
   
   "<subroutine name>_vec"

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

All of the above fortran routines can be called from c using

.. epigraph:
   
   "<fortran subroutine name>"_("<calling sequence>")


Note that all the variables in the calling sequence must be passed
as pointers. For sample code, see examples in 'c/'.  

.. container:: rttext

  `Back to top <fortran-c.html#fcexmp>`__

