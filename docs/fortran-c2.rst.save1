Fortran and C examples
======================

Helmholtz FMM
--------------

The Helmholtz FMM evaluates the following potential and it's
gradient 

.. math::
 
    u(x) = \sum_{j=1}^{N} c_{j} G_{k}(x-x_{j})  - \nabla  G_{k}(x-x_{j}) \cdot v_{j} \, .
     
Here $x_{j}$ are the source locations, $c_{j}$ are the 
charge strengths and $v_{j}$ are the dipole strengths.
     

The subroutine names take the following form:

.. highlights::  
   
   hfmm3dpart<1><2><3>

    - <1>: Collection of `x' where $u$ and it's gradient

        - stos: Evaluate $u$ and it's gradient at the source locations $x_{i}$ 
        - stot: Evaluate $u$ and it's gradient at $t_{i}$, a collection of target locations specified  by the user.
        - stost: Evaluate $u$ and it's gradient at both source and target locations $x_{i}$ and $t_{i}$.

    - <2>: kind of interaction (charges/dipoles/both). The charge interactions are given by $c_{j} G_{k}(x-x_{j})$, and the dipole interactions are given by $-\nabla G_{k}(x-x_{j}) \cdot v_{j}$

        - c: charges
        - d: dipoles
        - cd: charges + dipoles
 
    - <3>: Flag for evaluating potential or potential + gradient

        - p: on output only $u$ is evaluated
        - g: on output both $u$ and it's gradient are evaluated

These are all the single density routines. To get a vectorized version 
of any of the routines use

.. epigraph::  
   
   "<subroutine name>_vec"


.. include:: fortrandocs_helm_opt2.raw
