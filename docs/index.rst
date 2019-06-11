.. FMM3D documentation master file, created by
   sphinx-quickstart on Wed Nov  1 16:19:13 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Fast multipole methods in three dimensions (FMM3D)
==================================================

.. image:: FMM-logo.png
    :width: 60%
    :align: center
	    
`FMM3D <https://github.com/flatironinstitute/FMM3D>`_ 
is a set of libraries to compute N-body interactions 
for Laplace, and Helmholtz 
to a specified precision, in three dimensions,
on a multi-core shared-memory machine.
The library is written in Fortran,
and has wrappers to C, MATLAB, and Python.
As an example, given $M$ arbitrary points $y_j \in \mathbb{R}^{3}$ 
with corresponding real numbers $c_j$, and 
$N$ arbitrary points $x_{j} \in \mathbb{R}^{3}$, the Laplace FMM
evaluates the $N$ real numbers

.. math:: u_{\ell} = \sum_{j=1}^M \frac{c_j}{\| x_{\ell} - y_{j}\|} ~, 
   \qquad \mbox{ for } \; \ell=1,2,\ldots N ~.
   :label: lapcp

The $y_j$ can be interpreted as source locations, $c_j$
as charge strengths, and $u_{\ell}$ as the resulting potential at
target location $x_{\ell}$.

Such N-body interactions are needed in many applications in 
science and engineering, including molecular dynamics, astrophysics, 
rheology, and numerical solution of partial differential equations.
The naive CPU effort to evaluate :eq:`lapcp` is $O(NM)$.
The library approximates :eq:`lapcp` to a requested relative precision
$\epsilon$ with linear effort $O((M+N) \log^{3} (1/\epsilon))$.

The FMM relies on compressing the interactions between well-separated 
clusters of source and target points at a hierarchy of scales using
analytic outgoing, incoming, and plane-wave 
expansions of the interaction kernel and associated translation
operators. 
This library is an improved version of the `FMMLIB3D <https://github.com/zgimbutas/fmmlib3d>`_
software, Copyright (C) 2010-2012: Leslie Greengard and Zydrunas Gimbutas, released under the 
BSD license. The major changes are the following:

-  The use of plane wave expansions for diagonalizing the outgoing to incoming translation operators
-  Vectorization of the FMM, to apply the same kernel with same source and target locations on multiple 
   strength vectors.
-  A redesign of the adaptive tree data structure

For sources and targets distributed in the volume, this code is 4 times
faster than the previous generation on a single CPU core, and for
sources and targets distributed on a surface, this code is 2 times
faster.  

.. note::
   
   The plane wave expansions for the Helmholtz FMMs have only been incorporated
   for low frequency problems (problems less than 32 wavelengths in size in each dimension), 
   and real Helmholtz parameter. 


.. note::

   For very small repeated problems (less than 1000 input and output points),
   users should also consider a dense matrix-matrix multiplication against
   the $N-$body interaction matrix using BLAS3 (eg DGEMM,ZGEMM).
   This is currently work in progress.

   
.. toctree::
   :maxdepth: 2
	   
   install
   math
   fortran-c
   matlab
   python
   legacy
   ackn
   ref
   

   
