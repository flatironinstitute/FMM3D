.. FMM3D documentation master file, created by
   sphinx-quickstart on Wed Nov  1 16:19:13 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Flatiron Institute Fast multipole methods in three dimensions (FMM3D)
====================================================

.. image:: logo.png
    :width: 45%
.. image:: spreadpic.png
    :width: 54%
	    
`FMM3D <https://github.com/flatironinstitute/FMM3D>`_ 
is a set of libraries to compute efficiently $N-$ body interactions 
for Laplace, Helmholtz, Stokes, and Maxwell kernels 
to a specified precision, in three dimensions,
on a multi-core shared-memory machine.
The library has a very simple interface, 
is written in Fortran (using OpenMP),
and has wrappers to C, MATLAB, and python.
As an example, given $M$ arbitrary points $y_j \in \mathbb{R}^{3}$ 
and complex numbers $c_j$, with $j=1,\dots,M$, and 
$N$ arbitrary points $x_{k} \in \mathhbb{R}^{3}$, the Laplace FMM
evaluates the $N$ complex numbers

.. math:: u_k = \sum_{j=1}^M \frac{c_j}{\| x_{k} - y_{j} \|} ~, 
   \qquad \mbox{ for } \; k=1,2,\ldots N ~.
   :label: lapcp

The $y_j$ can be interpreted as source locations, $c_j$
as charge strengths, and $u_k$ then as the potential at
target location $x_{k}$.




Such N-body interactions are needed in many applications in 
science and engineering, including molecular dynmaics, astrophysics, 
rheology, and numerical partial differential equations.
The naive CPU effort to evaluate :eq:`lapcp` is $O(NM)$.
The library approximates :eq:`lapcp` to a requested relative precision
$\epsilon$ with linear effort $O((M+N) \log (1/\epsilon))$.

The FMM relies on compressing the interactions between well-separated 
cluster of source and target points at a heirarchy of scales using
analytic outgoing, incoming, and plane-wave 
expansions of the interaction kernel and associated translation
operators. 
This library is the a modified version of the FMM3D-library, with the
following additions:

#. Use of plane wave expansions for diagonalizing the outgoing to
incoming translation operators
#. Vectorizing the FMM, to apply the same kernel with same source
and target locations on multiple densities.

For sources and targets distributed in the volume, this code is 4 times
faster than the previous generation on a single CPU core, and for
sources and targets distributed on a surface, this code is 2 times
faster.  


.. note::

   For very small repeated problems (less than 10000 input and output points),
   users should also consider a dense matrix-matrix multiplication against
   the $N-$body interaction matrix using BLAS3 (eg ZGEMM).
   This is currently work in progress.

   
.. toctree::
   :maxdepth: 2
	   
   install
   math
   dirs
   usage
   matlab
   pythoninterface
   juliainterface
   related
   issues
   ackn
   refs
   

   
