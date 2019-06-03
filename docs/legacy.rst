.. _lega: 

FMMLIB3D Legacy interfaces
===========================

The current version of the FMM codes are backward compatible with the
previous version of this library: `FMMLIB3D
<https::/github.com/zgimbutas/fmmlib3d>`_.
On this page, we refer to these wrappers as the legacy wrappers.

.. note::
    
    The field associated with the potential returned in FMMLIB3D is negative of the gradient of the potential.

-  `Laplace wrappers <legacy.html#lap-lega>`__
-  `Helmholtz wrappers <legacy.html#helm-lega>`__

.. _lap-lega:

Laplace
########
The legacy Fortran Laplace wrappers are contained 
in ``src/Laplace/lfmm3dwrap_legacy.f`` and the legacy MATLAB Laplace
wrappers are contained in ``matlab/lfmm3dpart.m`` and
``matlab/l3dpartdirect.m``. 

Currently we have interfaces for the following four Fortran wrappers and
two matlab wrappers:

-   Two self evaluation wrappers (:ref:`lfmm3dpart`)
-   The main fmm wrapper and direct evaluation wrapper in fortran (:ref:`lfmm3dparttarg`)
-   :ref:`lap-lega-mat`


.. note::
   In the Laplace wrappers for FMMLIB3D, the charge strengths, dipole
   strengths, potentials, and fields are complex numbers as opposed
   to real numbers for the rest of the library.

.. note::
   lfmm3dpartself and lfmm3dpart are identical subroutines except for
   their names.
   

.. _lfmm3dpart:

lfmm3dpart and lfmm3dpartself
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Evaluation points: Sources
- Interaction kernel: Charges/Dipoles/Charges+Dipoles
- Outputs requested: Potential/Fields/Potential+Fields


.. code:: fortran

   subroutine lfmm3dpart(ier,iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld)

.. code:: fortran

   subroutine lfmm3dpartself(ier,iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld)

This subroutine evaluates the potential/field/potential
and field


  .. math::

      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - d_{j}\left(v_{j} \cdot \nabla \left( \frac{1}{\|x-x_{j}\|}\right) \right)

at the source locations $x=x_{j}$. When $x=x_{m}$, the term 
corresponding to $x_{m}$ is dropped from the sum.

Input arguments:

  -    iprec: integer
          | precision flag
          | iprec=-2 => tolerance = 0.5d0
          | iprec=-1 => tolerance = 0.5d-1
          | iprec=0 => tolerance = 0.5d-2
          | iprec=1 => tolerance = 0.5d-3
          | iprec=2 => tolerance = 0.5d-6
          | iprec=3 => tolerance = 0.5d-9
          | iprec=4 => tolerance = 0.5d-12
  -    nsource: integer
          Number of sources
  -    source: double precision(3,nsource)
          Source locations, $x_{j}$
  -    ifcharge: integer
          | charge computation flag 
          | ifcharge =1 => include charge contribution, otherwise do not
  -    charge: double complex(nsource)
          Charge strengths, $c_{j}$
  -    ifdipole: integer
          | dipole computation flag 
          | ifdipole =1 => include dipole contribution, otherwise do not
  -    dipstr: double complex(nsource)
          Dipole strengths, $d_{j}$
  -    dipvec: double precision(3,nsource)
          Dipole orientation vectors, $v_{j}$
  -    ifpot: integer
          | potential flag 
          | ifpot =1 => compute potential, otherwise do not
  -    iffld: integer
          | Field flag 
          | iffld =1 => compute field, otherwise do not

Output arguments:

  -    ier: integer
          error code, currently unused
  -    pot: double complex(nsource)
          Potential at source locations, if requested, $u(x_{j})$
  -    fld: double complex(3,nsource)
          Field at source locations, if requested, -$\nabla u(x_{j})$

.. container:: rttext

  `Back to Laplace legacy wrappers <legacy.html#lap-lega>`__

.. container:: rttext

  `Back to top <legacy.html#lega>`__


.. _lfmm3dparttarg:

lfmm3dparttarg and l3dpartdirect
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- Evaluation points: Sources/Targets/Sources+targets
- Interaction kernel: Charges/Dipoles/Charges+Dipoles
- Outputs requested: Potential/Fields/Potential+Fields


.. code:: fortran

   subroutine lfmm3dparttarg(ier,iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld,ntarg,targ,ifpottarg,pottarg,iffldtarg,fldtarg)

This subroutine evaluates the potential/field/potential
and field


  .. math::

      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - d_{j}\left(v_{j} \cdot \nabla \left( \frac{1}{\|x-x_{j}\|}\right) \right)

at the source locations $x=x_{j}$/target locations $x=t_{j}$/
source and target locations. When $x=x_{m}$, the term 
corresponding to $x_{m}$ is dropped from the sum.

Input arguments:

  -    iprec: integer
          | precision flag
          | iprec=-2 => tolerance = 0.5d0
          | iprec=-1 => tolerance = 0.5d-1
          | iprec=0 => tolerance = 0.5d-2
          | iprec=1 => tolerance = 0.5d-3
          | iprec=2 => tolerance = 0.5d-6
          | iprec=3 => tolerance = 0.5d-9
          | iprec=4 => tolerance = 0.5d-12
  -    nsource: integer
          Number of sources
  -    source: double precision(3,nsource)
          Source locations, $x_{j}$
  -    ifcharge: integer
          | charge computation flag 
          | ifcharge =1 => include charge contribution, otherwise do not
  -    charge: double complex(nsource)
          Charge strengths, $c_{j}$
  -    ifdipole: integer
          | dipole computation flag 
          | ifdipole =1 => include dipole contribution, otherwise do not
  -    dipstr: double complex(nsource)
          Dipole strengths, $d_{j}$
  -    dipvec: double precision(3,nsource)
          Dipole orientation vectors, $v_{j}$
  -    ifpot: integer
          | potential flag 
          | ifpot =1 => compute potential, otherwise do not
  -    iffld: integer
          | Field flag 
          | iffld =1 => compute field, otherwise do not
  -    ntarg: integer
          Number of targets 
  -    targ: double precision(3,ntarg)
          Source locations, $x_{j}$
  -    ifpottarg: integer
          | target potential flag 
          | ifpottarg =1 => compute potential, otherwise do not
  -    iffldtarg: integer
          | target field flag 
          | iffldtarg =1 => compute field, otherwise do not

Output arguments:

  -    ier: integer
          error code, currently unused
  -    pot: double complex(nsource)
          Potential at source locations, if requested, $u(x_{j})$
  -    fld: double complex(3,nsource)
          Field at source locations, if requested, -$\nabla u(x_{j})$
  -    pottarg: double complex(ntarg)
          Potential at target locations, if requested, $u(t_{j})$
  -    fld: double complex(3,ntarg)
          Field at source locations, if requested, -$\nabla u(t_{j})$

---------------------------------------------------------

Wrapper for direct evaluation of Laplace N-body interactions.

.. code:: fortran

   subroutine l3dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld,ntarg,targ,ifpottarg,pottarg,iffldtarg,fldtarg)

------------------------------------------------------------

Example:

-   see ``examples/lfmm3d_legacy_example.f``. The corresponding makefile
    is ``examples/lfmm3d_legacy_example.make``.


.. container:: rttext

  `Back to Laplace legacy wrappers <legacy.html#lap-lega>`__

.. container:: rttext

  `Back to top <legacy.html#lega>`__

.. _lap-lega-mat:

MATLAB wrappers
~~~~~~~~~~~~~~~~

- `matlab/lfmm3dpart.m`
- Evaluation points: Sources/Targets/Sources+targets
- Interaction kernel: Charges/Dipoles/Charges+Dipoles
- Outputs requested: Potential/Fields/Potential+Fields


.. code:: matlab

   function [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)

This subroutine evaluates the potential/field/potential
and field


  .. math::

      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - d_{j}\left(v_{j} \cdot \nabla \left( \frac{1}{\|x-x_{j}\|}\right) \right)

at the source locations $x=x_{j}$/target locations $x=t_{j}$/
source and target locations. When $x=x_{m}$, the term 
corresponding to $x_{m}$ is dropped from the sum.

See :ref:`lfmm3dparttarg` for a detailed description of input and 
output arguments. 
The output pot,pottarg,fld,fldtarg are contained in the output structure
U.

The function can be called in 4 different ways

.. code:: matlab

   function [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec)
   function [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld)
   function [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ)
   function [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)

The default argument for ifpot,iffld,ifpottarg,iffldtarg is 1, the
defaults for ntarg is 0, and targ is zeros(3,1)

---------------------------------------------------------

Wrapper for direct evaluation of Laplace N-body interactions.

- `matlab/l3dpartdirect.m`

.. code:: matlab

   function [U]=l3dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)

The function can be called in 4 different ways

.. code:: matlab

   function [U]=l3dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec)
   function [U]=l3dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld)
   function [U]=l3dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ)
   function [U]=l3dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)

------------------------------------------------------------

Example:

-   see ``matlab/test_lfmm3dpart_direct.m``. 

.. container:: rttext

  `Back to Laplace legacy wrappers <legacy.html#lap-lega>`__

.. container:: rttext

  `Back to top <legacy.html#lega>`__


.. _helm-lega:

Helmholtz
###################

The legacy Fortran Helmholtz wrappers are contained 
in ``src/Helmholtz/hfmm3dwrap_legacy.f`` and the legacy MATLAB Helmholtz
wrappers are contained in ``matlab/hfmm3dpart.m`` and
``matlab/h3dpartdirect.m``. 

Currently we have interfaces for the following four Fortran wrappers and
two matlab wrappers:

-   Two self evaluation wrappers (:ref:`hfmm3dpart`)
-   The main fmm wrapper and direct evaluation wrapper in fortran (:ref:`hfmm3dparttarg`)
-   :ref:`helm-lega-mat`


.. note::
   hfmm3dpartself and hfmm3dpart are identical subroutines except for
   their names.
   

.. _hfmm3dpart:

hfmm3dpart and lfmm3dpartself
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Evaluation points: Sources
- Interaction kernel: Charges/Dipoles/Charges+Dipoles
- Outputs requested: Potential/Fields/Potential+Fields


.. code:: fortran

   subroutine hfmm3dpart(ier,iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld)

.. code:: fortran

   subroutine hfmm3dpartself(ier,iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld)

This subroutine evaluates the potential/field/potential
and field


  .. math::

      u(x) = \sum_{j=1}^{N} c_{j} \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} - d_{j}\left(v_{j} \cdot \nabla \left( \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right) \right)

at the source locations $x=x_{j}$. When $x=x_{m}$, the term 
corresponding to $x_{m}$ is dropped from the sum.

Input arguments:

  -    iprec: integer
          | precision flag
          | iprec=-2 => tolerance = 0.5d0
          | iprec=-1 => tolerance = 0.5d-1
          | iprec=0 => tolerance = 0.5d-2
          | iprec=1 => tolerance = 0.5d-3
          | iprec=2 => tolerance = 0.5d-6
          | iprec=3 => tolerance = 0.5d-9
          | iprec=4 => tolerance = 0.5d-12
  -    zk: double complex
          Helmholtz parameter, k
  -    nsource: integer
          Number of sources
  -    source: double precision(3,nsource)
          Source locations, $x_{j}$
  -    ifcharge: integer
          | charge computation flag 
          | ifcharge =1 => include charge contribution, otherwise do not
  -    charge: double complex(nsource)
          Charge strengths, $c_{j}$
  -    ifdipole: integer
          | dipole computation flag 
          | ifdipole =1 => include dipole contribution, otherwise do not
  -    dipstr: double complex(nsource)
          Dipole strengths, $d_{j}$
  -    dipvec: double precision(3,nsource)
          Dipole orientation vectors, $v_{j}$
  -    ifpot: integer
          | potential flag 
          | ifpot =1 => compute potential, otherwise do not
  -    iffld: integer
          | Field flag 
          | iffld =1 => compute field, otherwise do not

Output arguments:

  -    ier: integer
          error code, currently unused
  -    pot: double complex(nsource)
          Potential at source locations, if requested, $u(x_{j})$
  -    fld: double complex(3,nsource)
          Field at source locations, if requested, -$\nabla u(x_{j})$

.. container:: rttext

  `Back to Helmholtz legacy wrappers <legacy.html#helm-lega>`__

.. container:: rttext

  `Back to top <legacy.html#lega>`__


.. _hfmm3dparttarg:

hfmm3dparttarg and h3dpartdirect
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- Evaluation points: Sources/Targets/Sources+targets
- Interaction kernel: Charges/Dipoles/Charges+Dipoles
- Outputs requested: Potential/Fields/Potential+Fields


.. code:: fortran

   subroutine hfmm3dparttarg(ier,iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld,ntarg,targ,ifpottarg,pottarg,iffldtarg,fldtarg)

This subroutine evaluates the potential/field/potential
and field


  .. math::

      u(x) = \sum_{j=1}^{N} c_{j} \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} - d_{j}\left(v_{j} \cdot \nabla \left( \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right) \right)

at the source locations $x=x_{j}$/target locations $x=t_{j}$/
source and target locations. When $x=x_{m}$, the term 
corresponding to $x_{m}$ is dropped from the sum.

Input arguments:

  -    iprec: integer
          | precision flag
          | iprec=-2 => tolerance = 0.5d0
          | iprec=-1 => tolerance = 0.5d-1
          | iprec=0 => tolerance = 0.5d-2
          | iprec=1 => tolerance = 0.5d-3
          | iprec=2 => tolerance = 0.5d-6
          | iprec=3 => tolerance = 0.5d-9
          | iprec=4 => tolerance = 0.5d-12
  -    zk: double complex
          Helmholtz parameter, k
  -    nsource: integer
          Number of sources
  -    source: double precision(3,nsource)
          Source locations, $x_{j}$
  -    ifcharge: integer
          | charge computation flag 
          | ifcharge =1 => include charge contribution, otherwise do not
  -    charge: double complex(nsource)
          Charge strengths, $c_{j}$
  -    ifdipole: integer
          | dipole computation flag 
          | ifdipole =1 => include dipole contribution, otherwise do not
  -    dipstr: double complex(nsource)
          Dipole strengths, $d_{j}$
  -    dipvec: double precision(3,nsource)
          Dipole orientation vectors, $v_{j}$
  -    ifpot: integer
          | potential flag 
          | ifpot =1 => compute potential, otherwise do not
  -    iffld: integer
          | Field flag 
          | iffld =1 => compute field, otherwise do not
  -    ntarg: integer
          Number of targets 
  -    targ: double precision(3,ntarg)
          Source locations, $x_{j}$
  -    ifpottarg: integer
          | target potential flag 
          | ifpottarg =1 => compute potential, otherwise do not
  -    iffldtarg: integer
          | target field flag 
          | iffldtarg =1 => compute field, otherwise do not

Output arguments:

  -    ier: integer
          error code, currently unused
  -    pot: double complex(nsource)
          Potential at source locations, if requested, $u(x_{j})$
  -    fld: double complex(3,nsource)
          Field at source locations, if requested, -$\nabla u(x_{j})$
  -    pottarg: double complex(ntarg)
          Potential at target locations, if requested, $u(t_{j})$
  -    fld: double complex(3,ntarg)
          Field at source locations, if requested, -$\nabla u(t_{j})$

---------------------------------------------------------

Wrapper for direct evaluation of Helmholtz N-body interactions.

.. code:: fortran

   subroutine h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld,ntarg,targ,ifpottarg,pottarg,iffldtarg,fldtarg)

------------------------------------------------------------

Example:

-   see ``examples/hfmm3d_legacy_example.f``. The corresponding makefile
    is ``examples/hfmm3d_legacy_example.make``.


.. container:: rttext

  `Back to Helmholtz legacy wrappers <legacy.html#helm-lega>`__

.. container:: rttext

  `Back to top <legacy.html#lega>`__

.. _helm-lega-mat:

MATLAB wrappers
~~~~~~~~~~~~~~~~

- `matlab/hfmm3dpart.m`
- Evaluation points: Sources/Targets/Sources+targets
- Interaction kernel: Charges/Dipoles/Charges+Dipoles
- Outputs requested: Potential/Fields/Potential+Fields


.. code:: matlab

   function [U]=hfmm3dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)

This subroutine evaluates the potential/field/potential
and field


  .. math::

      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - d_{j}\left(v_{j} \cdot \nabla \left( \frac{1}{\|x-x_{j}\|}\right) \right)

at the source locations $x=x_{j}$/target locations $x=t_{j}$/
source and target locations. When $x=x_{m}$, the term 
corresponding to $x_{m}$ is dropped from the sum.

See :ref:`hfmm3dparttarg` for a detailed description of input and 
output arguments. 
The output pot,pottarg,fld,fldtarg are contained in the output structure
U.

The function can be called in 4 different ways

.. code:: matlab

   function [U]=hfmm3dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec)
   function [U]=hfmm3dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld)
   function [U]=hfmm3dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ)
   function [U]=hfmm3dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)

The default argument for ifpot,iffld,ifpottarg,iffldtarg is 1, the
defaults for ntarg is 0, and targ is zeros(3,1)

---------------------------------------------------------

Wrapper for direct evaluation of Helmholtz N-body interactions.

- `matlab/h3dpartdirect.m`

.. code:: matlab

   function [U]=h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)

The function can be called in 4 different ways

.. code:: matlab

   function [U]=h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec)
   function [U]=h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld)
   function [U]=h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ)
   function [U]=h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)

------------------------------------------------------------

Example:

-   see ``matlab/test_hfmm3dpart_direct.m``. 

.. container:: rttext

  `Back to Helmholtz legacy wrappers <legacy.html#helm-lega>`__

.. container:: rttext

  `Back to top <legacy.html#lega>`__



