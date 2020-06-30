Installation
============

Obtaining FMM3D
***************

The source code can be downloaded from https://github.com/flatironinstitute/FMM3D 


Dependencies
************

This library is supported for unix/linux, Mac OSX, and Windows.

For the basic libraries

* Fortran compiler, such as ``gfortran`` packaged with GCC
* GNU make

Optional:

* for building Python wrappers you will need ``python3`` and ``pip3`` 
* for building standard MATLAB wrappers: MATLAB
* for modifying MATLAB wrappers (experts only): ``mwrap``

Quick install instructions
*********************************************

Make sure you have dependencies installed, and `cd` into your FMM3D
directory. 

-  For linux, run ``make install``.
-  For Mac OSX, run ``cp make.inc.macos.gnu make.inc`` followed by ``make install``.
-  For Windows, run ``cp make.inc.windows.mingw make.inc`` followed by ``make install`` 

This should compile the static library
in ``lib-static/``, the dynamic library in ``lib/`` and copy the dynamic 
library to ``$(HOME)/lib`` on Linux,  to
``/usr/local/lib`` on Mac OSX, and to ``C:\lib`` on Windows.
The location of the default installation directory can be changed by
running::

    make install PREFIX=(INSTALL_DIR)


In order to link against the dynamic library, you will have to update
the ``PATH`` environment variable on Windows, ``LD_LIBRARY_PATH`` environment
variable on Linux and ``DYLD_LIBRARY_PATH`` environment variable on Mac OSX
to the installation directory.
You may then link to the FMM library using the ``-lfmm3d`` option.

.. note :: 
   On MacOSX, /usr/local/lib is included by default in the
   DYLD_LIBRARY_PATH.


To verify successful compilation of the program, run ``make test``
which compiles some fortran test drivers in ``test/`` linked against
the static library, after which it
runs the test programs. The last 14 lines of the terminal output should be::

   cat print_testreshelm.txt
   Successfully completed 5 out of 5 tests in helmrouts3d testing suite
   Successfully completed 18 out of 18 tests in hfmm3d testing suite
   Successfully completed 6 out of 6 tests in hfmm3d scale testing suite
   Successfully completed 18 out of 18 tests in hfmm3d vec testing suite
   Successfully completed 1 out of 1 tests in helm3d_mps testing suite
   cat print_testreslap.txt
   Successfully completed 5 out of 5 tests in laprouts3d testing suite
   Successfully completed 18 out of 18 tests in lfmm3d testing suite
   Successfully completed 2 out of 2 tests in lfmm3d scale testing suite
   Successfully completed 18 out of 18 tests in lfmm3d vec testing suite
   rm print_testreshelm.txt
   rm print_testreslap.txt


To verify successful installation of the program, and the correct
setting for environment variables, run ``make test-dyn`` which compiles
some fortran test drivers in ``test/`` linked against the dynamic
library, after which it runs teh test prgram. The output ofshould be the
same as above.

.. note ::
   By default, ``make install`` creates the easy-to-install version of the library. To
   compile the library in its high-performance mode, append
   ``FAST_KER=ON`` to the make task. For instance ``make install`` should be replaced by 
   ``make install FAST_KER=ON``. See :ref:`custom-install` for
   other options.
   

If ``make test`` fails, see more detailed instructions below.

If ``make test-dyn`` fails with an error about not finding
``-llibfmm3d_dll`` or ``-lfmm3d`` or ``libfmm3d`` make sure that the
appropriate environment variables have been set. If it fails with other
issues, see more detailed instructions below.

Type ``make`` to see a list of other build options (language
interfaces, etc). Please see `Fortran and C interfaces <fortran-c.html>`__ and look in
``examples/`` for sample drivers.

If there is an error in testing on a standard set-up,
please file a bug report as a New Issue at https://github.com/flatironinstitute/FMM3D/issues

.. _custom-install:

Custom library compilation options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the (default) easy-to-install version,
the library is compiled  without using the optimized direct evaluation kernels.

In order to disable multi-threading, append ``OMP=OFF`` to the make task.

In order to use the optimized direct evaluation kernels (this
automatically turns on multithreading as well), append ``FAST_KER=ON`` to
the make task. This option is currently *not* supported on Windows.


All of these different libraries are
built with the same name, so you will have to move them to other
locations, or build a 2nd copy of the repo, if you want to keep both
versions.

You *must* do at least ``make objclean`` before changing to the openmp
/fast direct kernel evaluation options.


Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``make examples`` to compile and run the examples for calling from Fortran.
*  ``make c-examples`` to compile and run the examples for calling from C.

The ``examples`` directory is a good place to see usage 
examples for Fortran.
There are three sample Fortran drivers  
for both the Laplace and Helmholtz FMMs,
one which demonstrates the use of FMMs, one which demonstrates
the use of vectorized FMMs, and one which demonstrates the 
use of the same calling sequence as FMMLIB3D - so that legacy codes
are backward compatible with `FMMLIB3D <https://github.com/zgimbutas/fmmlib3d>`_.

The sample drivers for the Laplace FMM are
``lfmm3d_example.f``, ``lfmm3d_vec_example.f``, and
``lfmm3d_legacy_example.f``, and the corresponding makefiles
are ``lfmm3d_example.make``, ``lfmm3d_vec_example.make``, and
``lfmm3d_legacy_example.make``. These demonstrate how to link
to the dynamic library ``libfmm3d.so``.
The analogous Helmholtz drivers are ``hfmm3d_example.f``,
``hfmm3d_vec_example.f``, and ``hfmm3d_legacy_example.f``.
The corresponding makefiles are ``hfmm3d_example.make``, 
``hfmm3d_vec_example.make``, and ``hfmm3d_legacy_example.make``.


Analogous C sample drivers can be found in ``c/``.


Building Python wrappers
****************************

First make sure you have python (version 3 or higher) and pip installed. 

You may then execute ``make python`` (after copying over the
operating system specific make.inc.* file to make.inc) which calls
pip for the install and then runs some tests.

To rerun the tests, you may run ``pytest`` in ``python/`` 
or alternatively run ``python python/test_hfmm.py`` and 
``python python/test_lfmm.py``.

See ``python/hfmmexample.py`` and ``python/lfmmexample.py`` to see
usage examples for the Python wrappers.

.. note::
   On windows, you will need to update ``distutils.cfg`` located in 
   ``(PYTHON_INSTALL_DIR)\Lib\distutils`` and set it to::

       [build]
       compiler=mingw32

       [build_ext]
       compiler=mingw32

   which forces python to use the mingw compiler for building its
   modules. In case you wish to revert to using VC/C++ for building python
   modules, make sure to update distutils.cfg appropriately.


A few words about Python environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There can be confusion and conflicts between various versions of Python and installed packages. It is therefore a very good idea to use virtual environments. Here's a simple way to do it (after installing python-virtualenv)::

  Open a terminal
  virtualenv -p /usr/bin/python3 env1
  . env1/bin/activate

Now you are in a virtual environment that starts from scratch. All pip installed packages will go inside the env1 directory. (You can get out of the environment by typing ``deactivate``)



Building the MATLAB wrappers
****************************

First make sure you have MATLAB installed. 

Then run ``make matlab`` (after copying over the operating
system specific make.inc.* file to make.inc) which links the .m files to
the .c file in the matlab folder.

To run tests, you can run ``matlab test_hfmm3d.m`` and 
``matlab test_lfmm3d.m`` and it should return with $0$ crashes.

Example codes for demonstrating the Helmholtz and Laplace
interfaces are ``hfmm3d_example.m`` and ``lfmm3d_example.m``.


Tips for installing dependencies
**********************************

On Ubuntu linux
~~~~~~~~~~~~~~~~

On Ubuntu linux (assuming python3 as opposed to python)::

  sudo apt-get install make build-essential gfortran  


On Fedora/CentOS linux
~~~~~~~~~~~~~~~~~~~~~~~~

On a Fedora/CentOS linux system, these dependencies can be installed as 
follows::

  sudo yum install make gcc gcc-c++ gcc-gfortran libgomp 

.. _mac-inst:

On Mac OSX
~~~~~~~~~~~~~~~~~~~~~~~~

First setup Homebrew as follows. If you don't have Xcode, install
Command Line Tools by opening a terminal (from /Applications/Utilities/)
and typing::

  xcode-select --install

Then install Homebrew by pasting the installation command from
https://brew.sh

Then do::
  
  brew install gcc 


On Windows
~~~~~~~~~~~~~~~

Download 64 bit mingw (available `here <http://mingw-w64.org/doku.php>`_). 
Follow the install instructions and append to the environment variable ``PATH`` the
location of the bin directory of your mingw installation.

Download  and install ``make`` for windows 
(Available `here <http://gnuwin32.sourceforge.net/packages/make.htm>`_).

Download and install ``git`` for windows
(Available `here <https://git-scm.com/download/win>`_).

Tips for installing optional dependencies
******************************************

Installing python and pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Ubuntu linux
##################

::

  sudo apt-get install python3 python3-pip


On Mac OSX
############

Make sure you have homebrew installed. See `Tips for installing dependencies -> On Mac OSX <install.html#mac-inst>`__ 

::
  
  brew install python3


On Windows
###########

Download and install python3.7 from python.org.

Configuring MATLAB
~~~~~~~~~~~~~~~~~~~

On Windows
############

Update ``MINGW_LPATH`` in ``make.inc.windows.mingw`` to point to the
appropriate installation directory (it should be the one within the
``gcc`` folder).

To setup mingw as the C compiler on MATLAB run ``configuremingw.p``
(which can be downloaded from 
`here <https://www.mathworks.com/matlabcentral/answers/uploaded_files/88639/configuremingw.p>`_)
and choose the mingw directory. To verify successful setup run ``mex
-setup`` from matlab and it should be configured to compile with mingw.


Installing MWrap
~~~~~~~~~~~~~~~~~~

If you make any changes to the 
fortran code, you will need to regenerate the .c files
from the .mw files for which mwrap is required.
This is not needed for most users.
`MWrap <http://www.cs.cornell.edu/~bindel/sw/mwrap>`_
is a very useful MEX interface generator by Dave Bindel.

Make sure you have ``flex`` and ``bison`` installed.
Download version 0.33.5 or later from https://github.com/zgimbutas/mwrap, un-tar the package, cd into it, then::
  
  make
  sudo cp mwrap /usr/local/bin/


