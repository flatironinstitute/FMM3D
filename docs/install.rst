Installation
============

Obtaining FINUFFT
*****************

Go to the github page https://github.com/flatironinstitute/FMM3D and
follow instructions (eg see the green button).


Dependencies
************

This library is currently supported for unix/linux
and also tested on Mac OSX. 

For the basic libraries

* Fortran compiler, such as ``gfortran`` packaged with GCC
* GNU make

Optional:

* for matlab/octave wrappers: MATLAB, or octave and its development libraries
* for building new matlab/octave wrappers (experts only): ``mwrap``
* for the python wrappers you will need ``python3`` and ``pip3``. 


Tips for installing dependencies on various operating systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On a Fedora/CentOS linux system, these dependencies can be installed as 
follows::

  sudo yum install make gcc gcc-c++ gcc-gfortran libgomp octave octave-devel

then see below for ``mwrap``.

.. note::

   we are not exactly sure how to install python and pip using yum

On Ubuntu linux (assuming python3 as opposed to python)::

  sudo apt-get install make build-essential gfortran python3 python3-pip octave liboctave-dev

On Mac OSX:

Make sure you have ``make`` installed, eg via XCode.

Install gcc, for instance using pre-compiled binaries from
http://hpc.sourceforge.net/

(Note: we are not exactly sure how to install python3 and pip3 on mac)

Currently in Mac OSX, ``make lib`` fails to make the shared object library (.so);
however the static (.a) library is of reasonable size and works fine.


Installing MWrap
----------------

This is not needed for most users.
`MWrap <http://www.cs.cornell.edu/~bindel/sw/mwrap>`_
is a very useful MEX interface generator by Dave Bindel.
Make sure you have ``flex`` and ``bison`` installed.
Download version 0.33 or later from http://www.cs.cornell.edu/~bindel/sw/mwrap, un-tar the package, cd into it, then::
  
  make
  sudo cp mwrap /usr/local/bin/

Compilation
***********

We first describe compilation for default options 
(double precision, openmp) via GCC.
If you have a nonstandard unix environment (eg a Mac) or want to change 
the compiler, then place your compiler and linking options in a new file 
``make.inc``.
For example such files see ``make.inc.*``. See ``makefile`` for what can be overridden.

Compile and do a rapid (less than 1-second) test of FINUFFT via::

  make test

This should compile the main libraries then run tests which should report 
zero crashes and zero fails. 

Use ``make perftest`` for larger FMM tests taking 15-30 seconds.

Run ``make`` without arguments for full list of possible make tasks.

Note that the library includes the C and fortran interfaces
defined in ``src/fmm3d_c.h`` and ``fortran/fmm3d_f.h`` respectively.
If there is an error in testing on a standard set-up,
please file a bug report as a New Issue at https://github.com/flatironinstitute/FMM3D/issues

Custom library compilation options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may want to make the library for other data types. 
However, single-threaded vs multithreaded are
built with the same name, so you will have to move them to other
locations, or build a 2nd copy of the repo, if you want to keep both
versions.

You *must* do at least ``make objclean`` before changing openmp options.

**Single precision**: append ``PREC=SINGLE`` to the make task.
Single-precision saves half the RAM, and increases
speed slightly (<20%). The  C++, C, and fortran demos are all tested in
single precision. However, it will break matlab, octave, python interfaces.

**Single-threaded**: append ``OMP=OFF`` to the make task.


Building examples and wrappers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``make examples`` to compile and run the examples for calling from Fortran.

The ``examples`` and ``test`` directories are good places to see usage 
examples.

``make c`` to compile and run the C wrappers and examples.


Building the python wrappers
****************************

First make sure you have python3 and pip3 (or python and pip) installed. 

You may then do ``make python3`` which calls
pip3 for the install then runs some tests. 
An additional test you could do is::

  python3 run_speed_tests.py


A few words about python environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There can be confusion and conflicts between various versions of python and installed packages. It is therefore a very good idea to use virtual environments. Here's a simple way to do it (after installing python-virtualenv)::

  Open a terminal
  virtualenv -p /usr/bin/python3 env1
  . env1/bin/activate

Now you are in a virtual environment that starts from scratch. All pip installed packages will go inside the env1 directory. (You can get out of the environment by typing ``deactivate``)
