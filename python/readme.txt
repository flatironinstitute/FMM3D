To build the fortran module:
python setup.py build

To install in development mode:
python setup.py develop

To install in the system (some pitfalls):
pip install -e . (currently not working)

Alternate way to install (some pitfalls):
python setup.py develop


On compilation you can import the package fmm3dpy
which has 2 main subroutines

hfmm3d - for computing N-body Helmholtz interactions
lfmm3d - for computing N-body Laplace interactions

Examples for both of these are demonstrated in
hfmmexample.py and lfmmexample.py


