This folder contains the MATLAB wrappers, legacy wrappers for FMMLIB3D, 
and precompiled mex interfaces for MacOSX and Linux.

We strongly recommend recompiling the mex interfaces (see 
https://fmm3d.readthedocs.io/en/latest/install.html for detailed installation 
instructions).

Example codes for demonstrating the Helmholtz and Laplace interfaces are
``hfmm3d_example.m`` and ``lfmm3d_example.m``.

Example codes for demonstrating the legacy Helmholtz and Laplace wrappers are
``test_hfmm3dpart_direct.m`` and ``test_lfmm3dpart_direct.m``.

After editing [stfmm3d.f](../src/Stokes/stfmm3d.f), run from the repository root:

```sh
make matlab
```

This recompiles the Fortran object, updates the static library, and relinks
the MATLAB MEX binaries. The existing gateway is unchanged, so neither
`make mex` nor `mwrap` is needed. If MATLAB's `mex` is not on your path,
set `MEX=/path/to/MATLAB/bin/mex` in the command or in `make.inc`.
Then, in MATLAB:

```matlab
clear fmm3d
addpath('matlab') % from the repository root
stfmm3d_laplaceTest
```

[stfmm3d.m](./stfmm3d.m) performs recentering and velocity assembly in Fortran.
[stfmm3d_laplaceTest.m](./stfmm3d_laplaceTest.m) contains the MATLAB reference
composition and reports agreement, complete-call timings and translation
accuracy. It compares the updated Fortran routine with MATLAB composition
using recentering and without it (`laplace_unshifted`). All three methods
use the Fortran Laplace FMM. Edit the parameters at the top of the script
to change the problem sizes and repetitions. Benchmark the old binary
before rebuilding if a before/after comparison is needed.
