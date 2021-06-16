# FMM3D.jl

FMM3D.jl is a set of julia interfaces for
computing N-body interactions using the
Flatiron Institute's `FMM3D` library.

## Using FMM3D.jl

For now, FMM3D.jl can be obtained by downloading this
git repository and a recent version of julia
(1.6 or later is required). Then, from the
root directory (the parent directory of the
julia subdirectory) run the following in
julia:

```julia
 Pkg.add(url=".",subdir="julia")
 ```

You can test the package with

```julia
Pkg.test("FMM3D")
```

## Contributing

Contributions are welcome! Check out the issues
tab for feature requests/bugs.

It's pretty straightforward to write
a wrapper for any subroutine in the FMM3D library
because they are all exposed by the FMM3D_jll library
and can be executed using `ccall` (see existing wrappers
in src/FMM3D.jl for examples). Many of the
low-level routines have other uses and could
benefit from a user-friendly wrapper.

Requirements for a pull request to be approved:
* Any new wrapper should have associated tests
in test/runtests.jl.
* If the pull request is fixing a bug, that bug
should be tested in test/runtests.jl
* Any new wrapper which is exported by the module
should have documentattion. Typically, the documentation
from the original Fortran code is a good starting
point. See existing wrappers for style.
* Wrappers for the highest level routines (actual
FMMs and direct evaluation codes) should be added to
the list in the module documentation at the top of
src/FMM3D.jl. All other wrappers may be added to the
list in the documentation for the `lower_level_routs`
function in src/FMM3D.jl.