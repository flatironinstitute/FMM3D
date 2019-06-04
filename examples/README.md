  This folder contains the examples for the Fortran wrappers. 

  - The Helmholtz FMM example files are ``hfmm3d_example.f``
  ``hfmm3d_vec_example.f`` (for vectorized routines), and
  ``hfmm3d_legacy_example.f`` (for the legacy routines).
  To run these examples, run their corresponding makefiles, 
  ``hfmm3d_example.make``, ``hfmm3d_vec_example.make``, and
  ``hfmm3d_legacy_example.make``.

  - The Helmholtz FMM example files are ``lfmm3d_example.f``
  ``lfmm3d_vec_example.f`` (for vectorized routines), and
  ``lfmm3d_legacy_example.f`` (for the legacy routines).
  To run these examples, run their corresponding makefiles, 
  ``lfmm3d_example.make``, ``lfmm3d_vec_example.make``, and
  ``lfmm3d_legacy_example.make``.

  - In order to run any of the makefiles, you will need to have compiled
  the static library in the main folder by running
    
        make lib

    Then to run the makefile, run

        make -f "<makefile>"
 



