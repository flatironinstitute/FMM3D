  This folder contains the header files for the c wrappers,
  some utility files, and a few examples. 

  -  The Helmholtz wrappers are contained in ``hfmm3d_c.h``
  -  The Laplace wrappers are contained in ``lfmm3d_c.h``

  - The Helmholtz FMM example files are ``hfmm3d_example.c`` and
  ``hfmm3d_vec_example.c`` (for vectorized routines).
  To run these examples, run their corresponding makefiles, 
  ``hfmm3d_example.make`` and ``hfmm3d_vec_example.make``.

  - The Laplace FMM example files are ``lfmm3d_example.c`` and
  ``lfmm3d_vec_example.c`` (for vectorized routines).
  To run these examples, run their corresponding makefiles, 
  ``lfmm3d_example.make`` and ``lfmm3d_vec_example.make``.

  - In order to run any of the makefiles, you will need to have compiled
  the static library in the main folder by running
    
        make lib

    Then to run the makefile, run

        make -f "<makefile>"
 



