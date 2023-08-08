% FMM3D: MATLAB/Octave wrappers to 3D fast multipole methods.
%        Center for Computational Mathematics, Flatiron Institute.
%
% Functions:
%   lfmm3d   - FMM in 3D for Laplace (electrostatic) kernels.
%   l3ddir   - Direct (slow) 3D Laplace kernel sums (reference for LFMM3D).
%   hfmm3d   - FMM in 3D for Helmholtz (acoustic frequency-domain) kernel.
%   h3ddir   - Direct (slow) 3D Helmholtz kernel sums (reference for HFMM3D).
%   emfmm3d  - FMM in 3D for Maxwell (frequency-domain electromagnetic) kernels.
%   em3ddir  - Direct (slow) 3D Maxwell kernel sums (reference for EMFMM3D).
%   stfmm3d  - FMM in 3D for Stokes (viscous fluid hydrodynamic) kernels.
%   st3ddir  - Direct (slow) 3D Stokes kernel sums (reference for STFMM3D).
%
% For tester driver scripts see:
%   lfmm3dTest
%   hfmm3dTest
%   emfmm3dTest
%   stfmm3dTest
%   stfmm3dPerfTest
%
% For examples of use see:
%   lfmm3d_example
%   lfmm3d_big_example
%   hfmm3d_example
%   hfmm3d_big_example
%   emfmm3d_example
%   stfmm3d_example
%
% Legacy codes/interfaces (from the Gimbutas-Greengard CMCL 2012 library):
%   lfmm3dpart       - Laplace particle targ FMM in R^3.
%   l3dpartdirect    - Laplace interactions in R^3, direct (slow) evaluation.
%   hfmm3dpart       - Helmholtz particle targ FMM in R^3.
%   h3dpartdirect    - Helmholtz interactions in R^3, direct (slow) evaluation.
%   lfmm3dLegacyTest - Test Laplace particle FMMs in R^3
%   hfmm3dLegacyTest - Test Helmholtz particle FMMs in R^3
