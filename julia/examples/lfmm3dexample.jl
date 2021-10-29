using FMM3D
using Random
using LinearAlgebra

#  This is a sample code to demonstrate how to use
#  the fmm libraries
#

# sample with one density, sources to sources,
# charge interactions, and potential only
#
n = 2000
nd = 1
sources = rand(3,n)
eps = 10.0^(-5)

charges = rand(n) 
vals = lfmm3d(eps,sources,charges=charges,pg=1)


# sample with a vector of densities, sources to 
# sources and targets, dipole interactions, 
# potential and gradietns

nd = 3
nt = 1870
targ = rand(3,nt)
dipvecs = rand(nd,3,n) 
vals = lfmm3d(eps,sources,dipvecs=dipvecs,
             targets=targ,nd=nd,pg=2,pgt=2)
