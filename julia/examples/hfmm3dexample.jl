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

zk = 1.1 + im*0
charges = rand(n) + im*rand(n)
vals = hfmm3d(eps,zk,sources,charges=charges,pg=1)


# sample with a vector of densities, sources to 
# sources and targets, dipole interactions, 
# potential and gradietns

nd = 3
zk = 1.1 + im*0
nt = 1870
targ = rand(3,nt)
dipvecs = rand(nd,3,n) + im*rand(nd,3,n)
vals = hfmm3d(eps,zk,sources,dipvecs=dipvecs,
             targets=targ,nd=nd,pg=2,pgt=2)
