#!/usr/bin/env python

import fmm3dpy as fmm
import numpy as np


#
#  This is a sample code to demonstrate how to use
#  the fmm libraries
#

# sample with one density, sources to sources,
# charge interactions, and potential only
#
n = 2000
nd = 1
sources = np.random.uniform(0,1,(3,n))
eps = 10**(-5)

zk = 1.1 + 1j*0
charges = np.random.uniform(0,1,n) + 1j*np.random.uniform(0,1,n)
out = fmm.hfmm3d(eps=eps,zk=zk, sources=sources,charges=charges,pg=1)


# sample with a vector of densities, sources to 
# sources and targets, dipole interactions, 
# potential and gradietns

nd = 3
zk = 1.1 + 1j*0
nt = 1870
targ = np.random.uniform(0,1,(3,nt))
dipvecs = np.random.uniform(0,1,(nd,3,n)) + 1j*np.random.uniform(0,1,(nd,3,n))
out = fmm.hfmm3d(eps=eps,zk=zk, sources=sources,dipvec=dipvecs,\
    targets=targ,nd=nd,pg=2,pgt=2)
