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
n = 200000
nd = 1
sources = np.random.uniform(0,1,(3,n))
eps = 10**(-5)

charges = np.random.uniform(0,1,n) 
out = fmm.lfmm3d(eps=eps,sources=sources,charges=charges,pg=1)


# sample with a vector of densities, sources to 
# sources and targets, dipole interactions, 
# potential and gradietns

nd = 3
nt = 1870
targ = np.random.uniform(0,1,(3,nt))
dipvecs = np.random.uniform(0,1,(nd,3,n)) 
out = fmm.lfmm3d(eps=eps,sources=sources,dipvec=dipvecs,\
    targets=targ,nd=nd,pg=2,pgt=2)
