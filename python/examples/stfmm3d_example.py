#!/usr/bin/env python

import fmm3dpy as fmm
import numpy as np


#
#  This is a sample code to demonstrate how to use
#  the fmm libraries
#

# sample with one density, sources to sources,
# stokeslet interactions, and potential+pressure+gradient
#
n = 20000
nd = 1
sources = np.random.uniform(0,1,(3,n))
eps = 10**(-5)

stoklet = np.random.uniform(0,1,(3,n)) 
out = fmm.stfmm3d(eps=eps,sources=sources,stoklet=stoklet,ifppreg=3)


# sample with a vector of densities, sources to 
# sources and targets, stokeslet + stresslet interactions, 
# potential + pressure

nd = 2
nt = 1870
targ = np.random.uniform(0,1,(3,nt))
stoklet = np.random.uniform(0,1,(nd,3,n)) 
strsvec = np.random.uniform(0,1,(nd,3,n)) 
strslet = np.random.uniform(0,1,(nd,3,n)) 
out = fmm.stfmm3d(eps=eps,sources=sources,stoklet=stoklet,\
    strslet=strslet,strsvec=strsvec,\
    targets=targ,nd=nd,ifppreg=2,ifppregtarg=2)
