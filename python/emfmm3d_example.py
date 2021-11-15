#!/usr/bin/env python

import fmm3dpy as fmm
import numpy as np


#
#  This is a sample code to demonstrate how to use
#  the fmm libraries
#
#
# Source to target, magnetic current only, output field example
#
n = 20000
nd = 1
sources = np.random.uniform(0,1,(3,n))

nt = 4000
targ = np.random.uniform(0,1,(3,nt))
eps = 10**(-5)

zk = 1.1 + 1j*0
h_current = np.random.uniform(0,1,(3,n)) + 1j*np.random.uniform(0,1,(3,n))
out = fmm.emfmm3d(eps=eps,zk=zk, sources=sources,h_current=h_current,targets=targ,ifE = 1)

#
# Example 2: source to target, electric current + charge, multiple densities, field + curl
#
#

nd = 3
zk = 1.1 + 1j*0
e_charge = np.random.uniform(0,1,(nd,n)) + 1j*np.random.uniform(0,1,(nd,n))
e_current = np.random.uniform(0,1,(nd,3,n)) + 1j*np.random.uniform(0,1,(nd,3,n))
out = fmm.emfmm3d(eps=eps,zk=zk, sources=sources,e_charge=e_charge,\
    e_current=e_current,targets=targ,nd=nd,ifE = 1, ifcurlE = 1)
