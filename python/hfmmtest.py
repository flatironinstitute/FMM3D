#!/usr/bin/env python

import fmm3dpy as fmm
import numpy as np


n = 2000
nd = 1
sources = np.random.uniform(0,1,(3,n))
eps = 10**(-5)

zk = 1.1 + 1j*0
charges = np.random.uniform(0,1,n) + 1j*np.random.uniform(0,1,n)
out = fmm.hfmm3d(eps=eps,zk=zk, sources=sources,charges=charges,pg=1)
print(out.pot[0:5])

