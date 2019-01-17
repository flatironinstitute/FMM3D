import fmm3d_fortran as fmm
import numpy as np

def sample_high_level_call(zk=1+0.01j,eps=5e-3):
  #zk=1.1+0.1*1j
  ns=1000
  sources=np.random.uniform(0,1,(3,ns))
  charge=np.random.uniform(0,1,(ns))+np.random.uniform(0,1,(ns))*1j
  print('here is the part that might work:')
  area=fmm.hfmm3dpartstoscp(eps,zk,sources,charge,ns)
  print('it might have worked.')
  return area
