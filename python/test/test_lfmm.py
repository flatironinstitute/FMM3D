#!/usr/bin/env python

import fmm3dpy as fmm
import numpy as np
import numpy.linalg as la

def main():
    test_lfmm()

def test_lfmm():
    ntests = 54
    testres = np.zeros(ntests)
    #
    #  This is a testing code for making sure all the 
    #  fmm routines are accessible through fmm3d.py
    #

    n = 2000
    ntest = 10
    zk = 1.1 + 1j*0
    sources = np.random.uniform(0,1,(3,n))
    stmp = sources[:,0:ntest]

    nt = 1880
    targ = np.random.uniform(0,1,(3,nt))
    ttmp = targ[:,0:ntest]
    eps = 10**(-5)

    zk = 1.1 + 1j*0
    charges = np.random.uniform(0,1,n)
    dipvec = np.random.uniform(0,1,(3,n))

    outex=fmm.Output()

    itest = 0
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges,pg=1)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges,pgt=1)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot")

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,dipvec=dipvec,pg=1)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,dipvec=dipvec,pgt=1)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot")


    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges, \
        dipvec=dipvec,pg=1)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        dipvec=dipvec,pgt=1)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot")


    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges,pg=2)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=2)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot and grad")
        

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,dipvec=dipvec,pg=2)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,dipvec=dipvec, \
        pgt=2)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges, \
        dipvec=dipvec,pg=2)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges,
       dipvec=dipvec,pgt=2)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot and grad")
    

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges,pg=3)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=3)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    out2.hess = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=3)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot, grad, and hess")
        

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,dipvec=dipvec,pg=3)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,dipvec=dipvec, \
        pgt=3)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    out2.hess = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=3)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot, grad, and hess")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges, \
        dipvec=dipvec,pg=3)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges,
       dipvec=dipvec,pgt=3)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    out2.hess = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=3)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot and grad")
    

    itest=itest+1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=1)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,pgt=1)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot")

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
        dipvec=dipvec,pgt=1)
      
    out2=fmm.l3ddir(sources=sources,targets=ttmp,\
        dipvec=dipvec,pgt=1)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ, \
        charges=charges, \
        dipvec=dipvec,pgt=1)

    out2=fmm.l3ddir(sources=sources,targets=ttmp, \
        charges=charges, \
        dipvec=dipvec,pgt=1)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1)
    
    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=2)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,pgt=2)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot and grad")

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=2)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,\
    dipvec=dipvec,pgt=2)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipvec=dipvec,pgt=2)
    out2 =fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,\
        dipvec=dipvec,pgt=2)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot and grad")


    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=3)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,pgt=3)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=3)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot, grad, and hess")

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=3)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,\
    dipvec=dipvec,pgt=3)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=3)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot, grad, and hess")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipvec=dipvec,pgt=3)
    out2 =fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,\
        dipvec=dipvec,pgt=3)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=3)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot, grad, and hess")


    itest = itest+1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=1,pg=1)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges,pgt=1)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,pgt=1)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot")
    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
        dipvec=dipvec,pgt=1,pg=1)
    out2=fmm.l3ddir(sources=sources,targets=stmp, \
        dipvec=dipvec,pgt=1)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp, \
          dipvec=dipvec,pgt=1)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ, \
        charges=charges, \
        dipvec=dipvec,pgt=1,pg=1)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        dipvec=dipvec,pgt=1)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          dipvec=dipvec,pgt=1)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=2,pg=2)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=2)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          pgt=2)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot and grad")
    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=2,pg=2)
    out2=fmm.l3ddir(sources=sources,targets=stmp, \
        dipvec=dipvec,pgt=2)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,dipvec=dipvec, \
          pgt=2)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipvec=dipvec,pgt=2,pg=2)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        dipvec=dipvec,pgt=2)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          dipvec=dipvec,pgt=2)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot and grad")


    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=3,pg=3)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=3)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    outex.hess = out2.hesstarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          pgt=3)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    outex.hesstarg = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=3,pgt=3)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot, grad, hess")
    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=3,pg=3)
    out2=fmm.l3ddir(sources=sources,targets=stmp, \
        dipvec=dipvec,pgt=3)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    outex.hess = out2.hesstarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,dipvec=dipvec, \
          pgt=3)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    outex.hesstarg = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=3,pgt=3)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot, grad, hess")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipvec=dipvec,pgt=3,pg=3)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        dipvec=dipvec,pgt=3)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    outex.hess = out2.hesstarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          dipvec=dipvec,pgt=3)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    outex.hesstarg = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=3,pgt=3)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot and grad")



    nd = 2
    charges = np.random.uniform(0,1,(nd,n))
    dipvec = np.random.uniform(0,1,(nd,3,n))

    itest = itest+1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges,pg=1,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges,pgt=1,nd=nd)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot, vectorized")

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,dipvec=dipvec,pg=1,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,dipvec=dipvec,pgt=1,nd=nd)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,nd=nd)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot, vectorized")


    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges, \
        dipvec=dipvec,pg=1,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        dipvec=dipvec,pgt=1,nd=nd)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot, vectorized")


    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges,pg=2,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=2,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,nd=nd)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot and grad, vectorized")
        

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,dipvec=dipvec,pg=2,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,dipvec=dipvec, \
        pgt=2,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges, \
        dipvec=dipvec,pg=2,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges,
       dipvec=dipvec,pgt=2,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,nd=nd)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot and grad, vectorized")
    

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges,pg=3,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=3,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    out2.hess = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=3,nd=nd)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot, grad, hess vectorized")
        

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,dipvec=dipvec,pg=3,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,dipvec=dipvec, \
        pgt=3,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    out2.hess = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=3,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot, grad, hess vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,charges=charges, \
        dipvec=dipvec,pg=3,nd=nd)
    out2 = fmm.l3ddir(sources=sources,targets=stmp,charges=charges,
       dipvec=dipvec,pgt=3,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    out2.hess = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=3,nd=nd)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot, grad, hess vectorized")
    

    itest=itest+1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=1,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,pgt=1,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot, vectorized")

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
        dipvec=dipvec,pgt=1,nd=nd)
      
    out2=fmm.l3ddir(sources=sources,targets=ttmp,\
        dipvec=dipvec,pgt=1,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ, \
        charges=charges, \
        dipvec=dipvec,pgt=1,nd=nd)

    out2=fmm.l3ddir(sources=sources,targets=ttmp, \
        charges=charges, \
        dipvec=dipvec,pgt=1,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,nd=nd)
    
    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=2,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,pgt=2,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot and grad, vectorized")

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=2,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,\
    dipvec=dipvec,pgt=2,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipvec=dipvec,pgt=2,nd=nd)
    out2 =fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,\
        dipvec=dipvec,pgt=2,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot and grad, vectorized")


    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=3,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,pgt=3,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=3,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot, grad, hess, vectorized")

    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=3,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=ttmp,\
    dipvec=dipvec,pgt=3,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=3,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot, grad, hess vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipvec=dipvec,pgt=3,nd=nd)
    out2 =fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,\
        dipvec=dipvec,pgt=3,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=3,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot, grad, hess vectorized")

    itest = itest+1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=1,pg=1,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges,pgt=1,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges,pgt=1,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot, vectorized")
    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
        dipvec=dipvec,pgt=1,pg=1,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp, \
        dipvec=dipvec,pgt=1,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp, \
          dipvec=dipvec,pgt=1,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ, \
        charges=charges, \
        dipvec=dipvec,pgt=1,pg=1,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        dipvec=dipvec,pgt=1,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          dipvec=dipvec,pgt=1,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=2,pg=2,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=2,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          pgt=2,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot and grad, vectorized")
    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=2,pg=2,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp, \
        dipvec=dipvec,pgt=2,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,dipvec=dipvec, \
          pgt=2,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipvec=dipvec,pgt=2,pg=2,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        dipvec=dipvec,pgt=2,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          dipvec=dipvec,pgt=2,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot and grad, vectorized")


    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=3,pg=3,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=3,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    outex.hess = out2.hesstarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          pgt=3,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    outex.hesstarg = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=3,pgt=3,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot and grad, vectorized")
    itest = itest+1

    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=3,pg=3,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp, \
        dipvec=dipvec,pgt=3,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    outex.hess = out2.hesstarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,dipvec=dipvec, \
          pgt=3,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    outex.hesstarg = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=3,pgt=3,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.lfmm3d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipvec=dipvec,pgt=3,pg=3,nd=nd)
    out2=fmm.l3ddir(sources=sources,targets=stmp,charges=charges, \
        dipvec=dipvec,pgt=3,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    outex.hess = out2.hesstarg
    out2=fmm.l3ddir(sources=sources,targets=ttmp,charges=charges, \
          dipvec=dipvec,pgt=3,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    outex.hesstarg = out2.hesstarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=3,pgt=3,nd=nd)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot and grad, vectorized")

    if(sum(testres)==ntests):
        print("all lfmm tests succeeded")

if __name__ == "__main__":
    main()
