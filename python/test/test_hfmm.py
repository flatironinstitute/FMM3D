#!/usr/bin/env python

import fmm3dpy as fmm
import numpy as np

def main():
    test_hfmm()

def test_hfmm():
    ntest = 36
    testres = np.zeros(ntest)
    #
    #  This is a testing code for making sure all the 
    #  fmm routines are accessible through fmm3d.py
    #

    n = 2000
    zk = 1.1 + 1j*0
    sources = np.random.uniform(0,1,(3,n))

    nt = 1880
    targ = np.random.uniform(0,1,(3,nt))
    eps = 10**(-5)

    zk = 1.1 + 1j*0
    charges = np.random.uniform(0,1,n)+ 1j*np.random.uniform(0,1,n)
    dipvec = np.random.uniform(0,1,(3,n))+ 1j*np.random.uniform(0,1,(3,n))

    itest = 0
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,charges=charges,pg=1)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,dipvec=dipvec,pg=1)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,charges=charges, \
        dipvec=dipvec,pg=1)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,charges=charges,pg=2)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==(3,n) and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot and grad")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,dipvec=dipvec,pg=2)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==(3,n) and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,charges=charges, \
        dipvec=dipvec,pg=2)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==(3,n) and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot and grad")


    itest=itest+1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,pgt=1)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,\
        dipvec=dipvec,pgt=1)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ, \
        charges=charges, \
        dipvec=dipvec,pgt=1)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,pgt=2)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==(3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot and grad")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=2)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==(3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,\
         \
        dipvec=dipvec,pgt=2)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==(3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot and grad")

    itest = itest+1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,pgt=1,pg=1)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,\
        dipvec=dipvec,pgt=1,pg=1)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ, \
        charges=charges, \
        dipvec=dipvec,pgt=1,pg=1)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,pgt=2,pg=2)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==(3,n) and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==(3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot and grad")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=2,pg=2)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==(3,n) and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==(3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,\
         \
        dipvec=dipvec,pgt=2,pg=2)
    a = (np.shape(out.pot) == (n,) and np.shape(out.grad)==(3,n) and \
    np.shape(out.pottarg) == (nt,) and np.shape(out.gradtarg)==(3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot and grad")

    nd = 2
    charges = np.random.uniform(0,1,(nd,n))+ 1j*np.random.uniform(0,1,(nd,n))
    dipvec = np.random.uniform(0,1,(nd,3,n))+ 1j*np.random.uniform(0,1,(nd,3,n))

    itest = itest+1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,charges=charges,pg=1,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot, vectorized")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,dipvec=dipvec,pg=1,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,charges=charges, \
        dipvec=dipvec,pg=1,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,charges=charges,pg=2,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==(nd,3,n) and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot and grad, vectorized")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,dipvec=dipvec,pg=2,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==(nd,3,n) and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,charges=charges, \
        dipvec=dipvec,pg=2,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==(nd,3,n) and \
    np.shape(out.pottarg) == () and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot and grad, vectorized")


    itest=itest+1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,pgt=1,nd=nd)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot, vectorized")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,\
        dipvec=dipvec,pgt=1,nd=nd)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ, \
        charges=charges, \
        dipvec=dipvec,pgt=1,nd=nd)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,pgt=2,nd=nd)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==(nd,3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot and grad, vectorized")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=2,nd=nd)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==(nd,3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,\
         \
        dipvec=dipvec,pgt=2,nd=nd)
    a = (np.shape(out.pot) == () and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==(nd,3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot and grad, vectorized")

    itest = itest+1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,pgt=1,pg=1,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot, vectorized")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,\
        dipvec=dipvec,pgt=1,pg=1,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ, \
        charges=charges, \
        dipvec=dipvec,pgt=1,pg=1,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==() and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==())

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,pgt=2,pg=2,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==(nd,3,n) and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==(nd,3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot and grad, vectorized")
    itest = itest+1

    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,\
    dipvec=dipvec,pgt=2,pg=2,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==(nd,3,n) and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==(nd,3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.hfmm3d(eps=eps,zk=zk,sources=sources,targets=targ,charges=charges,\
         \
        dipvec=dipvec,pgt=2,pg=2,nd=nd)
    a = (np.shape(out.pot) == (nd,n) and np.shape(out.grad)==(nd,3,n) and \
    np.shape(out.pottarg) == (nd,nt) and np.shape(out.gradtarg)==(nd,3,nt))

    if(a):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot and grad, vectorized")



    if(sum(testres)==ntest):
        print("all tests succeeded")

if __name__ == "__main__":
    main()
