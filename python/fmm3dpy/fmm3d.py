import hfmm3d_fortran as hfmm
#import lfmm3d_fortran as lfmm
import numpy as np


class Output():
    pot = None
    grad = None
    pottarg = None
    gradtarg = None

def hfmm3d(*,eps,zk,sources,charges=None,dipoles=None,dipvec=None,
          targets=None,pg=0,pgt=0,nd=1):
    out = Output()
    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    ns = sources.shape[1]
    ifcharge = 0
    ifdipole = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.size == ns, "Charges must be same length as second dimension of sources"
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipoles is not None or dipvec is not None):
        assert dipoles is not None, "Dipole vectors set but no dipole strength set"
        assert dipvec is not None, "Dipole strengths set but no dipole vectors set"
        if nd == 1:
            assert dipoles.size == ns, "Dipoles must be of same length as second dimension of sources"
            assert dipvec.shape[0] == 3 and dipvec.shape[1] == ns, "dipole vectors must be of shape [3,number of sources]"
        if nd>1:
            assert dipoles.shape[0] == nd and dipoles.shape[1] == ns, "Dipoles must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 3 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 3, "The first dimension of targets must be 3"
    if(targets == None or pgt != 1 or pgt !=2):
        if(pg == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot = hfmm.hfmm3dpartstoscp_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pot = hfmm.hfmm3dpartstoscp(eps,zk,sources,charges)
        if(pg == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad = hfmm.hfmm3dpartstoscg_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad = hfmm.hfmm3dpartstoscg(eps,zk,sources,charges)
        if(pg == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot = hfmm.hfmm3dpartstosdp_vec(eps,zk,sources,dipoles,dipvec,nd)
            if(nd == 1):
                out.pot = hfmm.hfmm3dpartstosdp(eps,zk,sources,dipoles,dipvec)
                
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad = hfmm.hfmm3dpartstosdg_vec(eps,zk,sources,dipoles,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad = hfmm.hfmm3dpartstosdg(eps,zk,sources,dipoles,dipvec)
        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot = hfmm.hfmm3dpartstoscdp_vec(eps,zk,sources,charges,dipoles,dipvec,nd)
            if(nd == 1):
                out.pot = hfmm.hfmm3dpartstoscdp(eps,zk,sources,charges,dipoles,dipvec)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad = hfmm.hfmm3dpartstoscdg_vec(eps,zk,sources,charges,dipoles,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad = hfmm.hfmm3dpartstoscdg(eps,zk,sources,charges,dipoles,dipvec)
    
    if((pg !=1 or pg !=2) and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg = hfmm.hfmm3dpartstotcp_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pottarg = hfmm.hfmm3dpartstotcp(eps,zk,sources,charges)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotcg_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotcg(eps,zk,sources,charges)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg = hfmm.hfmm3dpartstotdp_vec(eps,zk,sources,dipoles,dipvec,nd)
            if(nd == 1):
                out.pottarg = hfmm.hfmm3dpartstotdp(eps,zk,sources,dipoles,dipvec)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotdg_vec(eps,zk,sources,dipoles,dipvec,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotdg(eps,zk,sources,dipoles,dipvec)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg = hfmm.hfmm3dpartstotcdp_vec(eps,zk,sources,charges,dipoles,dipvec,nd)
            if(nd == 1):
                out.pottarg = hfmm.hfmm3dpartstotcdp(eps,zk,sources,charges,dipoles,dipvec)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotcdg_vec(eps,zk,sources,charges,dipoles,dipvec,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotcdg(eps,zk,sources,charges,dipoles,dipvec)
    

    if((pg == 1 or pg == 2) and targets is not None):
        assert pg == pgt, "if both potential or potential at gradient are requested at sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostcp_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostcp(eps,zk,sources,charges)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostcg_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostcg(eps,zk,sources,charges)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostdp_vec(eps,zk,sources,dipoles,dipvec,nd)
            if(nd == 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostdp(eps,zk,sources,dipoles,dipvec)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostdg_vec(eps,zk,sources,dipoles,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostdg(eps,zk,sources,dipoles,dipvec)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostcdp_vec(eps,zk,sources,charges,dipoles,dipvec,nd)
            if(nd == 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostcdp(eps,zk,sources,charges,dipoles,dipvec)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostcdg_vec(eps,zk,sources,charges,dipoles,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostcdg(eps,zk,sources,charges,dipoles,dipvec)

    return out

