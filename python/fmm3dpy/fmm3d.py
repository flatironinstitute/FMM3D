import hfmm3d_fortran as hfmm
import lfmm3d_fortran as lfmm
import numpy as np


class Output():
    pot = None
    grad = None
    pottarg = None
    gradtarg = None

def hfmm3d(*,eps,zk,sources,charges=None,dipvec=None,
          targets=None,pg=0,pgt=0,nd=1):
    """
      This subroutine computes the N-body Helmholtz interactions
      in three dimensions where the interaction kernel is given by e^{ikr}/r 
      and its gradients. 

      u(x) = \sum_{j=1}^{N} c_{j} e^{ik |x-x_{j}|}/|x-x_{j}| + 
                   Grad( e^{ik |x-x_{j}|}/|x-x_{j}|) . v_{j} \, ,

      where c_{j} are the charge densities,  
      v_{j} are the dipole orientation vectors, and 
      x_{j} are the source locations.

      When x=x_{m}, the term corresponding to x_{m} is dropped from the
      sum


      Args:
        eps: float   
               precision requested
        zk: complex
               Helmholtz parameter - k
        sources: float(3,n)   
               source locations (x_{j})
        charges: complex(nd,n) or complex(n)
                charge densities (c_{j})
        dipole orientation vectors: complex(nd,3,n) or complex(3,n)
                dipole orientation vectors (v_{j})
        targets: float(3,nt)
                target locations (x)
        pg:  integer
               source eval flag
               potential at sources evaluated if pg = 1
               potenial and gradient at sources evaluated if pg=2

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
        
        nd:   integer
               number of densities

        Returns:
          out.pot  - potential at source locations if requested
          out.grad - gradient at source locations if requested
          out.pottarg  - potential at target locations if requested
          out.gradtarg - gradient at target locations if requested
              

    """
    out = Output()
    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    ns = sources.shape[1]
    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.size == ns, "Charges must be same length as second dimension of sources"
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None):
        if nd == 1:
            assert dipvec.shape[0] == 3 and dipvec.shape[1] == ns, "dipole vectors must be of shape [3,number of sources]"
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 3 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 3, "The first dimension of targets must be 3"
        iftarg = 1
    if(iftarg == 0 or pgt != 1 or pgt !=2):
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
                out.pot = hfmm.hfmm3dpartstosdp_vec(eps,zk,sources,dipvec,nd)
            if(nd == 1):
                out.pot = hfmm.hfmm3dpartstosdp(eps,zk,sources,dipvec)
                
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad = hfmm.hfmm3dpartstosdg_vec(eps,zk,sources,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad = hfmm.hfmm3dpartstosdg(eps,zk,sources,dipvec)
        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot = hfmm.hfmm3dpartstoscdp_vec(eps,zk,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot = hfmm.hfmm3dpartstoscdp(eps,zk,sources,charges,dipvec)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad = hfmm.hfmm3dpartstoscdg_vec(eps,zk,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad = hfmm.hfmm3dpartstoscdg(eps,zk,sources,charges,dipvec)
    
    if(pg !=1 and pg !=2 and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg = hfmm.hfmm3dpartstotcp_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg = hfmm.hfmm3dpartstotcp(eps,zk,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotcg_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotcg(eps,zk,sources,charges,targets)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg = hfmm.hfmm3dpartstotdp_vec(eps,zk,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg = hfmm.hfmm3dpartstotdp(eps,zk,sources,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotdg_vec(eps,zk,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotdg(eps,zk,sources,dipvec,targets)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg = hfmm.hfmm3dpartstotcdp_vec(eps,zk,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg = hfmm.hfmm3dpartstotcdp(eps,zk,sources,charges,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotcdg_vec(eps,zk,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = hfmm.hfmm3dpartstotcdg(eps,zk,sources,charges,dipvec,targets)
    

    if((pg == 1 or pg == 2) and targets is not None):
        assert pg == pgt, "if both potential or potential at gradient are requested at sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostcp_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostcp(eps,zk,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostcg_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostcg(eps,zk,sources,charges,targets)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostdp_vec(eps,zk,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostdp(eps,zk,sources,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostdg_vec(eps,zk,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostdg(eps,zk,sources,dipvec,targets)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostcdp_vec(eps,zk,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg = hfmm.hfmm3dpartstostcdp(eps,zk,sources,charges,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostcdg_vec(eps,zk,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = hfmm.hfmm3dpartstostcdg(eps,zk,sources,charges,dipvec,targets)

    return out



def lfmm3d(*,eps,sources,charges=None,dipvec=None,
          targets=None,pg=0,pgt=0,nd=1):
    """
      This subroutine computes the N-body Laplace interactions
      in three dimensions where the interaction kernel is given by 1/r 
      and its gradients. 

      u(x) = \sum_{j=1}^{N} c_{j} /|x-x_{j}| + 
                   Grad( 1/|x-x_{j}|) . v_{j} \, ,

      where c_{j} are the charge densities, 
      v_{j} are the dipole orientation vectors, and 
      x_{j} are the source locations.

      When x=x_{m}, the term corresponding to x_{m} is dropped from the
      sum


      Args:
        eps: float   
               precision requested
        sources: float(3,n)   
               source locations (x_{j})
        charges: float(nd,n) or float(n)
                charge densities (c_{j})
        dipole orientation vectors: float(nd,3,n) or float(3,n)
                dipole orientation vectors (v_{j})
        targets: float(3,nt)
                target locations (x)
        pg:  integer
               source eval flag
               potential at sources evaluated if pg = 1
               potenial and gradient at sources evaluated if pg=2

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
        
        nd:   integer
               number of densities

        Returns:
          out.pot  - potential at source locations if requested
          out.grad - gradient at source locations if requested
          out.pottarg  - potential at target locations if requested
          out.gradtarg - gradient at target locations if requested
              

    """

    out = Output()
    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    ns = sources.shape[1]
    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.size == ns, "Charges must be same length as second dimension of sources"
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None):
        if nd == 1:
            assert dipvec.shape[0] == 3 and dipvec.shape[1] == ns, "dipole vectors must be of shape [3,number of sources]"
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 3 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 3, "The first dimension of targets must be 3"
        iftarg = 1
    if(iftarg == 0 or pgt != 1 or pgt !=2):
        if(pg == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot = lfmm.rfmm3dpartstoscp_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot = lfmm.rfmm3dpartstoscp(eps,sources,charges)
        if(pg == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad = lfmm.rfmm3dpartstoscg_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad = lfmm.rfmm3dpartstoscg(eps,sources,charges)
        if(pg == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot = lfmm.rfmm3dpartstosdp_vec(eps,sources,dipvec,nd)
            if(nd == 1):
                out.pot = lfmm.rfmm3dpartstosdp(eps,sources,dipvec)
                
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad = lfmm.rfmm3dpartstosdg_vec(eps,sources,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad = lfmm.rfmm3dpartstosdg(eps,sources,dipvec)
        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot = lfmm.rfmm3dpartstoscdp_vec(eps,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot = lfmm.rfmm3dpartstoscdp(eps,sources,charges,dipvec)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad = lfmm.rfmm3dpartstoscdg_vec(eps,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad = lfmm.rfmm3dpartstoscdg(eps,sources,charges,dipvec)
    
    if(pg !=1 and pg !=2 and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg = lfmm.rfmm3dpartstotcp_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg = lfmm.rfmm3dpartstotcp(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg = lfmm.rfmm3dpartstotcg_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = lfmm.rfmm3dpartstotcg(eps,sources,charges,targets)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg = lfmm.rfmm3dpartstotdp_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg = lfmm.rfmm3dpartstotdp(eps,sources,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg = lfmm.rfmm3dpartstotdg_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = lfmm.rfmm3dpartstotdg(eps,sources,dipvec,targets)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg = lfmm.rfmm3dpartstotcdp_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg = lfmm.rfmm3dpartstotcdp(eps,sources,charges,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg = lfmm.rfmm3dpartstotcdg_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg = lfmm.rfmm3dpartstotcdg(eps,sources,charges,dipvec,targets)
    

    if((pg == 1 or pg == 2) and targets is not None):
        assert pg == pgt, "if both potential or potential at gradient are requested at sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg = lfmm.rfmm3dpartstostcp_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg = lfmm.rfmm3dpartstostcp(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = lfmm.rfmm3dpartstostcg_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = lfmm.rfmm3dpartstostcg(eps,sources,charges,targets)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg = lfmm.rfmm3dpartstostdp_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg = lfmm.rfmm3dpartstostdp(eps,sources,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = lfmm.rfmm3dpartstostdg_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = lfmm.rfmm3dpartstostdg(eps,sources,dipvec,targets)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg = lfmm.rfmm3dpartstostcdp_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg = lfmm.rfmm3dpartstostcdp(eps,sources,charges,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = lfmm.rfmm3dpartstostcdg_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg = lfmm.rfmm3dpartstostcdg(eps,sources,charges,dipvec,targets)

    return out

