from . import hfmm3d_fortran as hfmm
from . import lfmm3d_fortran as lfmm
from . import emfmm3d_fortran as emfmm
from . import stfmm3d_fortran as stfmm
import numpy as np
import numpy.linalg as la


class Output():
    pot = None
    grad = None
    hess = None
    pottarg = None
    gradtarg = None
    hesstarg = None
    E = None
    curlE = None
    divE = None
    Etarg = None
    curlEtarg = None
    divEtarg = None
    pre = None
    pretarg = None
    ier = 0

def hfmm3d(*,eps,zk,sources,charges=None,dipvec=None,
          targets=None,pg=0,pgt=0,nd=1):
    r"""
      This subroutine computes the N-body Helmholtz interactions
      in three dimensions where the interaction kernel is given by e^{ikr}/r 
      and its gradients. 

      .. math::

          u(x) = \sum_{j=1}^{N} c_{j} \\frac{e^{ik \|x-x_{j}\|}}{\|x-x_{j}\|} - v_{j} \cdot \\nabla \left( \\frac{e^{ik \|x-x_{j}\|}}{\|x-x_{j}\|} \\right)  \, ,

      where $c_{j}$ are the charge densities,  
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the
      sum

      Args:
        eps (float): precision requested
        zk (complex): Helmholtz parameter 
        sources (float(3,n)): source locations ($x_{j}$)
        charges (complex(nd,n) or complex(n)): charge densities ($c_{j}$)
        dipvec (complex(nd,3,n) or complex(3,n)): dipole orientation vectors ($v_{j}$)
        targets (float(3,nt)): target locations (x)
        pg (integer): source eval flag. Potential at sources evaluated if pg = 1. Potenial and gradient at sources evaluated if pg=2
        pgt (integer): target eval flag. Potential at targets evaluated if pgt = 1. Potenial and gradient at targets evaluated if pgt=2
        nd (integer): number of densities

      Returns:
        Returns an object of type Output (out) with the following variables

        out.pot: potential at source locations if requested
        out.grad: gradient at source locations if requested
        out.pottarg: potential at target locations if requested
        out.gradtarg: gradient at target locations if requested

      Example:
        see hmmexample.py
    r"""
    
    out = Output()
    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1

    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None):
        if nd == 1 and ns>1:
            assert dipvec.shape[0] == 3 and dipvec.shape[1] == ns, "dipole vectors must be of shape [3,number of sources]"
        if nd == 1 and ns==1:
            assert dipvec.shape[0] == 3, "dipole vectors must be of shape [3,number of sources]"
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 3 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 3, "The first dimension of targets must be 3"
        iftarg = 1
    if(iftarg == 0 or pgt != 1 or pgt !=2):
        if(pg == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.ier = hfmm.hfmm3d_s_c_p_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pot,out.ier = hfmm.hfmm3d_s_c_p(eps,zk,sources,charges)
        if(pg == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.ier = hfmm.hfmm3d_s_c_g_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = hfmm.hfmm3d_s_c_g(eps,zk,sources,charges)
        if(pg == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = hfmm.hfmm3d_s_d_p_vec(eps,zk,sources,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = hfmm.hfmm3d_s_d_p(eps,zk,sources,dipvec)
                
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = hfmm.hfmm3d_s_d_g_vec(eps,zk,sources,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = hfmm.hfmm3d_s_d_g(eps,zk,sources,dipvec)
        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = hfmm.hfmm3d_s_cd_p_vec(eps,zk,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = hfmm.hfmm3d_s_cd_p(eps,zk,sources,charges,dipvec)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = hfmm.hfmm3d_s_cd_g_vec(eps,zk,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = hfmm.hfmm3d_s_cd_g(eps,zk,sources,charges,dipvec)
    
    if(pg !=1 and pg !=2 and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.ier = hfmm.hfmm3d_t_c_p_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = hfmm.hfmm3d_t_c_p(eps,zk,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_t_c_g_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_t_c_g(eps,zk,sources,charges,targets)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = hfmm.hfmm3d_t_d_p_vec(eps,zk,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = hfmm.hfmm3d_t_d_p(eps,zk,sources,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_t_d_g_vec(eps,zk,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_t_d_g(eps,zk,sources,dipvec,targets)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = hfmm.hfmm3d_t_cd_p_vec(eps,zk,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = hfmm.hfmm3d_t_cd_p(eps,zk,sources,charges,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_t_cd_g_vec(eps,zk,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_t_cd_g(eps,zk,sources,charges,dipvec,targets)
    

    if((pg == 1 or pg == 2) and targets is not None):
        assert pg == pgt, "if both potential or potential at gradient are requested at sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm3d_st_c_p_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm3d_st_c_p(eps,zk,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_st_c_g_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_st_c_g(eps,zk,sources,charges,targets)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm3d_st_d_p_vec(eps,zk,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm3d_st_d_p(eps,zk,sources,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_st_d_g_vec(eps,zk,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_st_d_g(eps,zk,sources,dipvec,targets)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm3d_st_cd_p_vec(eps,zk,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm3d_st_cd_p(eps,zk,sources,charges,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_st_cd_g_vec(eps,zk,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm3d_st_cd_g(eps,zk,sources,charges,dipvec,targets)

    return out



def lfmm3d(*,eps,sources,charges=None,dipvec=None,
          targets=None,pg=0,pgt=0,nd=1):
    r"""
      This subroutine computes the N-body Laplace interactions
      in three dimensions where the interaction kernel is given by 1/r 
      and its gradients. 


      .. math:: 

          u(x) = \sum_{j=1}^{N} c_{j} / \|x-x_{j}\| + v_{j} \cdot \\nabla( 1/\|x-x_{j}\|)  \, ,

      where $c_{j}$ are the charge densities, 
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        eps: float
             precision requested

        sources: float(3,n)   
                 source locations (x_{j})
        charges: float(nd,n) or float(n)
                 charge densities (c_{j})
        dipole: float(nd,3,n) or float(3,n)
                                    dipole orientation vectors (v_{j})
        targets: float(3,nt)
                target locations (x)
        pg:  integer
               source eval flag
               potential at sources evaluated if pg = 1
               potenial and gradient at sources evaluated if pg=2
               potential, gradient and hessian at sources evaluated if pg=3

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potential, gradient and hessian at targets evaluated if pgt=3
        
        nd:   integer
               number of densities

      Returns:
        out.pot: potential at source locations if requested
        out.grad: gradient at source locations if requested
        out.hess: hessian at source locations if requested
        out.pottarg: potential at target locations if requested
        out.gradtarg: gradient at target locations if requested
        out.hesstarg: hessian at target locations if requested
      
      Example:
        see lfmmexample.py
    r"""

    out = Output()
    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None):
        if nd == 1 and ns>1:
            assert dipvec.shape[0] == 3 and dipvec.shape[1] == ns, "dipole vectors must be of shape [3,number of sources]"
        if nd == 1 and ns==1:
            assert dipvec.shape[0] == 3, "dipole vectors must be of shape [3,number of sources]"
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 3 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 3, "The first dimension of targets must be 3"
        iftarg = 1
#
# sources -> sources routines
#
    if(iftarg == 0 or pgt != 1 or pgt !=2 or pgt !=3):
        if(pg == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.ier = lfmm.lfmm3d_s_c_p_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.lfmm3d_s_c_p(eps,sources,charges)
        if(pg == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.lfmm3d_s_c_g_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.lfmm3d_s_c_g(eps,sources,charges)
        if(pg == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm3d_s_c_h_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm3d_s_c_h(eps,sources,charges)


        if(pg == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = lfmm.lfmm3d_s_d_p_vec(eps,sources,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.lfmm3d_s_d_p(eps,sources,dipvec)
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.lfmm3d_s_d_g_vec(eps,sources,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.lfmm3d_s_d_g(eps,sources,dipvec)
        if(pg == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm3d_s_d_h_vec(eps,sources,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm3d_s_d_h(eps,sources,dipvec)


        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = lfmm.lfmm3d_s_cd_p_vec(eps,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.lfmm3d_s_cd_p(eps,sources,charges,dipvec)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.lfmm3d_s_cd_g_vec(eps,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.lfmm3d_s_cd_g(eps,sources,charges,dipvec)
        if(pg == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm3d_s_cd_h_vec(eps,sources,charges,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm3d_s_cd_h(eps,sources,charges,dipvec)

#
# sources -> targets routines
#


    if(pg !=1 and pg !=2 and pg !=3 and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.lfmm3d_t_c_p_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.lfmm3d_t_c_p(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_t_c_g_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_t_c_g(eps,sources,charges,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_t_c_h_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_t_c_h(eps,sources,charges,targets)


        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.lfmm3d_t_d_p_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.lfmm3d_t_d_p(eps,sources,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_t_d_g_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_t_d_g(eps,sources,dipvec,targets)
        if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_t_d_h_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_t_d_h(eps,sources,dipvec,targets)


        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.lfmm3d_t_cd_p_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.lfmm3d_t_cd_p(eps,sources,charges,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_t_cd_g_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_t_cd_g(eps,sources,charges,dipvec,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_t_cd_h_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_t_cd_h(eps,sources,charges,dipvec,targets)
    
#
# sources to sources + targets
#
    if((pg == 1 or pg == 2 or pg == 3) and targets is not None):
        assert pg == pgt, "if output is requested at both sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm3d_st_c_p_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm3d_st_c_p(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_st_c_g_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_st_c_g(eps,sources,charges,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_st_c_h_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_st_c_h(eps,sources,charges,targets)


        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm3d_st_d_p_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm3d_st_d_p(eps,sources,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_st_d_g_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_st_d_g(eps,sources,dipvec,targets)
        if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_st_d_h_vec(eps,sources,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_st_d_h(eps,sources,dipvec,targets)


        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm3d_st_cd_p_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm3d_st_cd_p(eps,sources,charges,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_st_cd_g_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm3d_st_cd_g(eps,sources,charges,dipvec,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_st_cd_h_vec(eps,sources,charges,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm3d_st_cd_h(eps,sources,charges,dipvec,targets)

    return out

def emfmm3d(*,eps,zk,sources,h_current=None,e_current=None,e_charge=None,targets=None,ifE=0,ifcurlE=0,ifdivE=0,nd=1):
    r"""
      This function subrourine computes
          E = curl S_{k}[h_current] + S_{k}[e_current] + grad S_{k}[e_charge]  -- (1)
      using the vector Helmholtz fmm.
      The subroutine also computes divE, curlE
      with appropriate flags

      Remark: the subroutine uses a stabilized representation
      for computing the divergence by using integration by parts
      wherever possible. If the divergence is not requested, then the
      helmholtz fmm is called with 3*nd densities, while if the divergence
      is requested, then the helmholtz fmm is calld with 4*nd densities

      Args:
        eps (float): precision requested
        zk (complex): Helmholtz parameter
        sources (float(3,n)): source locations
        h_current (complex(nd,3,n) or complex(3,n)): a vector source
        e_current (complex(nd,3,n) or complex(3,n)): b vector source
        e_charge (complex(nd,n) or complex(n)): e_charge source
        targets (float(3,nt)): target locations
        ifE (integer): E is returned at the target locations if ifE = 1
        ifcurlE (integer): curl E is returned at the target locations if ifcurlE = 1
        ifdivE (integer): div E is returned at the target locations if ifdivE = 1
        nd (integer): number of densities

      Returns:
        Returns an object of type Output (out) with the following variables

        out.E: E field defined in (1) above at target locations if requested
        out.curlE: curl of E field at target locations if requested
        out.divE: divergence of E at target locations if requested

      Example:
        see emfmmexample.py
    r"""
    out = Output()

    if(targets is None):
        print("Nothing to compute, set targets")
        return out
    if(ifE == 0 and ifcurlE == 0 and ifdivE == 0):
        print("Nothing to compute, set either ifE, ifcurlE or ifdivE to non-zero")
        return out

    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    assert targets.shape[0] == 3, "The first dimension of targets must be 3"
    if(np.size(np.shape(targets))==2):
        nt = targets.shape[1]
    if(np.size(np.shape(targets))==1):
        nt = 1

    ifh_current = 0
    ife_current = 0
    ife_charge  = 0

    if(h_current is not None):
        if(nd == 1 and ns>1):
            assert h_current.shape[0] == 3 and h_current.shape[1] == ns, "h_current vectors must be of shape [3,number of sources]"
        if(nd == 1 and ns==1):
            assert h_current.shape[0] == 3, "h_current vectors must be of shape [3,number of sources]"
        if(nd>1):
            assert h_current.shape[0] == nd and h_current.shape[1] == 3 and h_current.shape[2] == ns, "h_current vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        h_current = h_current.reshape([nd,3,ns])
        ifh_current = 1
    else:
        h_current = np.zeros([nd,3,ns],dtype=complex)

    if(e_current is not None):
        if(nd == 1 and ns>1):
            assert e_current.shape[0] == 3 and e_current.shape[1] == ns, "e_current vectors must be of shape [3,number of sources]"
        if(nd == 1 and ns==1):
            assert e_current.shape[0] == 3, "e_current vectors must be of shape [3,number of sources]"
        if(nd>1):
            assert e_current.shape[0] == nd and e_current.shape[1] == 3 and e_current.shape[2] == ns, "e_current vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        e_current = e_current.reshape([nd,3,ns])
        ife_current = 1
    else:
        e_current = np.zeros([nd,3,ns],dtype=complex)

    if(e_charge is not None):
        if(nd == 1):
            assert e_charge.shape[0] == ns, "e_charge must be same length as second dimension of sources"
        if(nd>1):
            assert e_charge.shape[0] == nd and e_charge.shape[1]==ns, "e_charge must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        e_charge = e_charge.reshape([nd,ns])
        ife_charge = 1
    else:
        e_charge = np.zeros([nd,ns],dtype=complex)

    out.E,out.curlE,out.divE,out.ier = emfmm.emfmm3d(eps,zk,sources,ifh_current,h_current,ife_current,e_current,ife_charge,e_charge,targets,ifE,ifcurlE,ifdivE,nd,ns,nt)

    if(ifE==0):
        out.E = None
    if(ifcurlE==0):
        out.curlE = None
    if(ifdivE==0):
        out.divE = None

    return out

def stfmm3d(*,eps,sources,stoklet=None,strslet=None,strsvec=None,targets=None,ifppreg=0,ifppregtarg=0,nd=1):
    r"""
      Stokes FMM in R^{3}: evaluate all pairwise particle
      interactions (ignoring self-interactions) and
      interactions with targs.
 
      This routine computes sums of the form
 
        u(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
                 + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
 
      where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
      stresslet charge, and nu^{(m)} is the stresslet orientation
      (note that each of these is a 3 vector per source point y^{(m)}).
      For x a source point, the self-interaction in the sum is omitted.
 
      Optionally, the associated pressure p(x) and gradient grad u(x)
      are returned
 
        p(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
           + sum_m T_{ijk}(x,y^{(m)}) PI_{jk} mu^{(m)}_j nu^{(m)}_k
 
        grad u(x) = grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
                 + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]

      Args:
        eps: float   
               precision requested
        sources: float(3,n)   
               source locations
        stoklet: float(nd,3,n) or float(3,n)
               Stokeslet charge strengths (sigma vectors above)
        strslet: float(nd,3,n) or float(3,n)
               stresslet strengths (mu vectors above)
        strsvec: float(nd,3,n) or float(3,n)
               stresslet orientations (nu vectors above)
        targets: float(3,nt)
               target locations (x)

        ifppreg: integer
               flag for evaluating potential, gradient, and pressure
               at the sources
               ifppreg = 1, only potential
               ifppreg = 2, potential and pressure
               ifppreg = 3, potential, pressure, and gradient

        ifppregtarg: integer
               flag for evaluating potential, gradient, and pressure
               at the targets
               ifppregtarg = 1, only potential
               ifppregtarg = 2, potential and pressure
               ifppregtarg = 3, potential, pressure, and gradient

        nd:   integer
               number of densities

      Returns:
        out.pot: velocity at source locations if requested
        out.pre: pressure at source locations if requested
        out.grad: gradient of velocity at source locations if requested
        out.pottarg: velocity at target locations if requested
        out.pretarg: pressure at target locations if requested
        out.gradtarg: gradient of velocity at target locations if requested
              
      Example:
        see stfmmexample.py

    r"""
    out = Output()

    if(ifppreg == 0 and ifppregtarg == 0):
        print("Nothing to compute, set either ifppreg or ifppregtarg to non-zero")
        return out

    if(stoklet is None and strslet is None and strsvec is None):
        print("Nothing to compute, set either stoklet or strslet+strsvec to non-None")
        return out

    if(strslet is not None and strsvec is None):
        print("strslet and strsvec mush be both None or both not None")
        return out
    if(strslet is None and strsvec is not None):
        print("strslet and strsvec mush be both None or both not None")
        return out

    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    if(targets is not None):
        assert targets.shape[0] == 3, "The first dimension of targets must be 3"
        if(np.size(np.shape(targets))==2):
            nt = targets.shape[1]
        if(np.size(np.shape(targets))==1):
            nt = 1
    else:
        targets = np.zeros([3,0],dtype='double')
        nt = 0

    ifstoklet = 0
    ifstrslet = 0

    if(stoklet is not None):
        if(nd == 1 and ns>1):
            assert stoklet.shape[0] == 3 and stoklet.shape[1] == ns, "stoklet vectors must be of shape [3,number of sources]"
        if(nd == 1 and ns==1):
            assert stoklet.shape[0] == 3, "stoklet vectors must be of shape [3,number of sources]"
        if(nd>1):
            assert stoklet.shape[0] == nd and stoklet.shape[1] == 3 and stoklet.shape[2] == ns, "stoklet vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        stoklet = stoklet.reshape([nd,3,ns])
        ifstoklet = 1
    else:
        stoklet = np.zeros([nd,3,ns],dtype='double')

    if(strslet is not None and strsvec is not None):
        if(nd == 1 and ns>1):
            assert strslet.shape[0] == 3 and strslet.shape[1] == ns, "strslet vectors must be of shape [3,number of sources]"
            assert strsvec.shape[0] == 3 and strsvec.shape[1] == ns, "strsvec vectors must be of shape [3,number of sources]"
        if(nd == 1 and ns==1):
            assert strslet.shape[0] == 3, "strslet vectors must be of shape [3,number of sources]"
            assert strsvec.shape[0] == 3, "strsvec vectors must be of shape [3,number of sources]"
        if(nd>1):
            assert strslet.shape[0] == nd and strslet.shape[1] == 3 and strslet.shape[2] == ns, "strslet vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
            assert strsvec.shape[0] == nd and strsvec.shape[1] == 3 and strsvec.shape[2] == ns, "strsvec vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        strslet = strslet.reshape([nd,3,ns])
        strsvec = strsvec.reshape([nd,3,ns])
        ifstrslet = 1
    else:
        strslet = np.zeros([nd,3,ns],dtype='double')
        strsvec = np.zeros([nd,3,ns],dtype='double')

    out.pot,out.pre,out.grad,out.pottarg,out.pretarg,out.gradtarg,out.ier = stfmm.stfmm3d(eps,sources,ifstoklet,stoklet,ifstrslet,strslet,strsvec,ifppreg,targets,ifppregtarg,nd,ns,nt)

    if(ifppreg < 3):
        out.grad = None
    if(ifppregtarg < 3):
        out.gradtarg = None

    if(ifppreg < 2):
        out.pre = None
    if(ifppregtarg < 2):
        out.pretarg = None

    if(ifppreg < 1):
        out.pot = None
    if(ifppregtarg < 1):
        out.pottarg = None

    return out

def h3ddir(*,zk,sources,targets,charges=None,dipvec=None,
          pgt=0,nd=1,thresh=1e-16):
    r"""
      This subroutine computes the N-body Helmholtz interactions
      in three dimensions where the interaction kernel is given by $e^{ikr}/r$ 
      and its gradients. 


      .. math::

          u(x) = \sum_{j=1}^{N} c_{j} e^{ik |x-x_{j}|}/|x-x_{j}| - \\nabla( e^{ik |x-x_{j}|}/|x-x_{j}|) \cdot v_{j} \, ,

      where $c_{j}$ are the charge densities,  
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When |x-x_{m}| \leq thresh, the term corresponding to $x_{m}$ is dropped from the
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

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
        
        nd:   integer
               number of densities
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        out.pottarg  - potential at target locations if requested
        out.gradtarg - gradient at target locations if requested
              
      Example:
        see hfmmexample.py
    r"""

    out = Output()
    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    if(pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None):
        if nd == 1 and ns>1:
            assert dipvec.shape[0] == 3 and dipvec.shape[1] == ns, "dipole vectors must be of shape [3,number of sources]"
            dipvec=dipvec.reshape(1,3,ns)
        if nd == 1 and ns==1:
            assert dipvec.shape[0] == 3, "dipole vectors must be of shape [3,number of sources]"
            dipvec=dipvec.reshape(1,3,ns)
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 3 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1

    assert targets.shape[0] == 3, "The first dimension of targets must be 3"
    nt = targets.shape[1]
    if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
        out.pottarg = hfmm.h3ddirectcp(zk,sources,charges,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg = hfmm.h3ddirectcg(zk,sources,charges,targets,thresh)
    if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
        out.pottarg = hfmm.h3ddirectdp(zk,sources,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg = hfmm.h3ddirectdg(zk,sources,dipvec,targets,thresh)
    if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
        out.pottarg = hfmm.h3ddirectcdp(zk,sources,charges,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg = hfmm.h3ddirectcdg(zk,sources,charges,dipvec,targets,thresh)

    if(nd == 1):
        if(ifcharge==1):
            charges = charges.reshape(ns,)
        if(ifdipole==1):
            dipvec = dipvec.reshape(3,ns)
        if(pgt>0):
            out.pottarg = out.pottarg.reshape(nt,)
        if(pgt==2):
            out.gradtarg = out.gradtarg.reshape(3,nt)


    return out



def l3ddir(*,sources,targets,charges=None,dipvec=None,
          pgt=0,nd=1,thresh=1e-16):
    r"""
      This subroutine computes the N-body Laplace interactions
      in three dimensions where the interaction kernel is given by $1/r$ 
      and its gradients. 


      .. math::

          u(x) = \sum_{j=1}^{N} c_{j} /|x-x_{j}| -  \\nabla( 1/|x-x_{j}|) \cdot v_{j} \, ,

      where $c_{j}$ are the charge densities, 
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When |x-x_{m}|leq thresh, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        sources: float(3,n)   
               source locations (x_{j})
        charges: float(nd,n) or float(n)
                charge densities (c_{j})
        dipole orientation vectors: float(nd,3,n) or float(3,n)
                dipole orientation vectors (v_{j})
        targets: float(3,nt)
                target locations (x)

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potenial, gradient, and hessians at targets evaluated if pgt=3
        
        nd:   integer
               number of densities
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        out.pottarg  - potential at target locations if requested
        out.gradtarg - gradient at target locations if requested
        out.hesstarg - hessian at target locations if requested
              
      Example:
        see lfmmexample.py

    r"""

    out = Output()
    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    if(pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None):
        if nd == 1 and ns>1:
            assert dipvec.shape[0] == 3 and dipvec.shape[1] == ns, "dipole vectors must be of shape [3,number of sources]"
            dipvec=dipvec.reshape(1,3,ns)
        if nd == 1 and ns==1:
            assert dipvec.shape[0] == 3, "dipole vectors must be of shape [3,number of sources]"
            dipvec=dipvec.reshape(1,3,ns)
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 3 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1

    assert targets.shape[0] == 3, "The first dimension of targets must be 3"
    nt = targets.shape[1]
    if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
        out.pottarg = lfmm.l3ddirectcp(sources,charges,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg = lfmm.l3ddirectcg(sources,charges,targets,thresh)
    if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.l3ddirectch(sources,charges,targets,thresh)
    if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
        out.pottarg = lfmm.l3ddirectdp(sources,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg = lfmm.l3ddirectdg(sources,dipvec,targets,thresh)
    if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.l3ddirectdh(sources,dipvec,targets,thresh)
    if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
        out.pottarg = lfmm.l3ddirectcdp(sources,charges,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg = lfmm.l3ddirectcdg(sources,charges,dipvec,targets,thresh)
    if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.l3ddirectcdh(sources,charges,dipvec,targets,thresh)

    if(nd == 1):
        if(ifcharge == 1):
            charges = charges.reshape(ns,)
        if(ifdipole ==1): 
            dipvec = dipvec.reshape(3,ns)
        if(pgt>0):
            out.pottarg = out.pottarg.reshape(nt,)
        if(pgt==2):
            out.gradtarg = out.gradtarg.reshape(3,nt)
        if(pgt==3):
            out.hesstarg = out.hesstarg.reshape(6,nt)

    return out

def em3ddir(*,eps,zk,sources,h_current=None,e_current=None,e_charge=None,targets=None,ifE=0,ifcurlE=0,ifdivE=0,nd=1,thresh=1e-16):
    r"""
      This function subrourine computes
          E = curl S_{k}[h_current] + S_{k}[e_current] + grad S_{k}[e_charge]  -- (1)
      using the vector Helmholtz fmm.
      The subroutine also computes divE, curlE
      with appropriate flags

      Remark: the subroutine uses a stabilized representation
      for computing the divergence by using integration by parts
      wherever possible. If the divergence is not requested, then the
      helmholtz fmm is called with 3*nd densities, while if the divergence
      is requested, then the helmholtz fmm is calld with 4*nd densities

      Args:
        eps (float): precision requested
        zk (complex): Helmholtz parameter
        sources (float(3,n)): source locations
        h_current (complex(nd,3,n) or complex(3,n)): a vector source
        e_current (complex(nd,3,n) or complex(3,n)): b vector source
        e_charge (complex(nd,n) or complex(n)): e_charge source
        targets (float(3,nt)): target locations
        ifE (integer): E is returned at the target locations if ifE = 1
        ifcurlE (integer): curl E is returned at the target locations if ifcurlE = 1
        ifdivE (integer): div E is returned at the target locations if ifdivE = 1
        nd (integer): number of densities
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        Returns an object of type Output (out) with the following variables

        out.E: E field defined in (1) above at target locations if requested
        out.curlE: curl of E field at target locations if requested
        out.divE: divergence of E at target locations if requested

      Example:
        see emfmmexample.py
    r"""
    out = Output()

    if(targets == None):
        print("Nothing to compute, set targets")
        return out
    if(ifE == 0 and ifcurlE == 0 and ifdivE == 0):
        print("Nothing to compute, set either ifE, ifcurlE or ifdivE to non-zero")
        return out

    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    assert targets.shape[0] == 3, "The first dimension of targets must be 3"
    if(np.size(np.shape(targets))==2):
        nt = targets.shape[1]
    if(np.size(np.shape(targets))==1):
        nt = 1

    ifh_current = 0
    ife_current = 0
    ife_charge  = 0

    if(h_current is not None):
        if(nd == 1 and ns>1):
            assert h_current.shape[0] == 3 and h_current.shape[1] == ns, "h_current vectors must be of shape [3,number of sources]"
        if(nd == 1 and ns==1):
            assert h_current.shape[0] == 3, "h_current vectors must be of shape [3,number of sources]"
        if(nd>1):
            assert h_current.shape[0] == nd and h_current.shape[1] == 3 and h_current.shape[2] == ns, "h_current vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        h_current = h_current.reshape([nd,3,ns])
        ifh_current = 1
    else:
        h_current = np.zeros([nd,3,ns],dtype=complex)

    if(e_current is not None):
        if(nd == 1 and ns>1):
            assert e_current.shape[0] == 3 and e_current.shape[1] == ns, "e_current vectors must be of shape [3,number of sources]"
        if(nd == 1 and ns==1):
            assert e_current.shape[0] == 3, "e_current vectors must be of shape [3,number of sources]"
        if(nd>1):
            assert e_current.shape[0] == nd and e_current.shape[1] == 3 and e_current.shape[2] == ns, "e_current vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        e_current = e_current.reshape([nd,3,ns])
        ife_current = 1
    else:
        e_current = np.zeros([nd,3,ns],dtype=complex)

    if(e_charge is not None):
        if(nd == 1):
            assert e_charge.shape[0] == ns, "e_charge must be same length as second dimension of sources"
        if(nd>1):
            assert e_charge.shape[0] == nd and e_charge.shape[1]==ns, "e_charge must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        e_charge = e_charge.reshape([nd,ns])
        ife_charge = 1
    else:
        e_charge = np.zeros([nd,ns],dtype=complex)

    out.E,out.curlE,out.divE = emfmm.em3ddirect(eps,zk,sources,ifh_current,h_current,ife_current,e_current,ife_charge,e_charge,targets,ifE,ifcurlE,ifdivE,thresh,nd,ns,nt)

    if(ifE==0):
        out.E = None
    if(ifcurlE==0):
        out.curlE = None
    if(ifdivE==0):
        out.divE = None

    return out

def st3ddir(*,eps,sources,stoklet=None,strslet=None,strsvec=None,targets=None,ifppreg=0,ifppregtarg=0,nd=1,thresh=1e-16):
    r"""
      This subroutine evaluates all pairwise particle
      interactions (ignoring self-interactions) and
      interactions with targs.
 
      This routine computes sums of the form
 
        u(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
                 + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
 
      where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
      stresslet charge, and nu^{(m)} is the stresslet orientation
      (note that each of these is a 3 vector per source point y^{(m)}).
      For x a source point, the self-interaction in the sum is omitted.
 
      Optionally, the associated pressure p(x) and gradient grad u(x)
      are returned
 
        p(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
           + sum_m T_{ijk}(x,y^{(m)}) PI_{jk} mu^{(m)}_j nu^{(m)}_k
 
        grad u(x) = grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
                 + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]

      Args:
        eps: float   
               precision requested
        sources: float(3,n)   
               source locations
        stoklet: float(nd,3,n) or float(3,n)
               Stokeslet charge strengths (sigma vectors above)
        strslet: float(nd,3,n) or float(3,n)
               stresslet strengths (mu vectors above)
        strsvec: float(nd,3,n) or float(3,n)
               stresslet orientations (nu vectors above)
        targets: float(3,nt)
               target locations (x)

        ifppreg: integer
               flag for evaluating potential, gradient, and pressure
               at the sources
               ifppreg = 1, only potential
               ifppreg = 2, potential and pressure
               ifppreg = 3, potential, pressure, and gradient

        ifppregtarg: integer
               flag for evaluating potential, gradient, and pressure
               at the targets
               ifppregtarg = 1, only potential
               ifppregtarg = 2, potential and pressure
               ifppregtarg = 3, potential, pressure, and gradient

        nd:   integer
               number of densities
        
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        out.pot: velocity at source locations if requested
        out.pre: pressure at source locations if requested
        out.grad: gradient of velocity at source locations if requested
        out.pottarg: velocity at target locations if requested
        out.pretarg: pressure at target locations if requested
        out.gradtarg: gradient of velocity at target locations if requested
              
      Example:
        see stfmmexample.py

    r"""
    out = Output()

    if(ifppreg == 0 and ifppregtarg == 0):
        print("Nothing to compute, set either ifppreg or ifppregtarg to non-zero")
        return out

    if(stoklet == None and strslet == None and strsvec == None):
        print("Nothing to compute, set either stoklet or strslet+strsvec to non-None")
        return out

    if(strslet is not None and strsvec is None):
        print("strslet and strsvec mush be both None or both not None")
        return out
    if(strslet is None and strsvec is not None):
        print("strslet and strsvec mush be both None or both not None")
        return out

    assert sources.shape[0] == 3, "The first dimension of sources must be 3"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    if(targets is not None):
        assert targets.shape[0] == 3, "The first dimension of targets must be 3"
        if(np.size(np.shape(targets))==2):
            nt = targets.shape[1]
        if(np.size(np.shape(targets))==1):
            nt = 1
    else:
        targets = np.zeros([3,0],dtype='double')
        nt = 0

    ifstoklet = 0
    ifstrslet = 0

    if(stoklet is not None):
        if(nd == 1 and ns>1):
            assert stoklet.shape[0] == 3 and stoklet.shape[1] == ns, "stoklet vectors must be of shape [3,number of sources]"
        if(nd == 1 and ns==1):
            assert stoklet.shape[0] == 3, "stoklet vectors must be of shape [3,number of sources]"
        if(nd>1):
            assert stoklet.shape[0] == nd and stoklet.shape[1] == 3 and stoklet.shape[2] == ns, "stoklet vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        stoklet = stoklet.reshape([nd,3,ns])
        ifstoklet = 1
    else:
        stoklet = np.zeros([nd,3,ns],dtype='double')

    if(strslet is not None and strsvec is not None):
        if(nd == 1 and ns>1):
            assert strslet.shape[0] == 3 and strslet.shape[1] == ns, "strslet vectors must be of shape [3,number of sources]"
            assert strsvec.shape[0] == 3 and strsvec.shape[1] == ns, "strsvec vectors must be of shape [3,number of sources]"
        if(nd == 1 and ns==1):
            assert strslet.shape[0] == 3, "strslet vectors must be of shape [3,number of sources]"
            assert strsvec.shape[0] == 3, "strsvec vectors must be of shape [3,number of sources]"
        if(nd>1):
            assert strslet.shape[0] == nd and strslet.shape[1] == 3 and strslet.shape[2] == ns, "strslet vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
            assert strsvec.shape[0] == nd and strsvec.shape[1] == 3 and strsvec.shape[2] == ns, "strsvec vectors must be of shape [nd,3,ns] where nd is number of densities, and ns is number of sources"
        strslet = strslet.reshape([nd,3,ns])
        strsvec = strsvec.reshape([nd,3,ns])
        ifstrslet = 1
    else:
        strslet = np.zeros([nd,3,ns],dtype='double')
        strsvec = np.zeros([nd,3,ns],dtype='double')

    if(ifstoklet == 1 and ifstrslet == 0): 
        out.pot,out.pre,out.grad = stfmm.st3ddirectstokg(sources,stoklet,sources,thresh,nd,ns,nt)
        out.pottarg,out.pretarg,out.gradtarg = stfmm.st3ddirectstokg(sources,stoklet,targets,thresh,nd,ns,nt)
    else:
        out.pot,out.pre,out.grad = stfmm.st3ddirectstokstrsg(sources,stoklet,1,strslet,strsvec,sources,thresh,nd,ns,nt)
        out.pottarg,out.pretarg,out.gradtarg = stfmm.st3ddirectstokstrsg(sources,stoklet,1,strslet,strsvec,targets,thresh,nd,ns,nt)

    if(ifppreg < 3):
        out.grad = None
    if(ifppregtarg < 3):
        out.gradtarg = None

    if(ifppreg < 2):
        out.pre = None
    if(ifppregtarg < 2):
        out.pretarg = None

    if(ifppreg < 1):
        out.pot = None
    if(ifppregtarg < 1):
        out.pottarg = None

    return out

def comperr(*,ntest,out,outex,pg=0,pgt=0,nd=1):
    r = 0
    err = 0
    if(nd == 1):
        if(pg > 0):
            r = r+la.norm(outex.pot[0:ntest])**2
            err = err+la.norm(outex.pot[0:ntest]-out.pot[0:ntest])**2
        if(pg >= 2):
            g = out.grad[:,0:ntest].reshape(3*ntest,)
            gex = outex.grad[:,0:ntest].reshape(3*ntest,)
            r = r +la.norm(gex)**2
            err = err+la.norm(gex-g)**2
        if( pg >= 3):
            h = out.hess[:,0:ntest].reshape(6*ntest,)
            hhex = outex.hess[:,0:ntest].reshape(6*ntest,)
            r = r + la.norm(hhex)**2
            err = err + la.norm(hhex-h)**2
        if(pgt > 0):
            r = r+la.norm(outex.pottarg[0:ntest])**2
            err = err+la.norm(outex.pottarg[0:ntest]-out.pottarg[0:ntest])**2
        if(pgt >= 2):
            g = out.gradtarg[:,0:ntest].reshape(3*ntest,)
            gex = outex.gradtarg[:,0:ntest].reshape(3*ntest,)
            r = r +la.norm(gex)**2
            err = err+la.norm(gex-g)**2
        if( pgt >= 3):
            h = out.hesstarg[:,0:ntest].reshape(6*ntest,)
            hhex = outex.hesstarg[:,0:ntest].reshape(6*ntest,)
            r = r + la.norm(hhex)**2
            err = err + la.norm(hhex-h)**2
    if(nd > 1):
        if(pg > 0):
            p = out.pot[:,0:ntest].reshape(nd*ntest,)
            pex = outex.pot[:,0:ntest].reshape(nd*ntest,)
            r = r+la.norm(pex)**2
            err = err+la.norm(p-pex)**2
        if(pg >= 2):
            g = out.grad[:,:,0:ntest].reshape(3*nd*ntest,)
            gex = outex.grad[:,:,0:ntest].reshape(3*nd*ntest,)
            r = r +la.norm(gex)**2
            err = err+la.norm(gex-g)**2
        if( pg >= 3):
            h = out.hess[:,:,0:ntest].reshape(6*nd*ntest,)
            hhex = outex.hess[:,:,0:ntest].reshape(6*nd*ntest,)
            r = r + la.norm(hhex)**2
            err = err + la.norm(hhex-h)**2
        if(pgt > 0):
            p = out.pottarg[:,0:ntest].reshape(nd*ntest,)
            pex = outex.pottarg[:,0:ntest].reshape(nd*ntest,)
            r = r+la.norm(pex)**2
            err = err+la.norm(p-pex)**2
        if(pgt >= 2):
            g = out.gradtarg[:,:,0:ntest].reshape(3*nd*ntest,)
            gex = outex.gradtarg[:,:,0:ntest].reshape(3*nd*ntest,)
            r = r +la.norm(gex)**2
            err = err+la.norm(gex-g)**2
        if( pgt >= 3):
            h = out.hesstarg[:,:,0:ntest].reshape(6*nd*ntest,)
            hhex = outex.hesstarg[:,:,0:ntest].reshape(6*nd*ntest,)
            r = r + la.norm(hhex)**2
            err = err + la.norm(hhex-h)**2
    err = np.sqrt(err/r)
    return err
