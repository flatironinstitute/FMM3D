function [U] = hfmm3d(eps,zk,srcinfo,pg,targ,pgt)
%hfmm3d
%******
%
%   This subroutine computes the N-body Helmholtz
%   interactions and its gradients in three dimensions where 
%   the interaction kernel is given by e^{ikr}/r
% 
%   .. math::
%
%       u(x) = \sum_{j=1}^{N} c_{j} e^{ik|x-x_{j}|}/|x-x_{j}| - \nabla (e^{ik|x-x_{j}|}/|x-x_{j}|) \cdot v_{j} 
%
%   where $c_{j}$ are the charge densities
%   $v_{j}$ are the dipole orientation vectors, and
%   $x_{j}$ are the source locations.
%
%   When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped
%   from the sum
%
%   Args:
%
%   -    eps: double   
%              precision requested
%   -    zk: complex
%              Helmholtz parameter - k
%   -    srcinfo: structure
%               structure containing sourceinfo
%   -    srcinfo.sources: double(3,n)    
%              source locations ($x_{j}$)
%   -    srcinfo.nd: integer
%              number of charge/dipole densities (optional)
%              default - nd = 1
%   -    srcinfo.charges: complex(nd,n) 
%              charge densities ($c_{j}$) (optional)
%              default - term corresponding to charges dropped
%   -    srcinfo.dipoles: complex(nd,3,n) 
%               dipole orientation vectors ($v_{j}$) (optional)
%               default - term corresponding to dipoles dropped 
%   -    pg:  integer
%              source eval flag
%              potential at sources evaluated if pg = 1
%              potenial and gradient at sources evaluated if pg=2
%
%   -    targ: double(3,nt)
%               target locations (x) (optional)
%   -    pgt:  integer
%              target eval flag
%              potential at targets evaluated if pgt = 1
%              potenial and gradient at targets evaluated if pgt=2
%              (optional)
%
%   Returns:
%
%   -    U.pot  - potential at source locations if requested
%   -    U.grad - gradient at source locations if requested
%   -    U.pottarg  - potential at target locations if requested
%   -    U.gradtarg - gradient at target locations if requested
%              
%   Example:
%
%   -    see ``hfmmexample.m''

  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==3,'The first dimension of sources must be 3');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  pot = complex(zeros(nd,1)); 
  grad = complex(zeros(nd*3,1));
  

  if(pg>=1), pot = complex(zeros(nd,ns)); end;
  if(pg == 2), grad = complex(zeros(nd*3,ns)); end;

  pottarg = complex(zeros(nd,1));
  gradtarg = complex(zeros(nd*3,1));
  if( nargin == 4 )
    nt = 0;
    iftarg = 0;
    pgt = 0;
    targ = zeros(3,1);
  else
    [m,nt] = size(targ);
    iftarg = 1;
    assert(m==3,'First dimension of targets must be 3');
    if(pgt >=1), pottarg = complex(zeros(nd,nt)); end;
    if(pgt == 2), gradtarg = complex(zeros(nd*3,nt)); end;
  end

  if(pg ==0 && pgt ==0), disp('Nothing to compute, set eigher pg or pgt to 1 or 2'); return; end;

  if(isfield(srcinfo,'charges'))
    ifcharge = 1;
    charges = srcinfo.charges;
    if(nd==1), assert(length(charges)==ns,'Charges must be same length as second dimension of sources'); end;
    if(nd>1), [a,b] = size(charges); assert(a==nd && b==ns,'Charges must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end;
  else
    ifcharge = 0;
    charges = complex(zeros(nd,1));
  end

  if(isfield(srcinfo,'dipoles'))
    ifdipole = 1;
    dipoles = srcinfo.dipoles;
    if(nd == 1), [a,b] = size(squeeze(dipoles)); assert(a==3 && b==ns,'Dipoles must be of shape[3,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(dipoles); assert(a==nd && b==3 && c==ns, 'Dipoles must be of shape[nd,3,ns], where nd is number of densities, and ns is the number of sources'); end;
    dipoles = reshape(dipoles,[3*nd,ns]);
  else
    ifdipole = 0;
    dipoles = complex(zeros(nd*3,1));
  end

  nd3 = 3*nd;


  if(iftarg == 0 || (pgt ~=1 && pgt ~=2)) 
    if(pg == 1)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'hfmm3d_s_c_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], io dcomplex[xx])';
[pot] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, pot, 1, 1, 1, 1, 3, ns, nd, ns, nd, ns);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'hfmm3d_s_d_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], io dcomplex[xx])';
[pot] = fmm3d(mex_id_, nd, eps, zk, ns, sources, dipoles, pot, 1, 1, 1, 1, 3, ns, nd3, ns, nd, ns);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'hfmm3d_s_cd_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i dcomplex[xx], io dcomplex[xx])';
[pot] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, dipoles, pot, 1, 1, 1, 1, 3, ns, nd, ns, nd3, ns, nd, ns);
      end
      U.pot = pot;
    end
    if(pg == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'hfmm3d_s_c_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], io dcomplex[xx], io dcomplex[xx])';
[pot, grad] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, pot, grad, 1, 1, 1, 1, 3, ns, nd, ns, nd, ns, nd3, ns);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'hfmm3d_s_d_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], io dcomplex[xx], io dcomplex[xx])';
[pot, grad] = fmm3d(mex_id_, nd, eps, zk, ns, sources, dipoles, pot, grad, 1, 1, 1, 1, 3, ns, nd3, ns, nd, ns, nd3, ns);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'hfmm3d_s_cd_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i dcomplex[xx], io dcomplex[xx], io dcomplex[xx])';
[pot, grad] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, dipoles, pot, grad, 1, 1, 1, 1, 3, ns, nd, ns, nd3, ns, nd, ns, nd3, ns);
      end
      U.pot = pot;
      U.grad = squeeze(reshape(grad,[nd,3,ns]));
    end
  end
  if(iftarg == 1 && pg ~=1 && pg ~=2) 
    if(pgt == 1)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'hfmm3d_t_c_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], io dcomplex[xx])';
[pottarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, nt, targ, pottarg, 1, 1, 1, 1, 3, ns, nd, ns, 1, 3, nt, nd, nt);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'hfmm3d_t_d_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], io dcomplex[xx])';
[pottarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, dipoles, nt, targ, pottarg, 1, 1, 1, 1, 3, ns, nd3, ns, 1, 3, nt, nd, nt);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'hfmm3d_t_cd_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], io dcomplex[xx])';
[pottarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, dipoles, nt, targ, pottarg, 1, 1, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt);
      end
      U.pottarg = pottarg;
    end
    if(pgt == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'hfmm3d_t_c_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], io dcomplex[xx], io dcomplex[xx])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, nt, targ, pottarg, gradtarg, 1, 1, 1, 1, 3, ns, nd, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'hfmm3d_t_d_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], io dcomplex[xx], io dcomplex[xx])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, dipoles, nt, targ, pottarg, gradtarg, 1, 1, 1, 1, 3, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'hfmm3d_t_cd_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], io dcomplex[xx], io dcomplex[xx])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, dipoles, nt, targ, pottarg, gradtarg, 1, 1, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      U.pottarg = pottarg;
      U.gradtarg = squeeze(reshape(gradtarg,[nd,3,nt]));
    end
  end
  if(iftarg == 1 && (pg ==1 || pg ==2))
    assert(pg==pgt,'pg must be pgt');
    if(pgt == 1)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'hfmm3d_st_c_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], io dcomplex[xx], i int[x], i double[xx], io dcomplex[xx])';
[pot, pottarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, pot, nt, targ, pottarg, 1, 1, 1, 1, 3, ns, nd, ns, nd, ns, 1, 3, nt, nd, nt);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'hfmm3d_st_d_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], io dcomplex[xx], i int[x], i double[xx], io dcomplex[xx])';
[pot, pottarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, dipoles, pot, nt, targ, pottarg, 1, 1, 1, 1, 3, ns, nd3, ns, nd, ns, 1, 3, nt, nd, nt);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'hfmm3d_st_cd_p_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i dcomplex[xx], io dcomplex[xx], i int[x], i double[xx], io dcomplex[xx])';
[pot, pottarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, dipoles, pot, nt, targ, pottarg, 1, 1, 1, 1, 3, ns, nd, ns, nd3, ns, nd, ns, 1, 3, nt, nd, nt);
      end
      U.pot = pot;
      U.pottarg = pottarg;
    end
    if(pgt == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'hfmm3d_st_c_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i int[x], i double[xx], io dcomplex[xx], io dcomplex[xx])';
[pot, grad, pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, pot, grad, nt, targ, pottarg, gradtarg, 1, 1, 1, 1, 3, ns, nd, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'hfmm3d_st_d_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i int[x], i double[xx], io dcomplex[xx], io dcomplex[xx])';
[pot, grad, pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, dipoles, pot, grad, nt, targ, pottarg, gradtarg, 1, 1, 1, 1, 3, ns, nd3, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'hfmm3d_st_cd_g_vec(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i dcomplex[xx], i dcomplex[xx], io dcomplex[xx], io dcomplex[xx], i int[x], i double[xx], io dcomplex[xx], io dcomplex[xx])';
[pot, grad, pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, zk, ns, sources, charges, dipoles, pot, grad, nt, targ, pottarg, gradtarg, 1, 1, 1, 1, 3, ns, nd, ns, nd3, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      U.pot = pot;
      U.grad = squeeze(reshape(grad,[nd,3,ns]));
      U.pottarg = pottarg;
      U.gradtarg = squeeze(reshape(gradtarg,[nd,3,nt]));
    end
  end
end

% ---------------------------------------------------------------------
