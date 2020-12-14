function [U] = h3ddir(zk,srcinfo,targ,pgt)
%
%
%  This subroutine computes the N-body Helmholtz
%  interactions and its gradients in three dimensions where 
%  the interaction kernel is given by $e^{ikr}/r$
% 
%    u(x) = \sum_{j=1}^{N} c_{j} \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|} - 
%      v_{j} \cdot \nabla \left( \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)   
%
%  where $c_{j}$ are the charge densities
%  $v_{j}$ are the dipole orientation vectors, and
%  $x_{j}$ are the source locations.
%  When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
%  from the sum.
%  
%  The sum is evaluated directly - (slow code for testing)
% 
%  Args:
%
%  -  zk: complex
%        Helmholtz parameter, k
%  -  srcinfo: structure
%        structure containing sourceinfo
%     
%     *  srcinfo.sources: double(3,n)    
%           source locations, $x_{j}$
%     *  srcinfo.nd: integer
%           number of charge/dipole vectors (optional, 
%           default - nd = 1)
%     *  srcinfo.charges: complex(nd,n) 
%           charge densities, $c_{j}$ (optional, 
%           default - term corresponding to charges dropped)
%     *  srcinfo.dipoles: complex(nd,3,n) 
%           dipole orientation vectors, $v_{j}$ (optional
%           default - term corresponding to dipoles dropped) 
%  
%  -  targ: double(3,nt)
%        target locations, $t_{i}$ 
%  -  pgt: integer
%        | target eval flag 
%        | potential at targets evaluated if pgt = 1
%        | potenial and gradient at targets evaluated if pgt=2  
%  
%  Returns:
%  
%  -  U.pottarg: potential at target locations, if requested, $u(t_{i})$
%  -  U.gradtarg: gradient at target locations, if requested, $\nabla u(t_{i})$
 


  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==3,'The first dimension of sources must be 3');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  thresh = 1e-15;

  pottarg = complex(zeros(nd,1));
  gradtarg = complex(zeros(nd*3,1));
  [m,nt] = size(targ);
  assert(m==3,'First dimension of targets must be 3');
  if(pgt >=1), pottarg = complex(zeros(nd,nt)); end;
  if(pgt == 2), gradtarg = complex(zeros(nd*3,nt)); end;

  if(pgt ==0), disp('Nothing to compute, set pgt to 1 or 2'); return; end;

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

  if(pgt == 1)
    if(ifcharge==1 && ifdipole == 0)
      mex_id_ = 'h3ddirectcp(i int64_t[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, charges, ns, targ, nt, pottarg, thresh, 1, 1, 3, ns, nd, ns, 1, 3, nt, 1, nd, nt, 1);
    end
    if(ifcharge==0 && ifdipole == 1)
      mex_id_ = 'h3ddirectdp(i int64_t[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, targ, nt, pottarg, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, nt, 1, nd, nt, 1);
    end
    if(ifcharge==1 && ifdipole == 1)
      mex_id_ = 'h3ddirectcdp(i int64_t[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, targ, nt, pottarg, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, 1, nd, nt, 1);
    end
    U.pottarg = pottarg;
  end
  if(pgt == 2)
    if(ifcharge==1 && ifdipole == 0)
      mex_id_ = 'h3ddirectcg(i int64_t[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, charges, ns, targ, nt, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd, ns, 1, 3, nt, 1, nd, nt, nd3, nt, 1);
    end
    if(ifcharge==0 && ifdipole == 1)
      mex_id_ = 'h3ddirectdg(i int64_t[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, targ, nt, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, nt, 1, nd, nt, nd3, nt, 1);
    end
    if(ifcharge==1 && ifdipole == 1)
      mex_id_ = 'h3ddirectcdg(i int64_t[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int64_t[x], i double[xx], i int64_t[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, targ, nt, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, 1, nd, nt, nd3, nt, 1);
    end
    U.pottarg = pottarg;
    U.gradtarg = squeeze(reshape(gradtarg,[nd,3,nt]));
  end
end

% ---------------------------------------------------------------------
