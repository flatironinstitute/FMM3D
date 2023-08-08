function [U] = l3ddir(srcinfo,targ,pgt)
% L3DDIR    Direct (slow) 3D Laplace kernel sums (reference for LFMM3D).
%
% U = l3ddir(srcinfo,targ,pgt)
%
%  Laplace direct evaluation in R^3: evaluate all pairwise particle
%  interactions with targets. This is the slow O(N^2) direct code used
%  as a reference for testing the (fast) code lfmm3d.
%
%  Kernel definitions, input and outputs arguments are identical to
%  lfmm3d (see that function for all definitions), apart from:
%  1) the first argument (eps) is absent.
%  2) there are currently no outputs at sources, meaning that U.pot, U.grad,
%     and U.hess are missing (as if pg=0). In other words,
%     just targets for now, and targ is thus not an optional argument.
%
% See also: LFMM3D

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

  pottarg = zeros(nd,1);
  gradtarg = zeros(nd*3,1);
  hesstarg = zeros(nd*6,1);
  [m,nt] = size(targ);
  assert(m==3,'First dimension of targets must be 3');
  if(pgt >=1), pottarg = zeros(nd,nt); end;
  if(pgt >= 2), gradtarg = zeros(nd*3,nt); end;
  if(pgt >= 3), hesstarg = zeros(nd*6,nt); end;

  if(pgt ==0), disp('Nothing to compute, set eigher pgt to 1 or 2'); return; end;

  if(isfield(srcinfo,'charges'))
    ifcharge = 1;
    charges = srcinfo.charges;
    if(nd==1), assert(length(charges)==ns,'Charges must be same length as second dimension of sources'); end;
    if(nd>1), [a,b] = size(charges); assert(a==nd && b==ns,'Charges must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end;
  else
    ifcharge = 0;
    charges = zeros(nd,1);
  end

  if(isfield(srcinfo,'dipoles'))
    ifdipole = 1;
    dipoles = srcinfo.dipoles;
    if(nd == 1), [a,b] = size(squeeze(dipoles)); assert(a==3 && b==ns,'Dipoles must be of shape[3,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(dipoles); assert(a==nd && b==3 && c==ns, 'Dipoles must be of shape[nd,3,ns], where nd is number of densities, and ns is the number of sources'); end;
    dipoles = reshape(dipoles,[3*nd,ns]);
  else
    ifdipole = 0;
    dipoles = zeros(nd*3,1);
  end

  nd3 = 3*nd;
  nd6 = 6*nd;

  if(pgt == 1)
    if(ifcharge==1 && ifdipole == 0)
      mex_id_ = 'l3ddirectcp(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, sources, charges, ns, targ, nt, pottarg, thresh, 1, 3, ns, nd, ns, 1, 3, nt, 1, nd, nt, 1);
    end
    if(ifcharge==0 && ifdipole == 1)
      mex_id_ = 'l3ddirectdp(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, sources, dipoles, ns, targ, nt, pottarg, thresh, 1, 3, ns, nd3, ns, 1, 3, nt, 1, nd, nt, 1);
    end
    if(ifcharge==1 && ifdipole == 1)
      mex_id_ = 'l3ddirectcdp(i int[x], i double[xx], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, sources, charges, dipoles, ns, targ, nt, pottarg, thresh, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, 1, nd, nt, 1);
    end
    U.pottarg = pottarg;
  end
  if(pgt == 2)
    if(ifcharge==1 && ifdipole == 0)
      mex_id_ = 'l3ddirectcg(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, sources, charges, ns, targ, nt, pottarg, gradtarg, thresh, 1, 3, ns, nd, ns, 1, 3, nt, 1, nd, nt, nd3, nt, 1);
    end
    if(ifcharge==0 && ifdipole == 1)
      mex_id_ = 'l3ddirectdg(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, sources, dipoles, ns, targ, nt, pottarg, gradtarg, thresh, 1, 3, ns, nd3, ns, 1, 3, nt, 1, nd, nt, nd3, nt, 1);
    end
    if(ifcharge==1 && ifdipole == 1)
      mex_id_ = 'l3ddirectcdg(i int[x], i double[xx], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, sources, charges, dipoles, ns, targ, nt, pottarg, gradtarg, thresh, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, 1, nd, nt, nd3, nt, 1);
    end
    U.pottarg = pottarg;
    U.gradtarg = squeeze(reshape(gradtarg,[nd,3,nt]));
  end
  if(pgt == 3)
    if(ifcharge==1 && ifdipole == 0)
      mex_id_ = 'l3ddirectch(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pottarg, gradtarg, hesstarg] = fmm3d(mex_id_, nd, sources, charges, ns, targ, nt, pottarg, gradtarg, hesstarg, thresh, 1, 3, ns, nd, ns, 1, 3, nt, 1, nd, nt, nd3, nt, nd6, nt, 1);
    end
    if(ifcharge==0 && ifdipole == 1)
      mex_id_ = 'l3ddirectdh(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pottarg, gradtarg, hesstarg] = fmm3d(mex_id_, nd, sources, dipoles, ns, targ, nt, pottarg, gradtarg, hesstarg, thresh, 1, 3, ns, nd3, ns, 1, 3, nt, 1, nd, nt, nd3, nt, nd6, nt, 1);
    end
    if(ifcharge==1 && ifdipole == 1)
      mex_id_ = 'l3ddirectcdh(i int[x], i double[xx], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pottarg, gradtarg, hesstarg] = fmm3d(mex_id_, nd, sources, charges, dipoles, ns, targ, nt, pottarg, gradtarg, hesstarg, thresh, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, 1, nd, nt, nd3, nt, nd6, nt, 1);
    end
    U.pottarg = pottarg;
    U.gradtarg = squeeze(reshape(gradtarg,[nd,3,nt]));
    U.hesstarg = squeeze(reshape(hesstarg,[nd,6,nt]));
  end
end

% ---------------------------------------------------------------------
