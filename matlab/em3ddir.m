function [U] = em3ddir(zk,srcinfo,targ,ifE,ifcurlE,ifdivE)
% EM3DDIR    Slow direct Maxwell kernel sums (reference for EMFMM3D).
%
% U = em3ddir(zk,srcinfo,targ,ifE,ifcurlE,ifdivE)
%
%  Maxwell direct evaluation in R^3: evaluate all pairwise particle
%  interactions with targets. This is the slow O(N^2) direct code used
%  as a reference for testing the fast code emfmm3d.
%
%  Kernel definitions, input and outputs arguments are identical to
%  emfmm3d (see that function for all definitions), except that the first
%  argument (eps) is absent.
%
%  See also: EMFMM3D

  if(nargin<4)
    return;
  end
  if(nargin<5)
    ifcurlE = 0;
  end
  if(nargin<6)
    ifdivE = 0;
  end

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
  
  [m,nt] = size(targ);
  assert(m==3,'First dimension of targets must be 3');

  E = complex(zeros(nd*3,1)); nt_E = 1;
  curlE = complex(zeros(nd*3,1)); nt_curlE = 1;
  divE = complex(zeros(nd,1)); nt_divE = 1;
  
  if(ifE == 1), E = complex(zeros(nd*3,nt)); nt_E = nt; end;
  if(ifcurlE == 1), curlE = complex(zeros(nd*3,nt)); nt_curlE = nt; end;
  if(ifdivE == 1), divE = complex(zeros(nd,nt)); nt_divE = nt; end;

  if(ifE == 0 && ifcurlE == 0 && ifdivE == 0), disp('Nothing to compute, set eigher ifE, ifcurlE or ifdiv E to 1'); return; end;

  if(isfield(srcinfo,'e_charge'))
    ife_charge = 1;
    ns_e_charge = ns;
    e_charge = srcinfo.e_charge;
    if(nd==1), assert(length(e_charge)==ns,'Charges must be same length as second dimension of sources'); end;
    if(nd>1), [a,b] = size(e_charge); assert(a==nd && b==ns,'Charges must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end;
  else
    ife_charge = 0;
    ns_e_charge = 1;
    e_charge = complex(zeros(nd,1));
  end

  if(isfield(srcinfo,'h_current'))
    ifh_current = 1;
    ns_h_current = ns;
    h_current = srcinfo.h_current;
    if(nd == 1), [a,b] = size(squeeze(h_current)); assert(a==3 && b==ns,'h_current must be of shape[3,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(h_current); assert(a==nd && b==3 && c==ns, 'h_current must be of shape[nd,3,ns], where nd is number of densities, and ns is the number of sources'); end;
    h_current = reshape(h_current,[3*nd,ns]);
  else
    ifh_current = 0;
    ns_h_current = 1;
    h_current = complex(zeros(nd*3,1));
  end

  if(isfield(srcinfo,'e_current'))
    ife_current = 1;
    ns_e_current = ns;
    e_current = srcinfo.e_current;
    if(nd == 1), [a,b] = size(squeeze(e_current)); assert(a==3 && b==ns,'e_current must be of shape[3,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(e_current); assert(a==nd && b==3 && c==ns, 'e_current must be of shape[nd,3,ns], where nd is number of densities, and ns is the number of sources'); end;
    e_current = reshape(e_current,[3*nd,ns]);
  else
    ife_current = 0;
    ns_e_current = 1;
    e_current = complex(zeros(nd*3,1));
  end

  if(ife_charge == 0 && ife_current == 0 && ifh_current == 0), disp('Nothing to compute, set eigher e_charge, e_current or h_current'); return; end;

  nd3 = 3*nd;
  ier = 0;

  mex_id_ = 'em3ddirect(i int[x], i dcomplex[x], i int[x], i double[xx], i int[x], i dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i int[x], io dcomplex[xx], i int[x], io dcomplex[xx], i double[x])';
[E, curlE, divE] = fmm3d(mex_id_, nd, zk, ns, sources, ifh_current, h_current, ife_current, e_current, ife_charge, e_charge, nt, targ, ifE, E, ifcurlE, curlE, ifdivE, divE, thresh, 1, 1, 1, 3, ns, 1, nd3, ns_h_current, 1, nd3, ns_e_current, 1, nd, ns_e_charge, 1, 3, nt, 1, nd3, nt_E, 1, nd3, nt_curlE, 1, nd, nt_divE, 1);

  if(ifE == 1)
    U.E = squeeze(reshape(E,[nd,3,nt]));
  end
  if(ifcurlE == 1)
    U.curlE = squeeze(reshape(curlE,[nd,3,nt]));
  end
  if(ifdivE == 1)
    U.divE = divE;
  end
end


% ---------------------------------------------------------------------
