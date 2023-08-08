function [U] = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE)
% EMFMM3D   FMM in 3D for Maxwell (frequency-domain electromagnetic) kernels.
%
%  U = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE)
%
%  Frequency-domain Maxwell FMM in R^3: evaluate all pairwise particle
%  interactions (ignoring self-interactions) and
%  interactions with targets, using the fast multipole method
%  with precision eps.
%
%  Specifically, this subroutine computes a sum for the electric field
%
%      E(x) = sum_m curl G_k(x,y^{(m)}) h_current_m
%                 + G_k(x,y^{(m)}) e_current_m 
%                 + grad G_k(x,y^{(m)}) e_charge_m
%
%  for each requested evaluation point x, where h_current and e_current
%  are 3-vector densities and e_charge is a scalar density supplied
%  at each source point y^{(m)}. G_k is the Helmholtz Green function
%  without the 1/(4pi) scaling:
%
%      G_k(x,y) = e^(ik|x-y|)/|x-y|.
%
%  In contrast with other FMM routines in the library, this routine
%  has only 1 option for the evaluation points: they are specified
%  as targets. If a target x coincides with a source point y^{(m)}
%  that term in the sum is omitted. 
%
%  The electric field is, naturally, a 3-vector at each point x.
%  With appropriate input flags, the subroutine also computes divE,
%  curlE. 
%
%  Remark: the subroutine uses a stabilized representation
%  for computing the divergence by using integration by parts
%  wherever possible. If the divergence is not requested, then the
%  Helmholtz FMM is called with 3*nd densities, while if the divergence
%  is requested, then the Helmholtz FMM is called with 4*nd densities
% 
%  Args:
%
%  -  eps: double   
%        precision requested
%  -  zk: complex
%        Helmholtz parameter, k
%  -  srcinfo: structure
%        structure containing the following info about the sources:
%     *  srcinfo.sources: double(3,n)    
%           source locations, $x_{j}$
%     *  srcinfo.nd: integer
%           number of charge/dipole vectors (optional, 
%           default - nd = 1)
%     *  srcinfo.h_current: complex(nd,3,n) 
%           a vector source (optional,
%           default - term corresponding to h_current dropped) 
%     *  srcinfo.e_current: complex(nd,3,n) 
%           b vector source (optional,
%           default - term corresponding to e_current dropped) 
%     *  srcinfo.e_charge: complex(nd,n) 
%           e_charge source (optional, 
%           default - term corresponding to e_charge dropped)
%  -  targ: double(3,nt)
%        target locations, $t_{i}$
%  -  ifE: integer
%        E is returned at the target locations if ifE = 1
%  -  ifcurlE: integer
%        curl E is returned at the target locations if ifcurlE = 1
%  -  ifdivE: integer
%        div E is returned at the target locations if ifdivE = 1
%
%  Returns:
%  
%  -  U.E:     E field defined in above equation at targets if requested
%  -  U.curlE: curl of E field at target locations if requested
%  -  U.divE:  divergence of E at target locations if requested
%
% See also: HFMM3D, EM3DDIR

  if(nargin<5)
    return;
  end
  if(nargin<6)
    ifcurlE = 0;
  end
  if(nargin<7)
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

  mex_id_ = 'emfmm3d(i int[x], i double[x], i dcomplex[x], i int[x], i double[xx], i int[x], i dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i int[x], io dcomplex[xx], i int[x], io dcomplex[xx], io int[x])';
[E, curlE, divE, ier] = fmm3d(mex_id_, nd, eps, zk, ns, sources, ifh_current, h_current, ife_current, e_current, ife_charge, e_charge, nt, targ, ifE, E, ifcurlE, curlE, ifdivE, divE, ier, 1, 1, 1, 1, 3, ns, 1, nd3, ns_h_current, 1, nd3, ns_e_current, 1, nd, ns_e_charge, 1, 3, nt, 1, nd3, nt_E, 1, nd3, nt_curlE, 1, nd, nt_divE, 1);

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
