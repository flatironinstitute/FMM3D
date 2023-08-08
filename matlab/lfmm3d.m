function [U,varargout] = lfmm3d(eps,srcinfo,pg,varargin)
% LFMM3D    FMM in 3D for Laplace (electrostatic) kernels.
%
%  [U,varargout] = lfmm3d(eps,srcinfo,pg,varargin)
%
%  Laplace FMM in R^3: evaluate all pairwise particle
%  interactions (ignoring self-interactions) and possibly
%  interactions with targets, using the fast multipole method
%  with precision eps.
%
%  This subroutine computes the N-body Laplace
%  interactions and its gradients in three dimensions where 
%  the interaction kernel is given by $1/r$, namely
% 
%    u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - 
%      v_{j} \cdot \nabla \left( \frac{1}{\|x-x_{j}\|}\right)   
%
%  where $c_{j}$ are the charge densities,
%  $v_{j}$ are the dipole orientation vectors, and
%  $x_{j}$ are the source locations.
%
%  There are two options for the evaluation points $x$. These
%  can be the sources themselves (pg > 0; see below) or other
%  target points of interest (pgt > 0; see below). Both options
%  can be selected in one call.
%  When $x=x_{j}$, the term corresponding to $x_{j}$ is dropped
%  from the sum.
% 
%  Args:
%
%  -  eps: double   
%        precision requested
%  -  srcinfo: structure
%        structure containing the following info about the sources:     
%     *  srcinfo.sources: double(3,n)    
%           source locations, $x_{j}$
%     *  srcinfo.nd: integer
%           number of charge/dipole vectors (optional, 
%           default - nd = 1)
%     *  srcinfo.charges: double(nd,n) 
%           charge densities, $c_{j}$ (optional, 
%           default - term corresponding to charges dropped)
%     *  srcinfo.dipoles: double(nd,3,n) 
%           dipole orientation vectors, $v_{j}$ (optional
%           default - term corresponding to dipoles dropped) 
%  -  pg: integer
%        | source eval flag
%        | potential at sources evaluated if pg = 1
%        | potential and gradient at sources evaluated if pg=2
%        | potential, gradient and hessian at sources evaluated if pg=3
%        
%  Optional args
%  -  targ: double(3,nt)
%        target locations, $t_{i}$ 
%  -  pgt: integer
%        | target eval flag 
%        | potential at targets evaluated if pgt = 1
%        | potential and gradient at targets evaluated if pgt=2 
%  -  opts: options structure, values in brackets indicate default
%           values wherever applicable
%        opts.ndiv: set number of points for subdivision criterion
%        opts.idivflag: set subdivision criterion (0)
%           opts.idivflag = 0, subdivide on sources only
%           opts.idivflag = 1, subdivide on targets only
%           opts.idivflag = 2, subdivide on sources and targets
%        opts.ifnear: include near (list 1) interactions (true)
%  
%  Returns:
%  
%  -  U.pot:      potential at source locations, if requested, $u(x_{j})$
%  -  U.grad:     gradient at source locations, if requested, $\nabla u(x_{j})$
%  -  U.hess:     hessian at source locations, if requested, $\nabla^2 u(x_{j})$
%  -  U.pottarg:  potential at target locations, if requested, $u(t_{i})$
%  -  U.gradtarg: gradient at target locations, if requested, $\nabla u(t_{i})$
%  -  U.hesstarg: hessian at target locations, if requested, $\nabla^2 u(t_{i})$
%
%  - ier: error code for FMM run
%  - timeinfo: time taken in each step of the FMM
%       timeinfo(1): form multipole step
%       timeinfo(2): multipole->multipole translation step
%       timeinfo(3): multipole to local translation, form local + multipole eval step
%       timeinfo(4): local->local translation step
%       timeinfo(5): local eval step
%       timeinfo(6): direct evaluation step
%
%
%  Examples:
%  U = lfmm3d(eps,srcinfo,pg)
%     Call the FMM for sources only with default arguments
%  U = lfmm3d(eps,srcinfo,pg,targ,pgt)
%     Call the FMM for sources + targets with default arguments
%  U = lfmm3d(eps,srcinfo,pg,opts)
%     Call the FMM for sources only with user specified arguments
%  U = lfmm3d(eps,srcinfo,pg,targ,pgt)
%     Call the FMM for sources + targets with user specified arguments 
%  [U,ier] = lfmm3d(eps,srcinfo,pg)
%     Call the FMM for sources only with default arguments and returns
%     the error code for the FMM as well
%  [U,ier,timeinfo] = lfmm3d(eps,srcinfo,pg)
%     Call the FMM for sources only with default arguments, returns
%     the error code for the FMM as well and the time split
%      
% See also: L3DDIR


  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==3,'The first dimension of sources must be 3');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  pot = zeros(nd,ns); 
  grad = zeros(nd*3,ns);
  hess = zeros(nd*6,ns);
  
  if( nargin < 3)
    disp('Not enough input arguments, exiting\n');
    return;
  end
  if( nargin == 3 )
    nt = 0;
    pgt = 0;
    targ = zeros(3,1);
    opts = [];
  elseif (nargin == 4)
    nt = 0;
    pgt = 0;
    targ = zeros(3,1);
    opts = varargin{1};
  elseif (nargin == 5)
    targ = varargin{1};
    pgt = varargin{2};
    [m,nt] = size(targ);
    assert(m==3,'First dimension of targets must be 3');
    opts = [];
  elseif (nargin == 6)
    targ = varargin{1};
    pgt = varargin{2};
    [m,nt] = size(targ);
    assert(m==3,'First dimension of targets must be 3');
    opts = varargin{3};
  end
  ntuse = max(nt,1);
  pottarg = zeros(nd,ntuse);
  gradtarg = zeros(nd*3,ntuse);
  hesstarg = zeros(nd*6,ntuse);


  if((pg ==0 && pgt ==0) || (ns == 0)), disp('Nothing to compute, set eigher pg or pgt to 1 or 2'); return; end;

  if(isfield(srcinfo,'charges'))
    ifcharge = 1;
    charges = srcinfo.charges;
    if(nd==1), assert(length(charges)==ns,'Charges must be same length as second dimension of sources'); end;
    if(nd>1), [a,b] = size(charges); assert(a==nd && b==ns,'Charges must be of shape [nd,ns] where nd is the number of densities, and ns is the number of sources'); end;
  else
    ifcharge = 0;
    charges = zeros(nd,ns);
  end

  if(isfield(srcinfo,'dipoles'))
    ifdipole = 1;
    dipoles = srcinfo.dipoles;
    if(nd == 1), [a,b] = size(squeeze(dipoles)); assert(a==3 && b==ns,'Dipoles must be of shape[3,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(dipoles); assert(a==nd && b==3 && c==ns, 'Dipoles must be of shape[nd,3,ns], where nd is number of densities, and ns is the number of sources'); end;
    dipoles = reshape(dipoles,[3*nd,ns]);
  else
    ifdipole = 0;
    dipoles = zeros(nd*3,ns);
  end

  nd3 = 3*nd;
  nd6 = 6*nd;
  ier = 0;

  ndiv = 400;
  idivflag = 0;
  mex_id_ = 'lndiv(i double[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io int[x])';
[ndiv, idivflag] = fmm3d(mex_id_, eps, ns, nt, ifcharge, ifdipole, pg, pgt, ndiv, idivflag, 1, 1, 1, 1, 1, 1, 1, 1, 1);
  if(isfield(opts,'ndiv'))
    ndiv = opts.ndiv;
  end

  if(isfield(opts,'idivflag'))
    idivflag = opts.idivflag;
  end

  ifnear = 1;
  if(isfield(opts,'ifnear'))
    ifnear = opts.ifnear;
  end
  iper = 1;
  timeinfo = zeros(6,1);
  mex_id_ = 'lfmm3d_ndiv(i int[x], i double[x], i int[x], i double[xx], i int[x], i double[xx], i int[x], i double[xx], i int[x], i int[x], io double[xx], io double[xx], io double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i int[x], i int[x], i int[x], io double[x], io int[x])';
[pot, grad, hess, pottarg, gradtarg, hesstarg, timeinfo, ier] = fmm3d(mex_id_, nd, eps, ns, sources, ifcharge, charges, ifdipole, dipoles, iper, pg, pot, grad, hess, nt, targ, pgt, pottarg, gradtarg, hesstarg, ndiv, idivflag, ifnear, timeinfo, ier, 1, 1, 1, 3, ns, 1, nd, ns, 1, nd3, ns, 1, 1, nd, ns, nd3, ns, nd6, ns, 1, 3, ntuse, 1, nd, ntuse, nd3, ntuse, nd6, ntuse, 1, 1, 1, 6, 1);

  U.pot = [];
  U.grad = [];
  U.hess = [];
  U.pottarg = [];
  U.gradtarg = [];
  U.hesstarg = [];
  if(pg >= 1), U.pot = squeeze(reshape(pot,[nd,ns])); end;
  if(pg >= 2), U.grad = squeeze(reshape(grad,[nd,3,ns])); end;
  if(pg >= 3), U.hess = squeeze(reshape(hess,[nd,6,ns])); end;
  if(pgt >= 1), U.pottarg = squeeze(reshape(pottarg,[nd,nt])); end;
  if(pgt >= 2), U.gradtarg = squeeze(reshape(gradtarg,[nd,3,nt])); end;
  if(pgt >= 3), U.hesstarg = squeeze(reshape(hesstarg,[nd,6,nt])); end;

  varargout{1} = ier;
  varargout{2} = timeinfo;
end

% ---------------------------------------------------------------------
