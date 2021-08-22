function [U] = st3ddir(srcinfo,targ,ifppregtarg)
%
%
%  Stokes FMM in R^{3}: evaluate all pairwise particle
%  interactions (ignoring self-interactions) and
%  interactions with targs.
%
%  This routine computes sums of the form
%
%  u(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
%       + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
%
%  where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
%  stresslet charge, and nu^{(m)} is the stresslet orientation
%  (note that each of these is a 3 vector per source point y^{(m)}).
%  For x a source point, the self-interaction in the sum is omitted.
%
%  Optionally, the associated pressure p(x) and gradient grad u(x)
%  are returned
%
%    p(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
%         + sum_m T_{ijk}(x,y^{(m)}) PI_{jk} mu^{(m)}_j nu^{(m)}_k
%
%    grad u(x) = grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
%              + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
% 
%  Args:
%
%  -  srcinfo: structure
%        structure containing sourceinfo
%     
%     *  srcinfo.sources: double(3,n)    
%           source locations, $x_{j}$
%     *  srcinfo.nd: integer
%           number of densities (optional, 
%           default - nd = 1)
%     *  srcinfo.stoklet: double(nd,3,n) 
%           Stokeslet charge strengths, $sigma_{j}$ (optional, 
%           default - term corresponding to Stokeslet charge strengths dropped)
%     *  srcinfo.strslet: double(nd,3,n) 
%           stresslet strengths, $mu_{j}$ (optional
%           default - term corresponding to stresslet strengths dropped) 
%     *  srcinfo.strsvec: double(nd,3,n) 
%           stresslet orientations, $nu_{j}$ (optional
%           default - term corresponding to stresslet orientations dropped) 
%
%  -  targ: double(3,nt)
%        target locations, $t_{i}$ (optional)
%  -  ifppregtarg: integer
%        | target eval flag (optional)
%        | potential at targets evaluated if ifppregtarg = 1
%        | potential and pressure at targets evaluated if ifppregtarg = 2 
%        | potential, pressure and gradient at targets evaluated if ifppregtarg = 3
%  
%  Returns:
%  
%  -  U.pottarg: velocity at target locations if requested
%  -  U.pretarg: pressure at target locations if requested
%  -  U.gradtarg: gradient of velocity at target locations if requested

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

  if( nargin <= 1 )
    return;
  else
    if( nargin <= 2 ), ifppregtarg = 3; end;
    [m,nt] = size(targ);
    assert(m==3,'First dimension of targets must be 3');
    pottarg = zeros(nd*3,nt);
    pretarg = zeros(nd,nt);
    gradtarg = zeros(nd*9,nt);
  end

  if(ifppregtarg == 0), disp('Nothing to compute, set eigher ifppregtarg to 1 or 2 or 3'); return; end;

  if(isfield(srcinfo,'stoklet'))
    ifstoklet = 1;
    ns_stok = ns;
    stoklet = srcinfo.stoklet;
    if(nd == 1), [a,b] = size(squeeze(stoklet)); assert(a==3 && b==ns,'Stoklet must be of shape[3,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(stoklet); assert(a==nd && b==3 && c==ns, 'Stoklet must be of shape[nd,3,ns], where nd is number of densities, and ns is the number of sources'); end;
    stoklet = reshape(stoklet,[3*nd,ns]);
  else
    ifstoklet = 0;
    ns_stok = 1;
    stoklet = zeros(nd*3,1);
  end

  if(isfield(srcinfo,'strslet') && isfield(srcinfo,'strsvec'))
    ifstrslet = 1;
    ns_strs = ns;
    strslet = srcinfo.strslet;
    strsvec = srcinfo.strsvec;
    if(nd == 1), [a,b] = size(squeeze(strslet)); assert(a==3 && b==ns,'Strslet must be of shape[3,ns], where ns is the number of sources'); end;
    if(nd == 1), [a,b] = size(squeeze(strsvec)); assert(a==3 && b==ns,'Strsvec must be of shape[3,ns], where ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(strslet); assert(a==nd && b==3 && c==ns, 'Strslet must be of shape[nd,3,ns], where nd is number of densities, and ns is the number of sources'); end;
    if(nd>1), [a,b,c] = size(strsvec); assert(a==nd && b==3 && c==ns, 'Strsvec must be of shape[nd,3,ns], where nd is number of densities, and ns is the number of sources'); end;
    strslet = reshape(strslet,[3*nd,ns]);
    strsvec = reshape(strsvec,[3*nd,ns]);
  else
    ifstrslet = 0;
    ns_strs = 1;
    strslet = zeros(nd*3,1);
    strsvec = zeros(nd*3,1);
  end

  nd3 = 3*nd;
  nd9 = 9*nd;
  ier = 0;

  if(ifstoklet == 1 && ifstrslet == 0)
    mex_id_ = 'st3ddirectstokg(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pottarg, pretarg, gradtarg] = fmm3d(mex_id_, nd, sources, stoklet, ns, targ, nt, pottarg, pretarg, gradtarg, thresh, 1, 3, ns, nd3, ns_stok, 1, 3, nt, 1, nd3, nt, nd, nt, nd9, nt, 1);
  else
    istress = 1;
    mex_id_ = 'st3ddirectstokstrsg(i int[x], i double[xx], i double[xx], i int[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i double[x])';
[pottarg, pretarg, gradtarg] = fmm3d(mex_id_, nd, sources, stoklet, istress, strslet, strsvec, ns, targ, nt, pottarg, pretarg, gradtarg, thresh, 1, 3, ns, nd3, ns_stok, 1, nd3, ns_strs, nd3, ns_strs, 1, 3, nt, 1, nd3, nt, nd, nt, nd9, nt, 1);
  end

  U.pottarg = [];
  U.pretarg = [];
  U.gradtarg = [];
  if(ifppregtarg >= 1), U.pottarg = squeeze(reshape(pottarg,[nd,3,nt])); end;
  if(ifppregtarg >= 2), U.pretarg = pretarg; end;
  if(ifppregtarg >= 3), U.gradtarg = squeeze(reshape(gradtarg,[nd,3,3,nt])); end;
end
