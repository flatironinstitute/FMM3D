function [U] = st3ddir(srcinfo,targ,ifppregtarg)
% ST3DDIR    Direct (slow) 3D Stokes kernel sums (reference for STFMM3D).
%
%  U = st3ddir(srcinfo,targ,ifppregtarg)
%
%  Stokes direct evaluation in R^3: evaluate all pairwise particle
%  interactions with targets. This is the slow O(N^2) direct code used
%  as a reference for testing the (fast) code stfmm3d.
%
%  Kernel definitions, input and outputs arguments are identical to
%  stfmm3d (see that function for all definitions), apart from:
%  1) the first argument (eps) is absent.
%  2) there are currently no outputs at sources, meaning that U.pot, U.pre,
%     and U.grad are missing (as if ifppreg=0). In other words,
%     just targets for now, and targ is thus not an optional argument.
%
% See also: STFMM3D

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
