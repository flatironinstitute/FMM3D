function [U] = rfmm3d(eps,srcinfo,pg,targ,pgt)


  sources = srcinfo.sources;
  [m,ns] = size(sources);
  assert(m==3,'The first dimension of sources must be 3');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  pot = zeros(nd,1); 
  grad = zeros(nd*3,1);
  

  if(pg>=1), pot = zeros(nd,ns); end;
  if(pg == 2), grad = zeros(nd*3,ns); end;

  pottarg = zeros(nd,1);
  gradtarg = zeros(nd*3,1);
  if( nargin == 3 )
    nt = 0;
    iftarg = 0;
    pgt = 0;
    targ = zeros(3,1);
  else
    [m,nt] = size(targ);
    iftarg = 1;
    assert(m==3,'First dimension of targets must be 3');
    if(pgt >=1), pottarg = zeros(nd,nt); end;
    if(pgt == 2), gradtarg = zeros(nd*3,nt); end;
  end

  if(pg ==0 && pgt ==0), disp('Nothing to compute, set eigher pg or pgt to 1 or 2'); return; end;

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


  if(iftarg == 0 || (pgt ~=1 && pgt ~=2)) 
    if(pg == 1)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'rfmm3dpartstoscp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], io double[xx])';
[pot] = fmm3d(mex_id_, nd, eps, ns, sources, charges, pot, 1, 1, 1, 3, ns, nd, ns, nd, ns);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstosdp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], io double[xx])';
[pot] = fmm3d(mex_id_, nd, eps, ns, sources, dipoles, pot, 1, 1, 1, 3, ns, nd3, ns, nd, ns);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstoscdp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i double[xx], io double[xx])';
[pot] = fmm3d(mex_id_, nd, eps, ns, sources, charges, dipoles, pot, 1, 1, 1, 3, ns, nd, ns, nd3, ns, nd, ns);
      end
      U.pot = pot;
    end
    if(pg == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'rfmm3dpartstoscg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], io double[xx], io double[xx])';
[pot, grad] = fmm3d(mex_id_, nd, eps, ns, sources, charges, pot, grad, 1, 1, 1, 3, ns, nd, ns, nd, ns, nd3, ns);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstosdg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], io double[xx], io double[xx])';
[pot, grad] = fmm3d(mex_id_, nd, eps, ns, sources, dipoles, pot, grad, 1, 1, 1, 3, ns, nd3, ns, nd, ns, nd3, ns);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstoscdg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i double[xx], io double[xx], io double[xx])';
[pot, grad] = fmm3d(mex_id_, nd, eps, ns, sources, charges, dipoles, pot, grad, 1, 1, 1, 3, ns, nd, ns, nd3, ns, nd, ns, nd3, ns);
      end
      U.pot = pot;
      U.grad = squeeze(reshape(grad,[nd,3,ns]));
    end
  end
  if(iftarg == 1 && pg ~=1 && pg ~=2) 
    if(pgt == 1)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'rfmm3dpartstotcp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i int[x], i double[xx], io double[xx])';
[pottarg] = fmm3d(mex_id_, nd, eps, ns, sources, charges, nt, targ, pottarg, 1, 1, 1, 3, ns, nd, ns, 1, 3, nt, nd, nt);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstotdp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i int[x], i double[xx], io double[xx])';
[pottarg] = fmm3d(mex_id_, nd, eps, ns, sources, dipoles, nt, targ, pottarg, 1, 1, 1, 3, ns, nd3, ns, 1, 3, nt, nd, nt);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstotcdp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i double[xx], i int[x], i double[xx], io double[xx])';
[pottarg] = fmm3d(mex_id_, nd, eps, ns, sources, charges, dipoles, nt, targ, pottarg, 1, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt);
      end
      U.pottarg = pottarg;
    end
    if(pgt == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'rfmm3dpartstotcg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i int[x], i double[xx], io double[xx], io double[xx])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, ns, sources, charges, nt, targ, pottarg, gradtarg, 1, 1, 1, 3, ns, nd, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstotdg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i int[x], i double[xx], io double[xx], io double[xx])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, ns, sources, dipoles, nt, targ, pottarg, gradtarg, 1, 1, 1, 3, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstotcdg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i double[xx], i int[x], i double[xx], io double[xx], io double[xx])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, ns, sources, charges, dipoles, nt, targ, pottarg, gradtarg, 1, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      U.pottarg = pottarg;
      U.gradtarg = squeeze(reshape(gradtarg,[nd,3,nt]));
    end
  end
  if(iftarg == 1 && (pg ==1 || pg ==2))
    assert(pg==pgt,'pg must be pgt');
    if(pgt == 1)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'rfmm3dpartstostcp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], io double[xx], i int[x], i double[xx], io double[xx])';
[pot, pottarg] = fmm3d(mex_id_, nd, eps, ns, sources, charges, pot, nt, targ, pottarg, 1, 1, 1, 3, ns, nd, ns, nd, ns, 1, 3, nt, nd, nt);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstostdp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], io double[xx], i int[x], i double[xx], io double[xx])';
[pot, pottarg] = fmm3d(mex_id_, nd, eps, ns, sources, dipoles, pot, nt, targ, pottarg, 1, 1, 1, 3, ns, nd3, ns, nd, ns, 1, 3, nt, nd, nt);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstostcdp_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i double[xx], io double[xx], i int[x], i double[xx], io double[xx])';
[pot, pottarg] = fmm3d(mex_id_, nd, eps, ns, sources, charges, dipoles, pot, nt, targ, pottarg, 1, 1, 1, 3, ns, nd, ns, nd3, ns, nd, ns, 1, 3, nt, nd, nt);
      end
      U.pot = pot;
      U.pottarg = pottarg;
    end
    if(pgt == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'rfmm3dpartstostcg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], io double[xx], io double[xx], i int[x], i double[xx], io double[xx], io double[xx])';
[pot, grad, pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, ns, sources, charges, pot, grad, nt, targ, pottarg, gradtarg, 1, 1, 1, 3, ns, nd, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstostdg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], io double[xx], io double[xx], i int[x], i double[xx], io double[xx], io double[xx])';
[pot, grad, pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, ns, sources, dipoles, pot, grad, nt, targ, pottarg, gradtarg, 1, 1, 1, 3, ns, nd3, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'rfmm3dpartstostcdg_vec(i int[x], i double[x], i int[x], i double[xx], i double[xx], i double[xx], io double[xx], io double[xx], i int[x], i double[xx], io double[xx], io double[xx])';
[pot, grad, pottarg, gradtarg] = fmm3d(mex_id_, nd, eps, ns, sources, charges, dipoles, pot, grad, nt, targ, pottarg, gradtarg, 1, 1, 1, 3, ns, nd, ns, nd3, ns, nd, ns, nd3, ns, 1, 3, nt, nd, nt, nd3, nt);
      end
      U.pot = pot;
      U.grad = squeeze(reshape(grad,[nd,3,ns]));
      U.pottarg = pottarg;
      U.gradtarg = squeeze(reshape(gradtarg,[nd,3,nt]));
    end
  end
end

