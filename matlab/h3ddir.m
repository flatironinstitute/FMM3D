function [U] = h3ddir(eps,zk,ntest,srcinfo,pg,targ,pgt)

  sources = srcinfo.sources;
  stmp = sources(:,1:ntest);
  [m,ns] = size(sources);
  assert(m==3,'The first dimension of sources must be 3');
  if(~isfield(srcinfo,'nd'))
    nd = 1;
  end
  if(isfield(srcinfo,'nd'))
    nd = srcinfo.nd;
  end

  thresh = 1e-15;

  pot = complex(zeros(nd,1)); 
  grad = complex(zeros(nd*3,1));
  

  if(pg>=1), pot = complex(zeros(nd,ntest)); end;
  if(pg == 2), grad = complex(zeros(nd*3,ntest)); end;

  pottarg = complex(zeros(nd,1));
  gradtarg = complex(zeros(nd*3,1));
  if( nargin == 5 )
    nt = 0;
    iftarg = 0;
    pgt = 0;
    targ = zeros(3,1);
  else
    [m,nt] = size(targ);
    iftarg = 1;
    assert(m==3,'First dimension of targets must be 3');
    ttmp = targ(:,1:ntest);
    if(pgt >=1), pottarg = complex(zeros(nd,ntest)); end;
    if(pgt == 2), gradtarg = complex(zeros(nd*3,ntest)); end;
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
        mex_id_ = 'h3ddirectcp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pot] = fmm3d(mex_id_, nd, zk, sources, charges, ns, stmp, ntest, pot, thresh, 1, 1, 3, ns, nd, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'h3ddirectdp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pot] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, stmp, ntest, pot, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'h3ddirectcdp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pot] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, stmp, ntest, pot, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      U.pot = pot;
    end
    if(pg == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'h3ddirectcg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pot, grad] = fmm3d(mex_id_, nd, zk, sources, charges, ns, stmp, ntest, pot, grad, thresh, 1, 1, 3, ns, nd, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'h3ddirectdg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pot, grad] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, stmp, ntest, pot, grad, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'h3ddirectcdg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pot, grad] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, stmp, ntest, pot, grad, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      U.pot = pot;
      U.grad = squeeze(reshape(grad,[nd,3,ntest]));
    end
  end
  if(iftarg == 1 && pg ~=1 && pg ~=2) 
    if(pgt == 1)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'h3ddirectcp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, charges, ns, ttmp, ntest, pottarg, thresh, 1, 1, 3, ns, nd, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'h3ddirectdp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, ttmp, ntest, pottarg, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'h3ddirectcdp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, ttmp, ntest, pottarg, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      U.pottarg = pottarg;
    end
    if(pgt == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'h3ddirectcg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, charges, ns, ttmp, ntest, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'h3ddirectdg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, ttmp, ntest, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'h3ddirectcdg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, ttmp, ntest, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      U.pottarg = pottarg;
      U.gradtarg = squeeze(reshape(gradtarg,[nd,3,ntest]));
    end
  end
  if(iftarg == 1 && (pg ==1 || pg ==2))
    assert(pg==pgt,'pg must be pgt');
    if(pgt == 1)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'h3ddirectcp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pot] = fmm3d(mex_id_, nd, zk, sources, charges, ns, stmp, ntest, pot, thresh, 1, 1, 3, ns, nd, ns, 1, 3, ntest, 1, nd, ntest, 1);
        mex_id_ = 'h3ddirectcp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, charges, ns, ttmp, ntest, pottarg, thresh, 1, 1, 3, ns, nd, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'h3ddirectdp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pot] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, stmp, ntest, pot, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, 1);
        mex_id_ = 'h3ddirectdp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, ttmp, ntest, pottarg, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'h3ddirectcdp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pot] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, stmp, ntest, pot, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, 1);
        mex_id_ = 'h3ddirectcdp(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], i double[x])';
[pottarg] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, ttmp, ntest, pottarg, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, 1);
      end
      U.pot = pot;
      U.pottarg = pottarg;
    end
    if(pgt == 2)
      if(ifcharge==1 && ifdipole == 0)
        mex_id_ = 'h3ddirectcg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pot, grad] = fmm3d(mex_id_, nd, zk, sources, charges, ns, stmp, ntest, pot, grad, thresh, 1, 1, 3, ns, nd, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
        mex_id_ = 'h3ddirectcg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, charges, ns, ttmp, ntest, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      if(ifcharge==0 && ifdipole == 1)
        mex_id_ = 'h3ddirectdg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pot, grad] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, stmp, ntest, pot, grad, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
        mex_id_ = 'h3ddirectdg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, dipoles, ns, ttmp, ntest, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      if(ifcharge==1 && ifdipole == 1)
        mex_id_ = 'h3ddirectcdg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pot, grad] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, stmp, ntest, pot, grad, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
        mex_id_ = 'h3ddirectcdg(i int[x], i dcomplex[x], i double[xx], i dcomplex[xx], i dcomplex[xx], i int[x], i double[xx], i int[x], io dcomplex[xx], io dcomplex[xx], i double[x])';
[pottarg, gradtarg] = fmm3d(mex_id_, nd, zk, sources, charges, dipoles, ns, ttmp, ntest, pottarg, gradtarg, thresh, 1, 1, 3, ns, nd, ns, nd3, ns, 1, 3, ntest, 1, nd, ntest, nd3, ntest, 1);
      end
      U.pot = pot;
      U.grad = squeeze(reshape(grad,[nd,3,ntest]));
      U.pottarg = pottarg;
      U.gradtarg = squeeze(reshape(gradtarg,[nd,3,ntest]));
    end
  end
end

