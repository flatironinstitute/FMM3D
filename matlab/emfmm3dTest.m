clear srcinfo

ns = 4000;
srcinfo.sources = rand(3,ns);
srcinfo.e_charge = rand(1,ns)+1i*rand(1,ns);
srcinfo.e_current = rand(3,ns)+1i*rand(3,ns);
srcinfo.h_current = rand(3,ns)+1i*rand(3,ns);

nt = 3999;
targ = rand(3,nt);

eps = 1e-6;
ntests = 6;
ipass = zeros(ntests,1);
errs = zeros(ntests,1);
zk = complex(1.1);
ntest = 10;

ttmp = targ(:,1:ntest);

% nd = 1 tests
ifE = 1; ifcurlE = 0; ifdivE = 0;
U1 = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE);
U2 = em3ddir(zk,srcinfo,ttmp,ifE,ifcurlE,ifdivE);

err = norm(U1.E(:,1:ntest)-U2.E)^2;
ra = norm(U2.E)^2;
errs(1) = sqrt(err/ra);
assert(errs(1)<eps,'Failed: e_charge, e_current, h_current; E test');
ipass(1) = 1;

ifE = 1; ifcurlE = 1; ifdivE = 0;
U1 = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE);
U2 = em3ddir(zk,srcinfo,ttmp,ifE,ifcurlE,ifdivE);

err = norm(U1.E(:,1:ntest)-U2.E)^2 + norm(U1.curlE(:,1:ntest)-U2.curlE)^2;
ra = norm(U2.E)^2 + norm(U2.curlE)^2;
errs(2) = sqrt(err/ra);
assert(errs(2)<eps,'Failed: e_charge, e_current, h_current; E, curlE test');
ipass(2) = 1;

ifE = 1; ifcurlE = 1; ifdivE = 1;
U1 = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE);
U2 = em3ddir(zk,srcinfo,ttmp,ifE,ifcurlE,ifdivE);

err = norm(U1.E(:,1:ntest)-U2.E)^2 + norm(U1.curlE(:,1:ntest)-U2.curlE)^2 + norm(U1.divE(1:ntest)-U2.divE)^2;
ra = norm(U2.E)^2 + norm(U2.curlE)^2 + norm(U2.divE)^2;
errs(3) = sqrt(err/ra);
assert(errs(3)<eps,'Failed: e_charge, e_current, h_current; E, curlE, divE test');
ipass(3) = 1;

% nd = 2 tests
nd = 2;
srcinfo.nd = nd;
srcinfo.e_charge = rand(nd,ns)+1i*rand(nd,ns);
srcinfo.e_current = rand(nd,3,ns)+1i*rand(nd,3,ns);
srcinfo.h_current = rand(nd,3,ns)+1i*rand(nd,3,ns);

ifE = 1; ifcurlE = 0; ifdivE = 0;
U1 = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE);
U2 = em3ddir(zk,srcinfo,ttmp,ifE,ifcurlE,ifdivE);

ve1 = U1.E(:,:,1:ntest); ve1 = ve1(:);
ve2 = U2.E; ve2 = ve2(:);
err = norm(ve1-ve2)^2;
ra = norm(ve2)^2;
errs(4) = sqrt(err/ra);
assert(errs(4)<eps,'Failed: nd, e_charge, e_current, h_current; E test');
ipass(4) = 1;

ifE = 1; ifcurlE = 1; ifdivE = 0;
U1 = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE);
U2 = em3ddir(zk,srcinfo,ttmp,ifE,ifcurlE,ifdivE);

ve1 = U1.E(:,:,1:ntest); ve1 = ve1(:);
ve2 = U2.E; ve2 = ve2(:);
vce1 = U1.curlE(:,:,1:ntest); vce1 = vce1(:);
vce2 = U2.curlE; vce2 = vce2(:);
err = norm(ve1-ve2)^2 + norm(vce1-vce2)^2;
ra = norm(ve2)^2 + norm(vce2)^2;
errs(5) = sqrt(err/ra);
assert(errs(5)<eps,'Failed: nd, e_charge, e_current, h_current; E, curlE test');
ipass(5) = 1;

ifE = 1; ifcurlE = 1; ifdivE = 1;
U1 = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE);
U2 = em3ddir(zk,srcinfo,ttmp,ifE,ifcurlE,ifdivE);

ve1 = U1.E(:,:,1:ntest); ve1 = ve1(:);
ve2 = U2.E; ve2 = ve2(:);
vce1 = U1.curlE(:,:,1:ntest); vce1 = vce1(:);
vce2 = U2.curlE; vce2 = vce2(:);
vde1 = U1.divE(:,1:ntest); vde1 = vde1(:);
vde2 = U2.divE; vde2 = vde2(:);
err = norm(ve1-ve2)^2 + norm(vce1-vce2)^2 + norm(vde1-vde2)^2;
ra = norm(ve2)^2 + norm(vce2)^2 + norm(vde2)^2;
errs(6) = sqrt(err/ra);
assert(errs(6)<eps,'Failed: nd, e_charge, e_current, h_current; E, curlE, divE test');
ipass(6) = 1;

isum = sum(ipass);
fprintf("Successfully cleared %d out of 6 tests in maxwell testing suite\n",isum);
