clear srcinfo

ns = 4000;
srcinfo.sources = rand(3,ns);
srcinfo.charges = rand(1,ns)+1i*rand(1,ns);

nt = 3999;
targ = rand(3,nt);  

eps = 1e-5;
ntests = 36;
ipass = zeros(ntests,1);
errs = zeros(ntests,1);
zk = complex(1.1);
ntest = 10;

stmp = srcinfo.sources(:,1:ntest);
ttmp = targ(:,1:ntest);


% Test sources to sources, charge, pot
pg = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
errs(1) = norm(U1.pot(1:ntest)-U2.pot)/norm(U2.pot);
assert(errs(1)<eps,'Failed source to source, charge, pot test');
ipass(1) = 1;

% Test sources to sources, charge, grad
pg = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(2) = sqrt(err/ra);
assert(errs(2)<eps,'Failed source to source, charge, grad test');
ipass(2) = 1;

% Test sources to targets, charge, pot
pg = 0;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
errs(3) = norm(U1.pottarg(1:ntest)-U2.pottarg)/norm(U2.pottarg);
assert(errs(3)<eps,'Failed source to target, charge, pot test');
ipass(3) = 1;

% Test sources to targets, charge, grad
pg = 0;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(4) = sqrt(err/ra);
assert(errs(4)<eps,'Failed source to target, charge, grad test');
ipass(4) = 1;

% Test sources to sources+targets, charge, pot
pg = 1;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2 + norm(U1.pottarg(1:ntest)-U2.pottarg)^2;
ra = norm(U2.pot)^2 + norm(U2.pottarg)^2;
errs(5) = sqrt(err/ra);
assert(errs(5)<eps,'Failed source to source+target, charge, pot test');
ipass(5) = 1;

% Test sources to sources+targets, charge, grad
pg = 2;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp= U2.pottarg;
gradtmp = U2.gradtarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(6) = sqrt(err/ra);
assert(errs(6)<eps,'Failed source to source+target, charge, grad test');
ipass(6) = 1;

%%%%
% Testing for dipoles only
%%%%

srcinfo = rmfield(srcinfo,'charges');
srcinfo.dipoles = rand(3,ns) + 1i*rand(3,ns);

% Test sources to sources, dipole, pot
pg = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
errs(7) = norm(U1.pot(1:ntest)-U2.pot)/norm(U2.pot);
assert(errs(7)<eps,'Failed source to source, dipole, pot test');
ipass(7) = 1;

% Test sources to sources, dipole, grad
pg = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
errs(8) = sqrt(err/ra);
assert(errs(8)<eps,'Failed source to source, dipole, grad test');
ipass(8) = 1;

% Test sources to targets, dipole, pot
pg = 0;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
errs(9) = norm(U1.pottarg(1:ntest)-U2.pottarg)/norm(U2.pottarg);
assert(errs(9)<eps,'Failed source to target, dipole, pot test');
ipass(9) = 1;

% Test sources to target, dipole, grad
pg = 0;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(10) = sqrt(err/ra);
assert(errs(10)<eps,'Failed source to target, dipole, grad test');
ipass(10) = 1;

% Test sources to sources+targets, dipole, pot
pg = 1;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2 + norm(U1.pottarg(1:ntest)-U2.pottarg)^2;
ra = norm(U2.pot)^2 + norm(U2.pottarg)^2;
errs(11) = sqrt(err/ra);
assert(errs(11)<eps,'Failed source to source+target, dipole, pot test');
ipass(11) = 1;

% test sources to sources+target, dipole, grad
pg = 2;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(12) = sqrt(err/ra);
assert(errs(12)<eps,'failed source to source+target, dipole, grad test');
ipass(12) = 1;

%%%%
% Testing for charges + dipoles
%%%%

srcinfo.charges = rand(1,ns) + 1i*rand(1,ns);
% Test sources to sources, charge+dipole, pot
pg = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
errs(13) = norm(U1.pot(1:ntest)-U2.pot)/norm(U2.pot);
assert(errs(13)<eps,'Failed source to source, charge+dipole, pot test');
ipass(13) = 1;

% Test sources to sources, charge+dipole, grad
pg = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
errs(14) = sqrt(err/ra);
assert(errs(14)<eps,'Failed source to source, charge+dipole, grad test');
ipass(14) = 1;


% Test sources to target, charge+dipole, pot
pg = 0;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
errs(15) = norm(U1.pottarg(1:ntest)-U2.pottarg)/norm(U2.pottarg);
assert(errs(15)<eps,'Failed source to target, charge+dipole, pot test');
ipass(15) = 1;

% Test sources to target, charge+dipole, grad
pg = 0;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(16) = sqrt(err/ra);
assert(errs(16)<eps,'Failed source to target, charge+dipole, grad test');
ipass(16) = 1;

% Test sources to source+target, charge+dipole, pot
pg = 1;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2 + norm(U1.pottarg(1:ntest)-U2.pottarg)^2;
ra = norm(U2.pot)^2 + norm(U2.pottarg)^2;
errs(17) = sqrt(err/ra);
assert(errs(17)<eps,'Failed source to source+target, charge+dipole, pot test');
ipass(17) = 1;

% Test sources to sources+target, charge+dipole, grad
pg = 2;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(18) = sqrt(err/ra);
assert(errs(18)<eps,'Failed source to source+target, charge+dipole, grad test');
ipass(18) = 1;


%%%
%   Testing the vector routines
%%%%

srcinfo = rmfield(srcinfo,'charges');
srcinfo = rmfield(srcinfo,'dipoles');
nd = 2;
srcinfo.nd = nd;
srcinfo.charges = rand(nd,ns)+1i*rand(nd,ns);


% Test sources to sources, charge, pot
pg = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
errs(19) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(19)<eps,'Failed source to source, charge, pot vec test');
ipass(19) = 1;

% Test sources to sources, charge, grad
pg = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
tmpg = U1.grad(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.grad(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(20) = sqrt(err/ra);
assert(errs(20)<eps,'Failed source to source, charge, grad vec test');
ipass(20) = 1;

% Test sources to target, charge, pot
pg = 0;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
errs(21) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(21)<eps,'Failed source to target, charge, pot vec test');
ipass(21) = 1;

% Test sources to target, charge, grad
pg = 0;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(22) = sqrt(err/ra);
assert(errs(22)<eps,'Failed source to source, charge, grad vec test');
ipass(22) = 1;

% Test sources to source+target, charge, pot
pg = 1;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2;
ra = norm(p2)^2 + norm(tp2)^2;
errs(23) = sqrt(err/ra); 
assert(errs(23)<eps,'Failed source to source+target, charge, pot vec test');
ipass(23) = 1;

% Test sources to sources+target, charge, grad
pg = 2;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
g = U1.grad(:,:,1:ntest);
g = g(:);
g2 = U2.grad(:);
tg = U1.gradtarg(:,:,1:ntest);
tg = tg(:);
tg2 = U2.gradtarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2 + norm(g-g2)^2 + norm(tg-tg2)^2;
ra = norm(p2)^2 + norm(tp2)^2 + norm(g2)^2 + norm(tg2)^2;
errs(24) = sqrt(err/ra);
assert(errs(24)<eps,'Failed source to source+target, charge, grad vec test');
ipass(24) = 1;

%%%%
% Testing for dipoles only
%%%%

srcinfo = rmfield(srcinfo,'charges');
srcinfo.dipoles = rand(nd,3,ns) + 1i*rand(nd,3,ns);

% Test sources to sources, dipole, pot
pg = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
errs(25) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(25)<eps,'Failed source to source, dipole, pot vec test');
ipass(25) = 1;

% Test sources to sources, dipole, grad
pg = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
tmpg = U1.grad(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.grad(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(26) = sqrt(err/ra);
assert(errs(26)<eps,'Failed source to source, dipole, grad vec test');
ipass(26) = 1;


% Test sources to target, dipole, pot
pg = 0;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
errs(27) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(27)<eps,'Failed source to targ, dipole, pot vec test');
ipass(27) = 1;

% Test sources to target, dipole, grad
pg = 0; 
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(28) = sqrt(err/ra);
assert(errs(28)<eps,'Failed source to target, dipole, grad vec test');
ipass(28) = 1;

% Test sources to source+target, dipole, pot
pg = 1;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2;
ra = norm(p2)^2 + norm(tp2)^2;
errs(29) = sqrt(err/ra); 
assert(errs(29)<eps,'Failed source to source+targ, dipole, pot vec test');
ipass(29) = 1;

% Test sources to sources+target, dipole, grad
pg = 2; 
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
g = U1.grad(:,:,1:ntest);
g = g(:);
g2 = U2.grad(:);
tg = U1.gradtarg(:,:,1:ntest);
tg = tg(:);
tg2 = U2.gradtarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2 + norm(g-g2)^2 + norm(tg-tg2)^2;
ra = norm(p2)^2 + norm(tp2)^2 + norm(g2)^2 + norm(tg2)^2;
errs(30) = sqrt(err/ra);
assert(errs(30)<eps,'Failed source to source+target, dipole, grad vec test');
ipass(30) = 1;

%%%%
% Testing for charges + dipoles
%%%%

srcinfo.charges = rand(nd,ns) + 1i*rand(nd,ns);

% Test sources to sources, charge+dipole, pot
pg = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
errs(31) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(31)<eps,'Failed source to source, charge+dipole, pot vec test');
ipass(31) = 1;

% Test sources to sources, charge+dipole, grad
pg = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg);
U2 = h3ddir(zk,srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
tmpg = U1.grad(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.grad(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(32) = sqrt(err/ra);
assert(errs(32)<eps,'Failed source to source, charge+dipole, grad vec test');
ipass(32) = 1;

% Test sources to target, charge+dipole, pot
pg = 0;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
errs(33) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(33)<eps,'Failed source to target, charge+dipole, pot vec test');
ipass(33) = 1;

% Test sources to target, charge+dipole, grad
pg = 0;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(34) = sqrt(err/ra);
assert(errs(34)<eps,'Failed source to target, charge+dipole, grad vec test');
ipass(34) = 1;

% Test sources to source+target, charge+dipole, pot
pg = 1;
pgt = 1;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2;
ra = norm(p2)^2 + norm(tp2)^2;
errs(35) = sqrt(err/ra); 
assert(errs(35)<eps,'Failed source to target, charge+dipole, pot vec test');
ipass(35) = 1;

% Test sources to sources+target, charge+dipole, grad
pg = 2;
pgt = 2;
U1 = hfmm3d(eps,zk,srcinfo,pg,targ,pgt);
U2 = h3ddir(zk,srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
U2 = h3ddir(zk,srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
g = U1.grad(:,:,1:ntest);
g = g(:);
g2 = U2.grad(:);
tg = U1.gradtarg(:,:,1:ntest);
tg = tg(:);
tg2 = U2.gradtarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2 + norm(g-g2)^2 + norm(tg-tg2)^2;
ra = norm(p2)^2 + norm(tp2)^2 + norm(g2)^2 + norm(tg2)^2;
errs(36) = sqrt(err/ra);
assert(errs(36)<eps,'Failed source to target, charge+dipole, grad vec test');
ipass(36) = 1;

isum = sum(ipass);
fprintf("Successfully cleared %d out of 36 tests in Helmholtz testing suite\n",isum);


