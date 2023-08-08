clear srcinfo

ns = 4000;
srcinfo.sources = rand(3,ns);
srcinfo.charges = rand(1,ns);

nt = 3999;
targ = rand(3,nt);

eps = 1e-5;
ntests = 54;
ipass = zeros(ntests,1);
errs = zeros(ntests,1);
ntest = 10;

stmp = srcinfo.sources(:,1:ntest);
ttmp = targ(:,1:ntest);

itest = 1;


% Test sources to sources, charge, pot
pg = 1;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
errs(itest) = norm(U1.pot(1:ntest)-U2.pot)/norm(U2.pot);
assert(errs(itest)<eps,'Failed source to source, charge, pot test');
ipass(itest) = 1;


itest = itest + 1;
% Test sources to sources, charge, grad
pg = 2;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge, grad test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources, charge, hess
pg = 3;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
U2.hess = U2.hesstarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err + norm(U1.hess(:,1:ntest)-U2.hess)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg(:))^2;
ra = ra + norm(U2.hess(:))^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge, hess test');
ipass(itest) = 1;



itest = itest + 1;
% Test sources to targets, charge, pot
pg = 0;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
errs(itest) = norm(U1.pottarg(1:ntest)-U2.pottarg)/norm(U2.pottarg);
assert(errs(itest)<eps,'Failed source to target, charge, pot test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to targets, charge, grad
pg = 0;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to target, charge, grad test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to targets, charge, hess
pg = 0;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
err = err + norm(U1.hesstarg(:,1:ntest)-U2.hesstarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
ra = ra + norm(U2.hesstarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to target, charge, hess test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+targets, charge, pot
pg = 1;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2 + norm(U1.pottarg(1:ntest)-U2.pottarg)^2;
ra = norm(U2.pot)^2 + norm(U2.pottarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge, pot test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources+targets, charge, grad
pg = 2;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp= U2.pottarg;
gradtmp = U2.gradtarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge, grad test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+targets, charge, hess
pg = 3;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp= U2.pottarg;
gradtmp = U2.gradtarg;
hesstmp = U2.hesstarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
U2.hess = hesstmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
err = err+norm(U1.hess(:,1:ntest)-U2.hess)^2;
err = err+norm(U1.hesstarg(:,1:ntest)-U2.hesstarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
ra = ra+norm(U2.hess)^2 + norm(U2.hesstarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge, hess test');
ipass(itest) = 1;




%%%%
% Testing for dipoles only
%%%%

srcinfo = rmfield(srcinfo,'charges');
srcinfo.dipoles = rand(3,ns); 

itest = itest+1;


% Test sources to sources, dipole, pot
pg = 1;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
errs(itest) = norm(U1.pot(1:ntest)-U2.pot)/norm(U2.pot);
assert(errs(itest)<eps,'Failed source to source, dipole, pot test');
ipass(itest) = 1;


itest = itest + 1;
% Test sources to sources, dipole, grad
pg = 2;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, dipole, grad test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources, dipole, hess
pg = 3;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
U2.hess = U2.hesstarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err + norm(U1.hess(:,1:ntest)-U2.hess)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg(:))^2;
ra = ra + norm(U2.hess(:))^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, dipole, hess test');
ipass(itest) = 1;



itest = itest + 1;
% Test sources to targets, dipole, pot
pg = 0;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
errs(itest) = norm(U1.pottarg(1:ntest)-U2.pottarg)/norm(U2.pottarg);
assert(errs(itest)<eps,'Failed source to target, dipole, pot test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to targets, dipole, grad
pg = 0;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to target, dipole, grad test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to targets, dipole, hess
pg = 0;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
err = err + norm(U1.hesstarg(:,1:ntest)-U2.hesstarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
ra = ra + norm(U2.hesstarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to target, dipole, hess test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+targets, dipole, pot
pg = 1;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2 + norm(U1.pottarg(1:ntest)-U2.pottarg)^2;
ra = norm(U2.pot)^2 + norm(U2.pottarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, dipole, pot test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources+targets, dipole, grad
pg = 2;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp= U2.pottarg;
gradtmp = U2.gradtarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, dipole, grad test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+targets, dipole, hess
pg = 3;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp= U2.pottarg;
gradtmp = U2.gradtarg;
hesstmp = U2.hesstarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
U2.hess = hesstmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
err = err+norm(U1.hess(:,1:ntest)-U2.hess)^2;
err = err+norm(U1.hesstarg(:,1:ntest)-U2.hesstarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
ra = ra+norm(U2.hess)^2 + norm(U2.hesstarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, dipole, hess test');
ipass(itest) = 1;





%%%%
% Testing for charges + dipoles
%%%%

srcinfo.charges = rand(1,ns);


itest = itest+1;


% Test sources to sources, charge+dipole, pot
pg = 1;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
errs(itest) = norm(U1.pot(1:ntest)-U2.pot)/norm(U2.pot);
assert(errs(itest)<eps,'Failed source to source, charge+dipole, pot test');
ipass(itest) = 1;


itest = itest + 1;
% Test sources to sources, charge+dipole, grad
pg = 2;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge+dipole, grad test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources, charge+dipole, hess
pg = 3;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
U2.hess = U2.hesstarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err + norm(U1.hess(:,1:ntest)-U2.hess)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg(:))^2;
ra = ra + norm(U2.hess(:))^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge+dipole, hess test');
ipass(itest) = 1;



itest = itest + 1;
% Test sources to targets, charge+dipole, pot
pg = 0;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
errs(itest) = norm(U1.pottarg(1:ntest)-U2.pottarg)/norm(U2.pottarg);
assert(errs(itest)<eps,'Failed source to target, charge+dipole, pot test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to targets, charge+dipole, grad
pg = 0;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to target, charge+dipole, grad test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to targets, charge+dipole, hess
pg = 0;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
err = norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
err = err + norm(U1.hesstarg(:,1:ntest)-U2.hesstarg)^2;
ra = norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
ra = ra + norm(U2.hesstarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to target, charge+dipole, hess test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+targets, charge+dipole, pot
pg = 1;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2 + norm(U1.pottarg(1:ntest)-U2.pottarg)^2;
ra = norm(U2.pot)^2 + norm(U2.pottarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge+dipole, pot test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources+targets, charge+dipole, grad
pg = 2;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp= U2.pottarg;
gradtmp = U2.gradtarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge+dipole, grad test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+targets, charge+dipole, hess
pg = 3;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp= U2.pottarg;
gradtmp = U2.gradtarg;
hesstmp = U2.hesstarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
U2.hess = hesstmp;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
err = err+norm(U1.pottarg(1:ntest)-U2.pottarg)^2+norm(U1.gradtarg(:,1:ntest)-U2.gradtarg)^2;
err = err+norm(U1.hess(:,1:ntest)-U2.hess)^2;
err = err+norm(U1.hesstarg(:,1:ntest)-U2.hesstarg)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
ra = ra+norm(U2.pottarg)^2 + norm(U2.gradtarg)^2;
ra = ra+norm(U2.hess)^2 + norm(U2.hesstarg)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge+dipole, hess test');
ipass(itest) = 1;



%%%
%   Testing the vector routines
%%%%

srcinfo = rmfield(srcinfo,'charges');
srcinfo = rmfield(srcinfo,'dipoles');
nd = 2;
srcinfo.nd = nd;
srcinfo.charges = rand(nd,ns);

itest = itest+1;
% Test sources to sources, charge, pot
pg = 1;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
errs(itest) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(itest)<eps,'Failed source to source, charge, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources, charge, grad
pg = 2;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
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
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources, charge, hess
pg = 3;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
U2.hess = U2.hesstarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
tmpg = U1.grad(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.grad(:);
tmph = U1.hess(:,:,1:ntest);
tmph = tmph(:);
tmph2 = U2.hess(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2 + norm(tmph-tmph2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2 + norm(tmph2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge, hess vec test');
ipass(itest) = 1;

itest=itest+1;
% Test sources to target, charge, pot
pg = 0;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
errs(itest) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(itest)<eps,'Failed source to target, charge, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to target, charge, grad
pg = 0;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to target, charge, hess
pg = 0;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
tmph = U1.hesstarg(:,:,1:ntest);
tmph = tmph(:);
tmph2 = U2.hesstarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2+norm(tmph-tmph2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2+norm(tmph2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge, hess vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to source+target, charge, pot
pg = 1;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2;
ra = norm(p2)^2 + norm(tp2)^2;
errs(itest) = sqrt(err/ra); 
assert(errs(itest)<eps,'Failed source to source+target, charge, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources+target, charge, grad
pg = 2;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
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
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+target, charge, hess
pg = 3;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
hesstmp = U2.hesstarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
U2.hess = hesstmp;
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
h = U1.hess(:,:,1:ntest);
h = h(:);
h2 = U2.hess(:);
th = U1.hesstarg(:,:,1:ntest);
th = th(:);
th2 = U2.hesstarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2 + norm(g-g2)^2 + norm(tg-tg2)^2;
err = err + norm(h-h2)^2 + norm(th-th2)^2;
ra = norm(p2)^2 + norm(tp2)^2 + norm(g2)^2 + norm(tg2)^2;
ra = ra + norm(h2)^2+norm(th2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge, hess vec test');
ipass(itest) = 1;


%%%%
% Testing for dipoles only
%%%%

srcinfo = rmfield(srcinfo,'charges');
srcinfo.dipoles = rand(nd,3,ns); 



itest = itest+1;
% Test sources to sources, dipole, pot
pg = 1;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
errs(itest) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(itest)<eps,'Failed source to source, dipole, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources, dipole, grad
pg = 2;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
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
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, dipole, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources, dipole, hess
pg = 3;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
U2.hess = U2.hesstarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
tmpg = U1.grad(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.grad(:);
tmph = U1.hess(:,:,1:ntest);
tmph = tmph(:);
tmph2 = U2.hess(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2 + norm(tmph-tmph2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2 + norm(tmph2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, dipole, hess vec test');
ipass(itest) = 1;

itest=itest+1;
% Test sources to target, dipole, pot
pg = 0;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
errs(itest) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(itest)<eps,'Failed source to target, dipole, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to target, dipole, grad
pg = 0;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, dipole, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to target, dipole, hess
pg = 0;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
tmph = U1.hesstarg(:,:,1:ntest);
tmph = tmph(:);
tmph2 = U2.hesstarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2+norm(tmph-tmph2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2+norm(tmph2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, dipole, hess vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to source+target, dipole, pot
pg = 1;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2;
ra = norm(p2)^2 + norm(tp2)^2;
errs(itest) = sqrt(err/ra); 
assert(errs(itest)<eps,'Failed source to source+target, dipole, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources+target, dipole, grad
pg = 2;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
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
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, dipole, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+target, dipole, hess
pg = 3;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
hesstmp = U2.hesstarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
U2.hess = hesstmp;
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
h = U1.hess(:,:,1:ntest);
h = h(:);
h2 = U2.hess(:);
th = U1.hesstarg(:,:,1:ntest);
th = th(:);
th2 = U2.hesstarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2 + norm(g-g2)^2 + norm(tg-tg2)^2;
err = err + norm(h-h2)^2 + norm(th-th2)^2;
ra = norm(p2)^2 + norm(tp2)^2 + norm(g2)^2 + norm(tg2)^2;
ra = ra + norm(h2)^2+norm(th2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, dipole, hess vec test');
ipass(itest) = 1;



%%%%
% Testing for charges + dipoles
%%%%

srcinfo.charges = rand(nd,ns);


itest = itest+1;
% Test sources to sources, charge+dipole, pot
pg = 1;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
errs(itest) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(itest)<eps,'Failed source to source, charge+dipole, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources, charge+dipole, grad
pg = 2;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
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
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge+dipole, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources, charge+dipole, hess
pg = 3;
U1 = lfmm3d(eps,srcinfo,pg);
U2 = l3ddir(srcinfo,stmp,pg);
U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
U2.hess = U2.hesstarg;
tmpp = U1.pot(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pot(:);
tmpg = U1.grad(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.grad(:);
tmph = U1.hess(:,:,1:ntest);
tmph = tmph(:);
tmph2 = U2.hess(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2 + norm(tmph-tmph2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2 + norm(tmph2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge+dipole, hess vec test');
ipass(itest) = 1;

itest=itest+1;
% Test sources to target, charge+dipole, pot
pg = 0;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
errs(itest) = norm(tmpp-tmpp2)/norm(tmpp2);
assert(errs(itest)<eps,'Failed source to target, charge+dipole, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to target, charge+dipole, grad
pg = 0;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge+dipole, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to target, charge+dipole, hess
pg = 0;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,ttmp,pgt);
tmpp = U1.pottarg(:,1:ntest);
tmpp = tmpp(:);
tmpp2 = U2.pottarg(:);
tmpg = U1.gradtarg(:,:,1:ntest);
tmpg = tmpg(:);
tmpg2 = U2.gradtarg(:);
tmph = U1.hesstarg(:,:,1:ntest);
tmph = tmph(:);
tmph2 = U2.hesstarg(:);
err = norm(tmpp-tmpp2)^2+norm(tmpg-tmpg2)^2+norm(tmph-tmph2)^2;
ra = norm(tmpp2)^2 + norm(tmpg2)^2+norm(tmph2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source, charge+dipole, hess vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to source+target, charge+dipole, pot
pg = 1;
pgt = 1;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
p = U1.pot(:,1:ntest);
p = p(:);
p2 = U2.pot(:);
tp = U1.pottarg(:,1:ntest);
tp = tp(:);
tp2 = U2.pottarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2;
ra = norm(p2)^2 + norm(tp2)^2;
errs(itest) = sqrt(err/ra); 
assert(errs(itest)<eps,'Failed source to source+target, charge+dipole, pot vec test');
ipass(itest) = 1;

itest = itest+1;
% Test sources to sources+target, charge+dipole, grad
pg = 2;
pgt = 2;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
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
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge+dipole, grad vec test');
ipass(itest) = 1;


itest = itest+1;
% Test sources to sources+target, charge+dipole, hess
pg = 3;
pgt = 3;
U1 = lfmm3d(eps,srcinfo,pg,targ,pgt);
U2 = l3ddir(srcinfo,stmp,pg);
pottmp = U2.pottarg;
gradtmp = U2.gradtarg;
hesstmp = U2.hesstarg;
U2 = l3ddir(srcinfo,ttmp,pgt);
U2.pot = pottmp;
U2.grad = gradtmp;
U2.hess = hesstmp;
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
h = U1.hess(:,:,1:ntest);
h = h(:);
h2 = U2.hess(:);
th = U1.hesstarg(:,:,1:ntest);
th = th(:);
th2 = U2.hesstarg(:);
err = norm(p-p2)^2 + norm(tp-tp2)^2 + norm(g-g2)^2 + norm(tg-tg2)^2;
err = err + norm(h-h2)^2 + norm(th-th2)^2;
ra = norm(p2)^2 + norm(tp2)^2 + norm(g2)^2 + norm(tg2)^2;
ra = ra + norm(h2)^2+norm(th2)^2;
errs(itest) = sqrt(err/ra);
assert(errs(itest)<eps,'Failed source to source+target, charge+dipole, hess vec test');
ipass(itest) = 1;


isum = sum(ipass);
fprintf("Successfully cleared %d out of 54 tests in laplace testing suite\n",isum);


