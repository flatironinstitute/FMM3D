clear srcinfo

ns = 4000;
srcinfo.sources = rand(3,ns);
srcinfo.stoklet = rand(3,ns);
srcinfo.strslet = rand(3,ns);
srcinfo.strsvec = rand(3,ns);

nt = 3999;
targ = rand(3,nt);

eps = 1e-5;
ntests = 2;
ipass = zeros(ntests,1);
errs = zeros(ntests,1);
ntest = 10;

stmp = srcinfo.sources(:,1:ntest);
ttmp = targ(:,1:ntest);

itest = 1;


% Test source to source, target, stoklet, strslet, pot, pre, grad
ifppreg = 3;
ifppregtarg = 3;
U1 = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);
U2 = st3ddir(srcinfo,stmp,ifppreg);
U3 = st3ddir(srcinfo,ttmp,ifppregtarg);
s1 = U1.pot(:,1:ntest); s1 = s1(:);
s2 = U2.pottarg; s2 = s2(:);
s3 = U1.pre(1:ntest); s3 = s3(:);
s4 = U2.pretarg; s4 = s4(:);
s5 = U1.grad(:,:,1:ntest); s5 = s5(:);
s6 = U2.gradtarg; s6 = s6(:);
t1 = U1.pottarg(:,1:ntest); t1 = t1(:);
t2 = U3.pottarg; t2 = t2(:);
t3 = U1.pretarg(1:ntest); t3 = t3(:);
t4 = U3.pretarg; t4 = t4(:);
t5 = U1.gradtarg(:,:,1:ntest); t5 = t5(:);
t6 = U3.gradtarg; t6 = t6(:);

err = norm(s1-s2)^2 + norm(s3-s4)^2 + norm(s5-s6)^2 + norm(t1-t2)^2 + norm(t3-t4)^2 + norm(t5-t6)^2;
ra = norm(s2)^2 + norm(s4)^2 + norm(s6)^2 + norm(t2)^2 + norm(t4)^2 + norm(t6)^2;

errs(1) = sqrt(err/ra);
assert(errs(1)<eps,'Failed source to source, target, stoklet, strslet, pot, pre, grad test');
ipass(1) = 1;


% Test nd=2, source to source, target, stoklet, strslet, pot, pre, grad
nd = 2;
srcinfo.nd = nd;
srcinfo.stoklet = rand(nd,3,ns);
srcinfo.strslet = rand(nd,3,ns);
srcinfo.strsvec = rand(nd,3,ns);

U1 = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);
U2 = st3ddir(srcinfo,stmp,ifppreg);
U3 = st3ddir(srcinfo,ttmp,ifppregtarg);
s1 = U1.pot(:,:,1:ntest); s1 = s1(:);
s2 = U2.pottarg; s2 = s2(:);
s3 = U1.pre(:,1:ntest); s3 = s3(:);
s4 = U2.pretarg; s4 = s4(:);
s5 = U1.grad(:,:,:,1:ntest); s5 = s5(:);
s6 = U2.gradtarg; s6 = s6(:);
t1 = U1.pottarg(:,:,1:ntest); t1 = t1(:);
t2 = U3.pottarg; t2 = t2(:);
t3 = U1.pretarg(:,1:ntest); t3 = t3(:);
t4 = U3.pretarg; t4 = t4(:);
t5 = U1.gradtarg(:,:,:,1:ntest); t5 = t5(:);
t6 = U3.gradtarg; t6 = t6(:);

err = norm(s1-s2)^2 + norm(s3-s4)^2 + norm(s5-s6)^2 + norm(t1-t2)^2 + norm(t3-t4)^2 + norm(t5-t6)^2;
ra = norm(s2)^2 + norm(s4)^2 + norm(s6)^2 + norm(t2)^2 + norm(t4)^2 + norm(t6)^2;

errs(2) = sqrt(err/ra);
assert(errs(2)<eps,'Failed nd, source to source, target, stoklet, strslet, pot, pre, grad test');
ipass(2) = 1;

isum = sum(ipass);
fprintf("Successfully cleared %d out of 2 tests in stokes testing suite\n",isum);
