disp("This is an example matlab driver for the stokes fmm");
disp("On output the code prints the pot+pre+grad in example1 and, pot+pottarg+pre+pretarg in example2");
disp(" ");
disp(" ");
%%   source to source, stokeslet only, pot+pre+grad example
%

clear srcinfo ns U eps pg;
disp("Example 1: source to source, stokeslet, pot+pre+grad");
disp(" ");
disp(" ");

ns = 40000;
srcinfo.sources = rand(3,ns);
srcinfo.stoklet = rand(3,ns);

eps = 1e-5;
ifppreg = 3;
U = stfmm3d(eps,srcinfo,ifppreg);
disp("pot=");
disp(U.pot(:,1:9));
disp("pre=");
disp(U.pre(1:9));
disp("grad=");
disp(U.grad(:,:,1:3));

clear srcinfo ns U eps ifppreg;
%
%%    source to source+target, stokeslet+stresslet, pot + pressure example
%
disp(" ");
disp(" ");
disp("Example 2: source to source+target, stokeslet + stresslet, multiple densities, pot+pressure");
disp(" ");
disp(" ");

ns = 4000;
nt = 3999;

nd = 2;
srcinfo.nd = nd;

ifppreg = 2;
ifppregtarg = 2;

srcinfo.sources = rand(3,ns);
srcinfo.stoklet = rand(nd,3,ns);
srcinfo.strslet = rand(nd,3,ns);
srcinfo.strsvec = rand(nd,3,ns);

targ = rand(3,nt);
eps = 1e-5;

U = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);
disp("pot=");
disp(U.pot(:,:,1:3));
disp("pre=");
disp(U.pre(:,1:3));
disp("pottarg=");
disp(U.pottarg(:,:,1:3));
disp("pretarg=");
disp(U.pretarg(:,1:3));

clear srcinfo ns U eps ifppreg ifppregtarg nt targ nd;
