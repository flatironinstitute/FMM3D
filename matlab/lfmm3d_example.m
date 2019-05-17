%
%%   source to source, charges only, pot+grad example
%

ns = 2000;
srcinfo.sources = rand(3,ns);
srcinfo.charges = rand(1,ns);

eps = 1e-5;
pg = 2;
U = lfmm3d(eps,srcinfo,pg);


clear all;
%
%%    source to source+target, charges+dipoles, pot example
%
ns = 2000;
nt = 1999;

nd = 5;
srcinfo.nd = nd;

pg = 1;
pgt = 1;

srcinfo.sources = rand(3,ns);
srcinfo.charges = rand(nd,ns);
srcinfo.dipoles = rand(nd,3,ns);

targ = rand(3,nt);
eps = 1e-5;
zk = complex(1.1);

U = lfmm3d(eps,srcinfo,pg,targ,pgt);
