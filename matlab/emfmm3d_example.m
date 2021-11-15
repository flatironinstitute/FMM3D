
%
%
disp("This is an example matlab driver for the maxwell fmm");
disp("On output the code prints the field in example1 and, field+curl in example2");
disp(" ");
disp(" ");
%%   source to target, magnetic current only, field example
%

clear srcinfo targ ns U eps ifE ifcurlE ifdivE zk;
disp("Example 1: source to source, magnetic current, field only");
disp(" ");
disp(" ");

ns = 20000;
srcinfo.sources = rand(3,ns);
srcinfo.h_current = rand(3,ns)+1i*rand(3,ns);

nt = 4000;
targ = rand(3,nt);

eps = 1e-5;
ifE = 1; ifcurlE = 0; ifdivE=0;
zk = complex(1.1);
U = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE);
disp("pot=");
disp(U.E(:,1:3));

clear srcinfo targ ns U eps ifE ifcurlE ifdivE zk;
%
%%    source to target, electric current + charge, field + curl example
%
disp(" ");
disp(" ");
disp("Example 2: source to target, electric current+charge, multiple densities, field+curl");
disp(" ");
disp(" ");

ns = 20000;
nt = 19990;

nd = 3;
srcinfo.nd = nd;

ifE = 1; ifcurlE = 1; ifdivE = 0;

srcinfo.sources = rand(3,ns);
srcinfo.e_charge = rand(nd,ns)+1i*rand(nd,ns);
srcinfo.e_current = rand(nd,3,ns)+1i*rand(nd,3,ns);

targ = rand(3,nt);
eps = 1e-5;
zk = complex(1.1);

U = emfmm3d(eps,zk,srcinfo,targ,ifE,ifcurlE,ifdivE);

disp("field=");
disp(U.E(:,:,1:3));

disp("curl=");
disp(U.curlE(:,:,1:3));
clear srcinfo ns U eps ifE ifcurlE ifdivE nt targ nd zk;
