%
%
disp("This is an example matlab driver");
disp("On output the code prints the pot+grad in example1 and, pot+pottarg in example2");
disp(" ");
disp(" ");
%%   source to source, charges only, pot+grad example
%

clear srcinfo ns U eps pg zk;
disp("Example 1: source to source, charge, pot+grad");
disp(" ");
disp(" ");

ns = 300000000
srcinfo.sources = rand(3,ns);
srcinfo.charges = rand(1,ns)+1i*rand(1,ns);

eps = 1e-5;
pg = 2;
zk = complex(1.1);
U1 = hfmm3d(eps,zk,srcinfo,pg);

ntest = 10;
targ = srcinfo.sources(:,1:ntest);
U2 = h3ddir(zk,srcinfo,targ,pg);

U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
errs = sqrt(err/ra)
exit;
