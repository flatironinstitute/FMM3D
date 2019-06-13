%
%
disp("This is an example matlab driver");
disp("On output the code prints the pot+grad in example1, and computes the error at 10 source locations");
disp(" ");
disp(" ");
%%   source to source, charges only, pot+grad example
%

clear srcinfo ns U1 U2 eps pg;
disp("Example 1: source to source, charge, pot+grad");
disp(" ");
disp(" ");

ns = 300000000
srcinfo.sources = rand(3,ns);
srcinfo.charges = rand(1,ns);

eps = 1e-5;
pg = 2;
U1 = lfmm3d(eps,srcinfo,pg);

ntest = 10;
targ = srcinfo.sources(:,1:ntest);
U2 = l3ddir(srcinfo,targ,pg);

U2.pot = U2.pottarg;
U2.grad = U2.gradtarg;
err = norm(U1.pot(1:ntest)-U2.pot)^2+norm(U1.grad(:,1:ntest)-U2.grad)^2;
ra = norm(U2.pot)^2 + norm(U2.grad)^2;
errs = sqrt(err/ra)
exit;
