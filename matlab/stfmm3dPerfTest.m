disp("Perf and math test matlab driver for FMM3D/Stokes")
disp("Source to target, stokeslet+stresslet (1-vec), vel pot only...");
% Barnett 8/4/23

clear
ns = 100000;
nt = 100000;
eps = 1e-3;
ifstrs = 1;       % include stresslets?
ifppreg = 0;      % no eval at sources
ifppregtarg = 1;  % just vel out

rng(0)
srcinfo.sources = rand(3,ns);
srcinfo.stoklet = rand(3,ns);   % stokeslet strength vectors
if ifstrs
  srcinfo.strslet = rand(3,ns);   % optional stresslet strengths and n-vecs
  srcinfo.strsvec = rand(3,ns);
end
targ = rand(3,nt);         % targs and srcs in same unit cube

tic;
U = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);
t = toc;
fprintf("%d to %d points, eps=%.3g, done in %.3g s (%.3g tot pts/sec)\n",ns,nt,eps,t,(ns+nt)/t)
u = U.pottarg;
i = randi(nt);   % which targ
fprintf("velpot targ(i=%d) =\n",i);
disp(u(:,i))

% check u (velocity potential) vs direct formula
ui = zeros(3,1);
for j=1:ns
  R = targ(:,i) - srcinfo.sources(:,j);  % x-y with std sign (targ-src).
  r = sqrt(sum(R.^2));                   % dist
  f = srcinfo.stoklet(:,j);              % strength
  ui = ui + (1/r)*f + (1/r^3)*R*dot(f,R);
  if ifstrs
    mu = srcinfo.strslet(:,j);   % this stresslet strength
    nu = srcinfo.strsvec(:,j);   % this stresslet source vector (normal)
    ui = ui - (6/r^5)*R*dot(mu,R)*dot(nu,R);  % apply T_{ijk}, sign matches
  end
end
ui = 0.5 * ui;   % FMM3D is 1/4pi off from true 1/8pi prefactor
fprintf("rel err vs direct at ith targ: %.3g\n",norm(u(:,i)-ui)/norm(ui))


