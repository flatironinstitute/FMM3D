disp("Perf and math test matlab driver for FMM3D/Stokes")
% Barnett 8/4/23. Added traction output test, 4/12/24.

clear
ns = 100000;
nt = 100000;
eps = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Source to target, stokeslet+stresslet (1-vec), vel pot only...");
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
i = randi(nt);   % which targ to check
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
    ui = ui + (6/r^5)*R*dot(mu,R)*dot(nu,R);  % apply T_{ijk}, sign matches
  end
end
ui = 0.125/pi * ui;   % FMM3D is now true 1/8pi prefactor
fprintf("rel err vs direct at ith targ: %.3g\n\n",norm(u(:,i)-ui)/norm(ui))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Stokeslet sources to targets: compute vel pot, traction at oriented targs...");
% just give the changes from the above setup...
ifppregtarg = 3;  % request everything
srcinfo = rmfield(srcinfo,{'strslet','strsvec'});
targnor = rand(3,nt);   % normal vectors for targets (eg, surface normals)

tic;
U = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);
t = toc;
fprintf("%d to %d points with u,p,gradu: done in %.3g s (%.3g tot pts/sec)\n",ns,nt,t,(ns+nt)/t)
p = U.pretarg;         % pressure (1*ntarg)
gradu = U.gradtarg;    % grad vel (3*3*ntarg)
shearstress = gradu + permute(gradu,[2 1 3]);        % gradu + gradu^T
% stress sigma = -pI + mu*(gradu+gradu^T),   then T = sigma.n    We assume mu=1.
T = -(ones(3,1)*p).*targnor + squeeze(sum(shearstress .* reshape(kron(ones(3,1),targnor),3,3,nt),1));   % compute all tractions (painful tensor contraction in matlab)
i = randi(nt);   % which targ to check
%fprintf("targ(i=%d): p=%.15g,   gradu = \n",i,p(i)); disp(gradu(:,:,i))
fprintf("targ(i=%d) traction T = \n",i);
disp(T(:,i))

% check T (target tractions) vs direct formula
Ti = zeros(3,1);
nori = targnor(:,i);
for j=1:ns
  R = targ(:,i) - srcinfo.sources(:,j);  % x-y with std sign (targ-src).
  r = sqrt(sum(R.^2));                   % dist
  f = srcinfo.stoklet(:,j);              % strength
  Ti = Ti - (3/(4*pi))*(1/r^5)*R*dot(nori,R)*dot(f,R);  %  true T_{ijk} n_j f_k
end
%Ti = Ti * (4*pi);   % FMM3D has 1/4pi missing in prefactor
fprintf("rel err vs direct at ith targ: %.3g\n",norm(T(:,i)-Ti)/norm(Ti))

