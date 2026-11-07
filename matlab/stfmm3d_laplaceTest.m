% Stokeslet velocities: Laplace composition in Fortran versus MATLAB.
% All methods use the Fortran Laplace FMM. STFMM3D packs the four charge
% densities, recenters the points and assembles velocities in Fortran.
% The local MATLAB helper does the same work, with recentering on or off.
% Only Stokeslet velocities are tested. 

clear srcinfo
rng(1)

tol = 1e-8;
source_counts = [4096 32768 131072];
target_counts = source_counts;
nrepeat = 7;
shift_ns = 4096;
shift_nt = 32;
shift_tol = 1e-9;
shift_offsets = [0 1e2 1e4 1e6 1e8 1e10 1e12];

%% 1. Agreement at sources and targets, with one and two force densities
% The two recentered compositions should agree within the FMM tolerance.
ns = 513;
nt = 79;
srcinfo.sources = rand(3,ns)-0.5;
srcinfo.stoklet = randn(3,ns);
targ = rand(3,nt)-0.5;
errs = zeros(4,1);

% Test 1: source velocities, one density
U1 = stfmm3d(tol,srcinfo,1);
U2 = laplace_stokeslet(tol,srcinfo,1,zeros(3,0),0,true);
errs(1) = norm(U1.pot(:)-U2.pot(:))/norm(U1.pot(:));

% Test 2: target velocities, one density
U1 = stfmm3d(tol,srcinfo,0,targ,1);
U2 = laplace_stokeslet(tol,srcinfo,0,targ,1,true);
errs(2) = norm(U1.pottarg(:)-U2.pottarg(:))/norm(U1.pottarg(:));

% Test 3: source and target velocities, one density
U1 = stfmm3d(tol,srcinfo,1,targ,1);
U2 = laplace_stokeslet(tol,srcinfo,1,targ,1,true);
err = norm(U1.pot(:)-U2.pot(:))^2 + norm(U1.pottarg(:)-U2.pottarg(:))^2;
ra = norm(U1.pot(:))^2 + norm(U1.pottarg(:))^2;
errs(3) = sqrt(err/ra);

% Test 4: source and target velocities, two densities
srcinfo.nd = 2;
srcinfo.stoklet = randn(2,3,ns);
U1 = stfmm3d(tol,srcinfo,1,targ,1);
U2 = laplace_stokeslet(tol,srcinfo,1,targ,1,true);
err = norm(U1.pot(:)-U2.pot(:))^2 + norm(U1.pottarg(:)-U2.pottarg(:))^2;
ra = norm(U1.pot(:))^2 + norm(U1.pottarg(:))^2;
errs(4) = sqrt(err/ra);

fprintf('\nAgreement of Fortran and MATLAB compositions, tol=%.3g:\n',tol)
agreement_table = table({'Sources'; 'Targets'; 'Both'; 'Both, two densities'}, ...
    errs,'VariableNames',{'evaluation','relative_difference'});
disp(agreement_table)
% Allow differences between the two FMM approximations.
assert(all(errs < 100*tol),'Fortran and MATLAB compositions disagree.');

%% 2. Timings of complete target-velocity calls
% Columns: Fortran (internally centered), MATLAB centered, MATLAB unshifted.
% Warm up each method, then rotate execution order between repetitions.
% Packing, coordinate shifting, the FMM call and assembly are all timed.
ncase = numel(source_counts);
summary = zeros(ncase,5);
raw_total_seconds = zeros(nrepeat,3,ncase);
for problem = 1:ncase
    ns = source_counts(problem);
    nt = target_counts(problem);
    srcinfo = struct('sources',rand(3,ns)-0.5,'stoklet',randn(3,ns));
    targ = rand(3,nt)-0.5;

    stfmm3d(tol,srcinfo,0,targ,1);
    laplace_stokeslet(tol,srcinfo,0,targ,1,true);
    laplace_stokeslet(tol,srcinfo,0,targ,1,false);

    total_seconds = zeros(nrepeat,3);
    velocities = cell(1,3);
    for repeat = 1:nrepeat
        order = circshift(1:3,repeat-1);
        for method = order
            timer = tic;
            if method == 1
                U = stfmm3d(tol,srcinfo,0,targ,1);
            elseif method == 2
                U = laplace_stokeslet(tol,srcinfo,0,targ,1,true);
            else
                U = laplace_stokeslet(tol,srcinfo,0,targ,1,false);
            end
            total_seconds(repeat,method) = toc(timer);
            velocities{method} = U.pottarg;
        end
    end
    raw_total_seconds(:,:,problem) = total_seconds;
    summary(problem,:) = [ns nt median(total_seconds,1)];

    % Near the origin, both MATLAB paths should agree with Fortran.
    for method = 2:3
        relative_error = norm(velocities{method}(:)-velocities{1}(:)) ...
            /norm(velocities{1}(:));
        assert(relative_error < 100*tol,'Timed methods disagree.');
    end
end

fprintf('\nComplete-call timings in seconds, tol=%.3g; median of %d repeats:\n',tol,nrepeat)
timing_table = array2table(summary,'VariableNames', ...
    {'sources','targets','fortran_s','matlab_centered_s','matlab_unshifted_s'});
timing_table.matlab_over_fortran = summary(:,4)./summary(:,3);
disp(timing_table)
disp('matlab_over_fortran > 1 means the Fortran composition was faster.')

%% 3. Accuracy when the same cloud is moved far from the origin
% Translation by (a,-a,a) leaves the exact velocity unchanged. Use grid
% coordinates so even the largest translations preserve all separations.
% Direct sums on the original cloud give an independent reference.
rng(73)
sources = round(4096*(rand(3,shift_ns)-0.5))/4096;
targ = round(4096*(rand(3,shift_nt)-0.5))/4096;
forces = randn(3,shift_ns);
srcinfo = struct('sources',sources,'stoklet',forces);
Udirect = st3ddir(srcinfo,targ,1);
reference = Udirect.pottarg;
shift_errors = zeros(numel(shift_offsets),3);

for k = 1:numel(shift_offsets)
    offset = shift_offsets(k)*[1;-1;1];
    srcinfo.sources = sources+offset;
    translated_targets = targ+offset;

    U1 = stfmm3d(shift_tol,srcinfo,0,translated_targets,1);
    U2 = laplace_stokeslet(shift_tol,srcinfo,0,translated_targets,1,true);
    U3 = laplace_stokeslet(shift_tol,srcinfo,0,translated_targets,1,false);
    shift_errors(k,1) = norm(U1.pottarg(:)-reference(:))/norm(reference(:));
    shift_errors(k,2) = norm(U2.pottarg(:)-reference(:))/norm(reference(:));
    shift_errors(k,3) = norm(U3.pottarg(:)-reference(:))/norm(reference(:));
end

fprintf('\nTranslation accuracy: %d sources, %d targets, tol=%.3g\n', ...
    shift_ns,shift_nt,shift_tol)
disp('Relative errors against direct sums; offset a means translation by (a,-a,a).')
translation_error_table = array2table([shift_offsets(:),shift_errors], ...
    'VariableNames',{'offset_a','fortran_centered','matlab_centered','laplace_unshifted'});
disp(translation_error_table)
% Only recentered methods must meet the tolerance. The unshifted column
% demonstrates cancellation errors growing with the coordinate offset.
assert(all(all(shift_errors(:,1:2) < 100*shift_tol)), ...
    'Recentering did not preserve the velocity accuracy.');
fprintf('Agreement, timing and recentering checks passed.\n');

function U = laplace_stokeslet(tol,srcinfo,ifppreg,targ,ifppregtarg,recenter)
% Test helper for nonempty Stokeslet data; flags 0/1 select velocity outputs.
% K(x,y) = 1/(4*pi*|x-y|) is the Laplace charge kernel. Form four potentials:
%   psi_j = sum (f_j/2)*K, j=1,2,3;   psi_4 = sum ((y dot f)/2)*K.
% Then u_i = psi_i - sum_j x_j*d_i(psi_j) + d_i(psi_4).
% The factors 1/2 give the Stokeslet normalization 1/(8*pi).
% With recentering, x and y are measured from the common bounding-box
% midpoint. Displacements are unchanged, but the intermediate terms shrink.
sources = srcinfo.sources;
ns = size(sources,2);
nd = 1;
if isfield(srcinfo,'nd'), nd = srcinfo.nd; end
forces = reshape(srcinfo.stoklet,nd,3,ns);
if ifppregtarg == 0, targ = zeros(3,0); end

if recenter
    lower = min(sources,[],2);
    upper = max(sources,[],2);
    if ~isempty(targ)
        lower = min(lower,min(targ,[],2));
        upper = max(upper,max(targ,[],2));
    end
    center = lower/2+upper/2;
    sources = sources-center;
    targ = targ-center;
end

% Density index varies fastest, followed by the four scalar potentials.
lapinfo.sources = sources;
lapinfo.nd = 4*nd;
moment = sum(reshape(sources,1,3,ns).*forces,2);
lapinfo.charges = reshape(cat(2,forces/2,moment/2),4*nd,ns);
% Stokes velocities need Laplace potentials and gradients (Laplace flag 2).
if ifppregtarg == 1
    [L,ier] = lfmm3d(tol,lapinfo,2*ifppreg,targ,2);
else
    [L,ier] = lfmm3d(tol,lapinfo,2*ifppreg);
end
assert(ier == 0,'Laplace FMM failed.');

U = struct('pot',[],'pottarg',[]);
if ifppreg == 1
    U.pot = assemble_velocity(L.pot,L.grad,sources,nd);
end
if ifppregtarg == 1
    U.pottarg = assemble_velocity(L.pottarg,L.gradtarg,targ,nd);
end
end

function u = assemble_velocity(pot,grad,x,nd)
% grad(d,j,i,m) is d_i(psi_j) for density d at evaluation point m.
n = size(x,2);
pot = reshape(pot,nd,4,n);
grad = reshape(grad,nd,12,n);
x1 = reshape(x(1,:),1,1,n);
x2 = reshape(x(2,:),1,1,n);
x3 = reshape(x(3,:),1,1,n);
u = zeros(nd,3,n);
for i = 1:3
    k = 4*(i-1);
    u(:,i,:) = pot(:,i,:) + grad(:,k+4,:) ...
        - x1.*grad(:,k+1,:) - x2.*grad(:,k+2,:) - x3.*grad(:,k+3,:);
end
u = squeeze(u);
end
