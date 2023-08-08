%
%  Test Helmholtz particle FMMs in R^3
%

zk = complex(1.1);

N = maxNumCompThreads

nsource = 2000;

source = zeros(3,nsource);

idist=1;

if( idist == 1 ),
theta=rand(1,nsource)*pi;
phi=rand(1,nsource)*2*pi;
source(1,:)=.5*cos(phi).*sin(theta);
source(2,:)=.5*sin(phi).*sin(theta);
source(3,:)=.5*cos(theta);
end

if( idist == 2 ),
source(1,:)=rand(1,nsource);
source(2,:)=rand(1,nsource);
source(3,:)=rand(1,nsource);
end

%
%  timings
%

ifcharge=1;
charge = rand(1,nsource);
disp(class(charge))

ifdipole=1;
dipstr = complex(rand(1,nsource));
dipvec = rand(3,nsource);




ifcharge
ifdipole
ifpot = 1
iffld = 1

ntarget = min(10,nsource);
target = source(:,1:nsource);
target(1,:) = target(1,:) + 0.2;
[ndim,ntarget] = size(target);
%%%ntarget = 0;

%%%ntarget = 0;

ntarget
ifpottarg = 1
iffldtarg = 1

if( ntarget == 0 ),
ifpottarg = 0
iffldtarg = 0
end


disp('')
'Helmholtz particle target FMM in R^3'

tic
iprec=2
[U]=hfmm3dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
total_time=toc

'Helmholtz particle direct evaluation in R^3'

tic
[F]=h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
total_time=toc


if( ifpot ), U.pot=U.pot/(4*pi); end
if( iffld ), U.fld=U.fld/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( iffld ), F.fld=F.fld/(4*pi); end

if( ifpot ),
%rms_pot = norm((F.pot),2)/sqrt(nsource)
%rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
rel_error_pot = norm((U.pot - F.pot),2)/norm((F.pot),2)
end

if( iffld ),
%rms_fld = norm(F.fld,2)/sqrt(nsource)
%rms_error_fld = norm(U.fld - F.fld,2)/sqrt(nsource)
rel_error_fld = norm(U.fld - F.fld,2)/norm(F.fld,2)
end
%%%break;

if( ifpottarg ), U.pottarg=U.pottarg/(4*pi); end
if( iffldtarg ), U.fldtarg=U.fldtarg/(4*pi); end

if( ifpottarg ), F.pottarg=F.pottarg/(4*pi); end
if( iffldtarg ), F.fldtarg=F.fldtarg/(4*pi); end

if( ifpottarg ),
%rms_pottarg = norm((F.pottarg),2)/sqrt(nsource)
%rms_error_pottarg = norm((U.pottarg - F.pottarg),2)/sqrt(ntarget)
%norm_pottarg = norm((F.pottarg),2)
rel_error_pottarg = norm((U.pottarg - F.pottarg),2)/norm((F.pottarg),2)
end

if( iffldtarg ),
%rms_fldtarg = norm(F.fldtarg,2)/sqrt(ntarget)
%rms_error_fldtarg = ...
%    norm(U.fldtarg - F.fldtarg,2)/sqrt(ntarget)
rel_error_fldtarg = ...
    norm(U.fldtarg - F.fldtarg,2)/ ...
    norm(F.fldtarg,2)
end
%%%break;
