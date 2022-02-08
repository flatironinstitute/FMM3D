function [U]=h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarg,targ,ifpottarg,iffldtarg)
%HFMM3DPARTDIRECT Helmholtz interactions in R^3, direct evaluation.
%
% Helmholtz FMM in R^3: evaluate all pairwise particle
% interactions (ignoring self-interactions) and interactions with targs.
%
% [U]=H3DPARTDIRECT(ZK,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC);
%
% [U]=H3DPARTDIRECT(ZK,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD);
%
% [U]=H3DPARTDIRECT(ZK,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         Ntarg,targ);
%
% [U]=H3DPARTDIRECT(ZK,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         Ntarg,targ,IFPOTTARG,IFFLDTARG);
%
%
% This subroutine evaluates the Helmholtz potential and field due
% to a collection of charges and dipoles. We use (exp(ikr)/r) for the 
% Green's function, without the (1/4 pi) scaling. 
% Self-interactions are not-included.
%
% Input parameters:
% 
% zk - complex Helmholtz parameter
% nsource - number of sources
% source - real (3,nsource): source locations
% ifcharge - charge computation flag
%
%         0 => do not compute
%         1 => include charge contribution
% 
% charge - complex (nsource): charge strengths 
% ifdipole - dipole computation flag
%
%         0 => do not compute
%         1 => include dipole contributions
% 
% dipole - complex (nsource): dipole strengths
% dipvec - real (3,source): dipole orientation vectors
%
% ifpot - potential computation flag, 1 => compute the potential, otherwise no
% iffld - field computation flag, 1 => compute the field, otherwise no
%
% ntarg - number of targs
% targ - real (3,ntarg): targ locations
%
% ifpottarg - targ potential computation flag, 
%      1 => compute the targ potential, otherwise no
% iffldtarg - targ field computation flag, 
%      1 => compute the targ field, otherwise no
%
% Output parameters: 
%
% U.pot - complex (nsource) - potential at source locations
% U.fld - complex (3,nsource) - field (i.e. -gradient) at source locations
% U.pottarg - complex (ntarg) - potential at targ locations
% U.fldtarg - complex (3,ntarg) - field (i.e. -gradient) at targ locations
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%


if( nargin == 8 ) 
  ifpot = 1;
  iffld = 1;
  ntarg = 0;
  targ = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 10 ) 
  ntarg = 0;
  targ = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 12 ) 
  ifpottarg = 1;
  iffldtarg = 1;
end

ifcharge = double(ifcharge); ifdipole = double(ifdipole);
ifpot = double(ifpot); iffld = double(iffld);
ifpottarg = double(ifpottarg); iffldtarg = double(iffldtarg);

pot=0;
fld=zeros(3,1);
pottarg=0;
fldtarg=zeros(3,1);

if( ifpot == 1 ), pot=complex(zeros(1,nsource)); end;
if( iffld == 1 ), fld=complex(zeros(3,nsource)); end;
if( ifpottarg == 1 ), pottarg=complex(zeros(1,ntarg)); end;
if( iffldtarg == 1 ), fldtarg=complex(zeros(3,ntarg)); end;

ier=0;

mex_id_ = 'h3dpartdirect(i dcomplex[x], i int64_t[x], i double[xx], i int64_t[x], i dcomplex[], i int64_t[x], i dcomplex[], i double[xx], i int64_t[x], io dcomplex[], i int64_t[x], io dcomplex[], i int64_t[x], i double[], i int64_t[x], io dcomplex[], i int64_t[x], io dcomplex[])';
[pot, fld, pottarg, fldtarg] = fmm3d_legacy(mex_id_, zk, nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, iffld, fld, ntarg, targ, ifpottarg, pottarg, iffldtarg, fldtarg, 1, 1, 3, nsource, 1, 1, 3, nsource, 1, 1, 1, 1, 1);


if( ifpot == 1 ), U.pot=pot; end
if( iffld == 1 ), U.fld=fld; end
if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( iffldtarg == 1 ), U.fldtarg=fldtarg; end
U.ier=ier;

end

