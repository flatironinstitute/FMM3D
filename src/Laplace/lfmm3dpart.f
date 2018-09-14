cc Copyright (C) 2010-2011: Leslie Greengard and
cc and Manas Rachh
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c        This file contains the main FMM routines and some related
c        subroutines for evaluating Laplace potentials and fields
c        due ot point charges and dipoles
c
c        lfmm3dpartstos - Laplace FMM in R^{3}: evaluate all
c                    pairwise particle interactions (ignoring self-
c                    interaction); Complex charges, dipoles, potentials
c                    and fields
c
c        lfmm3dpartstot - Laplace FMM in R^{3}: evaluate all 
c                    pairwise particle interactions between sources
c                    and targets; complex charges, dipoles, potentials,
c                    and fields
c 
c        lfmm3dpartstost - Laplace FMM in R^{3}: evaluate all pairwise
c                    particle interactions (ignoring self-interaction)
c                    + interaction with targets; complex charges,
c                    dipoles, potentials, and fields
c     
c        rfmm3dpartstos - Laplace FMM in R^{3}: evaluate all
c                    pairwise particle interactions (ignoring self-
c                    interaction); real charges, dipoles, potentials
c                    and fields
c
c        rfmm3dpartstot - Laplace FMM in R^{3}: evaluate all 
c                    pairwise particle interactions between sources
c                    and targets; real charges, dipoles, potentials,
c                    and fields
c 
c        rfmm3dpartstost - Laplace FMM in R^{3}: evaluate all pairwise
c                    particle interactions (ignoring self-interaction)
c                    + interaction with targets; real charges,
c                    dipoles, potentials, and fields
c     
c-----------------------------------------------------------


        subroutine rfmm3dpartstos(ier,iprec,nsource,source,ifcharge,
     $     charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld)
c
c         Laplace FMM in R^{3}: evaluate all pairwise particle
c         interactions (ignoring self-interaction). 
c         We use (1/r) for the Green's function, without the 
c         1/(4 \pi) scaling.
c
c         On input, the charges, dipole strengths, 
c         and on output, the potentials and the field strengths
c         are all double precision
c
c         The main FMM routine permits both evaluation at
c         sources and at a collection of targets.
c         This subroutine is used to simply the user interface
c         by setting the number of targets to zero and calling
c         the more general FMM
c
c         See lfmm3dpartstost for explanation of calling sequence
c         arguments
c
        implicit double precision (a-h,o-z)
        dimension source(3,1)
        double precision charge(1),dipstr(1)
        dimension dipvec(3,1)
        double precision pot(1), fld(3,1)

        double complex, allocatable :: ccharge(:),cdipstr(:)
        double complex, allocatable :: cpot(:),cfld(:,:)

        dimension target(3,1)
        double complex pottarg(1)
        double complex fldtarg(3,1)

        if(ifcharge.eq.1) then
           allocate(ccharge(nsource)) 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              ccharge(i) = charge(i)
           enddo
C$OMP END PARALLEL DO           
        endif

        if(ifcharge.ne.1) then
           allocate(ccharge(1)) 
           ccharge(1) = 0
        endif

        if(ifdipole.eq.1) then
           allocate(cdipstr(nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              cdipstr(i) = dipstr(i)
           enddo
C$OMP END PARALLEL DO           
        endif

        if(ifdipole.ne.1) then
            allocate(cdipstr(1))
            cdipstr(1) = 0
        endif

        if(ifpot.eq.1) allocate(cpot(nsource))
        if(ifpot.ne.1) allocate(cpot(1))

        if(iffld.eq.1) allocate(cfld(3,nsource))
        if(iffld.ne.1) allocate(cfld(3,1))

        ntarget = 0
        ifpottarg = 0
        iffldtarg = 0
        

        call lfmm3dpartstost(ier,iprec,nsource,source,ifcharge,ccharge,
     $         ifdipole,cdipstr,dipvec,ifpot,cpot,iffld,cfld,ntarget,
     $         target,ifpottarg,pottarg,iffldtarg,fldtarg)

        if(ifpot.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              pot(i) = real(cpot(i))
           enddo
C$OMP END PARALLEL DO           
        endif

        if(iffld.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              fld(1,i) = real(cfld(1,i))
              fld(2,i) = real(cfld(2,i))
              fld(3,i) = real(cfld(3,i))
           enddo
C$OMP END PARALLEL DO           
        endif

        return
        end
c--------------------------------------------------------

        subroutine rfmm3dpartstot(ier,iprec,nsource,source,ifcharge,
     $     charge,ifdipole,dipstr,dipvec,ntarget,ifpottarg,pottarg,
     $     iffldtarg,fldtarg)
c
c         Laplace FMM in R^{3}: evaluate all pairwise particle
c         interactions at a given collection of targets. 
c         We use (1/r) for the Green's function, without the 
c         1/(4 \pi) scaling.
c
c         On input, the charges, dipole strengths, 
c         and on output, the potentials and the field strengths
c         are all double precision
c
c         The main FMM routine permits both evaluation at
c         sources and at a collection of targets.
c         This subroutine is used to simply the user interface
c         by setting the number of targets to zero and calling
c         the more general FMM
c
c         See lfmm3dpartstost for explanation of calling sequence
c         arguments
c
        implicit double precision (a-h,o-z)
        dimension source(3,1)
        double precision charge(1),dipstr(1)
        dimension dipvec(3,1)
        double complex pot(1), fld(3,1)

        double complex, allocatable :: ccharge(:),cdipstr(:)

        dimension target(3,1)
        double precision pottarg(1)
        double precision fldtarg(3,1)
        double complex, allocatable :: cpottarg(:),cfldtarg(:,:)

        if(ifcharge.eq.1) then
           allocate(ccharge(nsource)) 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              ccharge(i) = charge(i)
           enddo
C$OMP END PARALLEL DO           
        endif

        if(ifcharge.ne.1) then
           allocate(ccharge(1)) 
           ccharge(1) = 0
        endif

        if(ifdipole.eq.1) then
           allocate(cdipstr(nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              cdipstr(i) = dipstr(i)
           enddo
C$OMP END PARALLEL DO           
        endif

        if(ifdipole.ne.1) then
            allocate(cdipstr(1))
            cdipstr(1) = 0
        endif

        if(ifpottarg.eq.1) allocate(cpottarg(ntarget))
        if(ifpottarg.ne.1) allocate(cpottarg(1))

        if(iffldtarg.eq.1) allocate(cfldtarg(3,ntarget))
        if(iffldtarg.ne.1) allocate(cfldtarg(3,1))

        ifpot = 0
        iffld = 0

        call lfmm3dpartstost(ier,iprec,nsource,source,ifcharge,ccharge,
     $         ifdipole,cdipstr,dipvec,ifpot,pot,iffld,fld,ntarget,
     $         target,ifpottarg,cpottarg,iffldtarg,cfldtarg)

        if(ifpottarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,ntarget
              pottarg(i) = real(cpottarg(i))
           enddo
C$OMP END PARALLEL DO           
        endif

        if(iffldtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,ntarget
              fldtarg(1,i) = real(cfldtarg(1,i))
              fldtarg(2,i) = real(cfldtarg(2,i))
              fldtarg(3,i) = real(cfldtarg(3,i))
           enddo
C$OMP END PARALLEL DO           
        endif

        return
        end
c-------------------------------------------------------------        

        subroutine rfmm3dpartstost(ier,iprec,nsource,source,ifcharge,
     $     charge,ifdipole,dipstr,dipvec,ntarget,target,ifpot,pot,
     $     iffld,fld,ifpottarg,pottarg,iffldtarg,fldtarg)
c
c         Laplace FMM in R^{3}: evaluate all pairwise particle
c         interactions (ignoring self-interactions) and interactions
c         with targets.
c
c         We use (1/r) for the Green's function, without the 
c         1/(4 \pi) scaling. 
c
c         On input, the charges, dipole strengths, 
c         and on output, the potentials and the field strengths
c         are all double precision.
c
c
c         The main FMM routine permits both evaluation at
c         sources and at a collection of targets.
c         This subroutine is used to simply the user interface
c         by setting the number of targets to zero and calling
c         the more general FMM
c
c         See lfmm3dpartstost for explanation of calling sequence
c         arguments
c
        implicit double precision (a-h,o-z)
        dimension source(3,1)
        double precision charge(1),dipstr(1)
        dimension dipvec(3,1)
        double precision pot(1), fld(3,1)

        double complex, allocatable :: ccharge(:),cdipstr(:)
        double complex, allocatable :: cpot(:),cfld(:,:)

        dimension target(3,1)
        double precision pottarg(1)
        double precision fldtarg(3,1)
        double complex, allocatable :: cpottarg(:),cfldtarg(:,:)

        if(ifcharge.eq.1) then
           allocate(ccharge(nsource)) 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              ccharge(i) = charge(i)
           enddo
C$OMP END PARALLEL DO           
        endif

        if(ifcharge.ne.1) then
           allocate(ccharge(1)) 
           ccharge(1) = 0
        endif

        if(ifdipole.eq.1) then
           allocate(cdipstr(nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              cdipstr(i) = dipstr(i)
           enddo
C$OMP END PARALLEL DO           
        endif

        if(ifdipole.ne.1) then
            allocate(cdipstr(1))
            cdipstr(1) = 0
        endif

        if(ifpot.eq.1) allocate(cpot(nsource))
        if(ifpot.ne.1) allocate(cpot(1))

        if(iffld.eq.1) allocate(cfld(3,nsource))
        if(iffld.ne.1) allocate(cfld(3,1))

        if(ifpottarg.eq.1) allocate(cpottarg(ntarget))
        if(ifpottarg.ne.1) allocate(cpottarg(1))

        if(iffldtarg.eq.1) allocate(cfldtarg(3,ntarget))
        if(iffldtarg.ne.1) allocate(cfldtarg(3,1))

        call lfmm3dpartstost(ier,iprec,nsource,source,ifcharge,ccharge,
     $         ifdipole,cdipstr,dipvec,ifpot,cpot,iffld,cfld,ntarget,
     $         target,ifpottarg,cpottarg,iffldtarg,cfldtarg)

        if(ifpot.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              pot(i) = real(cpot(i))
           enddo
C$OMP END PARALLEL DO           
        endif

        if(iffld.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,nsource
              fld(1,i) = real(cfld(1,i))
              fld(2,i) = real(cfld(2,i))
              fld(3,i) = real(cfld(3,i))
           enddo
C$OMP END PARALLEL DO           
        endif


        if(ifpottarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,ntarget
              pottarg(i) = real(cpottarg(i))
           enddo
C$OMP END PARALLEL DO           
        endif

        if(iffldtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)           
           do i=1,ntarget
              fldtarg(1,i) = real(cfldtarg(1,i))
              fldtarg(2,i) = real(cfldtarg(2,i))
              fldtarg(3,i) = real(cfldtarg(3,i))
           enddo
C$OMP END PARALLEL DO           
        endif

        return
        end
c-------------------------------------------------------------        
        subroutine lfmm3dpartstos(ier,iprec,nsource,source,ifcharge,
     $     charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld)
c
c         Laplace FMM in R^{3}: evaluate all pairwise particle
c         interactions (ignoring self-interaction). 
c         We use (1/r) for the Green's function, without the 
c         1/(4 \pi) scaling. 
c
c         On input, the charges, dipole strengths, 
c         and on output, the potentials and the field strengths
c         are all double complex.
c
c         The main FMM routine permits both evaluation at
c         sources and at a collection of targets.
c         This subroutine is used to simply the user interface
c         by setting the number of targets to zero and calling
c         the more general FMM
c
c         See lfmm3dpartstost for explanation of calling sequence
c         arguments
c
        implicit double precision (a-h,o-z)
        dimension source(3,1)
        double complex charge(1),dipstr(1)
        dimension dipvec(3,1)
        double complex pot(1), fld(3,1)

        dimension target(3,1)
        double complex pottarg(1)
        double complex fldtarg(3,1)

        ntarget = 0
        ifpottarg = 0
        iffldtarg = 0
        

        call lfmm3dpartstost(ier,iprec,nsource,source,ifcharge,charge,
     $         ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld,ntarget,
     $         target,ifpottarg,pottarg,iffldtarg,fldtarg)

        return
        end
c--------------------------------------------------
        subroutine lfmm3dpartstot(ier,iprec,nsource,source,ifcharge,
     $     charge,ifdipole,dipstr,dipvec,ntarget,target,
     $     ifpottarg,pottarg,iffldtarg,fldtarg)
c
c         Laplace FMM in R^{3}: evaluate all pairwise particle
c         interactions at a given collection of targets. 
c         We use (1/r) for the Green's function, without the 
c         1/(4 \pi) scaling.
c
c         On input, the charges, dipole strengths, 
c         and on output, the potentials and the field strengths
c         are all double complex
c
c         The main FMM routine permits both evaluation at
c         sources and at a collection of targets.
c         This subroutine is used to simply the user interface
c         by setting the number of targets to zero and calling
c         the more general FMM
c
c         See lfmm3dpartstost for explanation of calling sequence
c         arguments
c
        implicit double precision (a-h,o-z)
        dimension source(3,1)
        double complex charge(1),dipstr(1)
        dimension dipvec(3,1)
        double complex pot(1), fld(3,1)

        dimension target(3,1)
        double complex pottarg(1)
        double complex fldtarg(3,1)

        ifpot = 0
        iffld = 0
        

        call lfmm3dpartstost(ier,iprec,nsource,source,ifcharge,charge,
     $         ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld,ntarget,
     $         target,ifpottarg,pottarg,iffldtarg,fldtarg)

        return
        end
c
c-----------------------------------------------
c
       subroutine lfmm3dpartstost(ier,iprec,nsource,source,ifcharge,
     $    charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld,ntarget,
     $    target,ifpottarg,pottarg,iffldtarg,fldtarg)
c
c        Laplace FMM in R^{3}: evaluate all pairwise particle
c        interactions (ignoring self-interactions) and interactions
c        with targets.
c
c        We use (1/r) for the Green's function, without the
c        1/(4 \pi) scaling.
c
c        This is primarily a memory management code.
c        The actual work is carried out in
c        subroutine lfmm3dparttargmain,
c
c        Input parameters:
c
c   iprec  in: integer  
c                 FMM precision flag
c
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c
c   nsource in: integer  
c                number of sources
c
c   source  in: double precision (3,nsource)
c                source(k,j) is the kth component of the jth
c                source locations
c
c   ifcharge  in: integer  
c             charge computation flag
c              ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c 
c   charge    in: double complex (nsource) 
c              charge strengths
c
c   ifdipole   in: integer
c              dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c
c   dipstr   in:double complex (nsource)
c              dipole strengths
c
c   dipvec   in: double precision (3,nsource) 
c              dipole orientation vectors
c
c   ifpot   in: integer
c              flag for evaluating potential at the sources
c              The potential will be evaluated if ifpot = 1
c
c   iffld   in:integer
c                flag for evaluating the gradient at the sources
c                gradient will be evaluated at the targets if
c                iffldtarg = 1
c
c   ntarget  in: integer  
c                 number of targets 
c
c   target  in: double precision (3,ntarget)
c               target(k,j) is the kth component of the jth
c               target location
c
c   ifpottarg   in: integer
c              flag for evaluating potential at the targets
c              The potential will be evaluated if ifpottarg = 1
c
c   iffldtarg   in:integer
c                flag for evaluating the gradient at the targets
c                gradient will be evaluated at the targets if
c                iffldtarg = 1
c
c     OUTPUT parameters:
c   ier        out: integer
c              Error return code
c              ier =0    => normal execution
c              ier = 4   => cannot allocate tree workspace
c              ier = 16  => cannot allocate mpole epxansion
c                           workspace in FMM
c
c   pot:    out: double complex(nsource) 
c               potential at the source locations
c
c   fld:   out: double complex(3,nsource)
c               gradient at the source locations
c
c   pottarg:    out: double complex(ntarget) 
c               potential at the target locations
c
c   fldtarg:   out: double complex(3,ntarget)
c               gradient at the target locations
     
       implicit none
cf2py   intent(out) ier
cf2py   intent(in) iprec
cf2py   intent(in) nsource, source
cf2py   intent(in) ifcharge,charge
cf2py   check(!ifcharge || (shape(charge,0) == nsource))  charge
cf2py   depend(nsource)  charge
cf2py   intent(in) ifdipole,dipvec,dipstr
cf2py   check(!ifdipole || (shape(dipstr,0) == nsource))  dipstr
cf2py   depend(nsource)  dipstr
cf2py   intent(in) ifpot,iffld
cf2py   intent(out) pot,fld
cf2py   intent(in) ifpottarg, iffldtarg
cf2py   intent(in) target
cf2py   intent(in) ntarget
cf2py   check((!ifpottarg && !iffldtarg) || (shape(target,0)==3 && shape(target,1) == ntarget))  target
cf2py   check((!ifpottarg) || (shape(pottarg,0)==ntarget))  pottarg
       integer ier,iprec,ifcharge,ifdipole
       integer ifpot,iffld,ifpottarg,iffldtarg

       integer ntarget,nsource


       double precision source(3,1),target(3,1)
       double complex charge(1),dipstr(1)
       double precision dipvec(3,1)

       double complex pot(1),fld(3,1),pottarg(1),fldtarg(3,1)

c
cc       tree variables
c
       integer ltree,mhung,idivflag,ndiv,isep,nboxes,nbmax,nlevels
       integer nlmax
       integer mnbors,mnlist1,mnlist2,mnlist3,mnlist4
       integer ipointer(32)
       integer, allocatable :: itree(:)
       double precision, allocatable :: treecenters(:,:),boxsize(:)

c
cc       temporary sorted arrays
c
       double precision, allocatable :: sourcesort(:,:),targetsort(:,:)
       double precision, allocatable :: radsrc(:)
       double complex, allocatable :: chargesort(:),dipstrsort(:)
       double precision, allocatable :: dipvecsort(:,:)

       integer, allocatable :: flagsort(:)
       double complex, allocatable :: potsort(:),fldsort(:,:)
       double complex, allocatable :: pottargsort(:),fldtargsort(:,:)
c
cc        temporary fmm arrays
c
       double precision epsfmm
       integer, allocatable :: nterms(:),iaddr(:,:)
       double precision, allocatable :: scales(:)
       double precision, allocatable :: rmlexp(:)

       integer lmptemp,nmax,lmptot
       double complex, allocatable :: mptemp(:),mptemp2(:)

c
cc       temporary variables not main fmm routine but
c        not used in particle code
       double precision expc(3),texpssort(100),scjsort,radexp
       double precision expcsort(3)
       integer ntj,nexpc,nadd

c
cc         other temporary variables
c
        integer i,iert,time1,time2,ifprint,ilev,omp_get_wtime,second

c
cc        figure out tree structure
c
c
cc         set criterion for box subdivision
c 
       if(iprec.eq.1) ndiv = 100
       if(iprec.eq.2) ndiv = 100
       if(iprec.eq.3) ndiv = 100
       if(iprec.eq.4) ndiv = 100


c
cc      set tree flags
c 
       isep = 1
       nlmax = 200
       nlevels = 0
       nboxes = 0
       mhung = 0
       ltree = 0

       nexpc = 0
       nadd = 0
       ntj = 0

       idivflag = 0

       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       nbmax = 0


       allocate(radsrc(nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)       
       do i=1,nsource
          radsrc(i) = 0
       enddo
C$OMP END PARALLEL DO  


       allocate(flagsort(ntarget))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
       do i=1,ntarget
          flagsort(i) = -1
       enddo
C$OMP END PARALLEL DO

       radexp = 0

c
cc     memory management code for contructing level restricted tree
        iert = 0
        call mklraptreemem(iert,source,nsource,radsrc,target,ntarget,
     1        expc,nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,
     2        nlevels,nboxes,mnbors,mnlist1,mnlist2,mnlist3,
     3        mnlist4,mhung,ltree)


        if(iert.ne.0) then
           ier = 4
           call prin2('Error in allocating tree memory, ier=*',ier,1)
           stop
        endif


        allocate(itree(ltree))
        allocate(boxsize(0:nlevels))
        allocate(treecenters(3,nboxes))

c       Call tree code
        call mklraptree(source,nsource,radsrc,target,ntarget,expc,
     1               nexpc,radexp,idivflag,ndiv,isep,mhung,mnbors,
     2               mnlist1,mnlist2,mnlist3,mnlist4,nlevels,
     2               nboxes,treecenters,boxsize,itree,ltree,ipointer)

       call prinf('nboxes=*',nboxes,1)

c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
      ifprint=0
c
c     set fmm tolerance based on iprec flag.
c       
      if( iprec .eq. 1 ) epsfmm=.5d-3
      if( iprec .eq. 2 ) epsfmm=.5d-6
      if( iprec .eq. 3 ) epsfmm=.5d-9
      if( iprec .eq. 4 ) epsfmm=.5d-12
c      
      if(ifprint .eq. 1) call prin2('epsfmm=*',epsfmm,1)

c     Allocate sorted source and target arrays      

      allocate(sourcesort(3,nsource))
      allocate(targetsort(3,ntarget))
      if(ifcharge.eq.1) allocate(chargesort(nsource))

      if(ifdipole.eq.1) then
         allocate(dipstrsort(nsource),dipvecsort(3,nsource))
      endif


      if(ifpot.eq.1) allocate(potsort(nsource))
      if(iffld.eq.1) allocate(fldsort(3,nsource))
      
      if(ifpottarg.eq.1) allocate(pottargsort(ntarget))
      if(iffldtarg.eq.1) allocate(fldtargsort(3,ntarget))


c     scaling factor for multipole and local expansions at all levels
c
      allocate(scales(0:nlevels),nterms(0:nlevels))
      do ilev = 0,nlevels
          scales(ilev) = boxsize(ilev)
      enddo

c
cc      initialize potential and field at source
c       locations
c

      if(ifpot.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nsource
         potsort(i) = 0.0d0
      enddo
C$OMP END PARALLEL DO
      endif

      if(iffld.eq.1) then 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nsource
         fldsort(1,i) = 0.0d0
         fldsort(2,i) = 0.0d0
         fldsort(3,i) = 0.0d0
      enddo
C$OMP END PARALLEL DO
      endif


c
cc       initialize potential and field at target
c        locations
c

      if(ifpottarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarget
         pottargsort(i) = 0.0d0
      enddo
C$OMP END PARALLEL DO
      endif

      if(iffldtarg.eq.1) then 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarget
         fldtargsort(1,i) = 0.0d0
         fldtargsort(2,i) = 0.0d0
         fldtargsort(3,i) = 0.0d0
      enddo
C$OMP END PARALLEL DO
      endif


c     Compute length of expansions at each level      
      nmax = 0
      do i=0,nlevels
         call l3dterms(epsfmm,nterms(i),ier)
         if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))
      lmptemp = (nmax+1)*(2*nmax+1)*2
      allocate(mptemp(lmptemp),mptemp2(lmptemp))

c
cc     reorder sources 
c
      call dreorderf(3,nsource,source,sourcesort,itree(ipointer(5)))
      if(ifcharge.eq.1) call dreorderf(2,nsource,charge,chargesort,
     1                     itree(ipointer(5)))

      if(ifdipole.eq.1) then
         call dreorderf(2,nsource,dipstr,dipstrsort,itree(ipointer(5)))
         call dreorderf(3,nsource,dipvec,dipvecsort,itree(ipointer(5)))
      endif

c
cc      reorder targets
c
      call dreorderf(3,ntarget,target,targetsort,itree(ipointer(6)))
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
      call mpalloc(itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
      if(ifprint.ge. 1) call prinf(' lmptot is *',lmptot,1)


      allocate(rmlexp(lmptot),stat=ier)
      if(ier.ne.0) then
         call prinf('Cannot allocate mpole expansion workspace,
     1              lmptot is *', lmptot,1)
         ier = 16
         return
      endif

c     Memory allocation is complete. 
c     Call main fmm routine


      time1=second()
C$      time1=omp_get_wtime()
      call lfmm3dmain(ier,iprec,
     $   nsource,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipstrsort,dipvecsort,
     $   ntarget,targetsort,nexpc,expcsort,
     $   epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $   itree,ltree,ipointer,isep,ndiv,nlevels,
     $   nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $   scales,treecenters,itree(ipointer(1)),nterms,
     $   ifpot,potsort,iffld,fldsort,
     $   ifpottarg,pottargsort,iffldtarg,fldtargsort,flagsort,ntj,
     $   texpssort,scjsort,nadd)

      time2=second()
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)

c
c     parameter ier from targmain routine is currently 
c     meaningless, reset to 0
      if( ier .ne. 0 ) ier = 0

      if(ifpot.eq.1) call dreorderi(2,nsource,potsort,pot,
     1                 itree(ipointer(5)))
      if(iffld.eq.1) call dreorderi(6,nsource,fldsort,fld,
     1                 itree(ipointer(5)))

      if(ifpottarg.eq.1) call dreorderi(2,ntarget,pottargsort,pottarg,
     1                 itree(ipointer(6)))
      if(iffldtarg.eq.1) call dreorderi(6,ntarget,fldtargsort,fldtarg,
     1                 itree(ipointer(6)))

   
       return
       end

c       
c---------------------------------------------------------------
c
      subroutine lfmm3dmain(ier,iprec,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ntarget,targetsort,nexpc,expcsort,
     $     epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $     itree,ltree,ipointer,isep,ndiv,nlevels, 
     $     nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $     scales,centers,laddr,nterms,
     $     ifpot,pot,iffld,fld,
     $     ifpottarg,pottarg,iffldtarg,fldtarg,flagsort,ntj,
     $     tsort,scjsort,nadd)
      implicit none

      integer ier,iprec
      integer nsource,ntarget,nexpc
      integer ndiv,nlevels

      integer ifcharge,ifdipole
      integer ifpot,iffld
      integer ifpottarg,iffldtarg
      double precision epsfmm

      double precision sourcesort(3,nsource)

      double complex chargesort(1)
      double complex dipstrsort(1)
      double precision dipvecsort(3,1)

      double precision targetsort(3,ntarget)
      integer flagsort(ntarget)

      double complex pot(1)
      double complex fld(3,1)
      double complex pottarg(1)
      double complex fldtarg(3,1)

      integer ntj
      double precision expcsort(3,nexpc)
      double complex tsort(0:ntj,-ntj:ntj,nexpc)
      double precision scjsort(nexpc)

      integer iaddr(2,nboxes), lmptot, lmptemp
      double precision rmlexp(lmptot)
      double precision mptemp(lmptemp)
      double precision mptemp2(lmptemp)
       
      double precision timeinfo(10)
      double precision centers(3,nboxes)

      integer isep, ltree, nadd
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer ipointer(32)
      integer itree(ltree)
      integer nboxes
      double precision scales(0:nlevels)
      double precision boxsize(0:nlevels)

      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer uall(200),dall(200),nall(120),sall(120),eall(72),wall(72)
      integer u1234(36),d5678(36),n1256(24),s3478(24)
      integer e1357(16),w2468(16),n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(20),e3(5),e5(5),e7(5),w2(5),w4(5),w6(5),w8(5)

c     temp variables
      integer i,j,k,l,ii,jj,kk,ll,m
      integer ibox,jbox,ilev,npts
      integer nchild,nlist1,nlist2,nlist3,nlist4

      integer istart,iend,istarts,iends
      integer istartt,iendt,istarte,iende
      integer isstart,isend,jsstart,jsend
      integer jstart,jend

      integer ifprint

      integer ifhesstarg
      double precision d,time1,time2,second,omp_get_wtime
      double complex pottmp,fldtmp(3),hesstmp(3)

c     PW variables
      integer ntmax, nexpmax, nlams, nmax, nthmax, nphmax,nmax2
      parameter (ntmax = 40)
      parameter (nexpmax = 1600)
      double precision, allocatable :: carray(:,:), dc(:,:)
      double precision, allocatable :: cs(:,:),fact(:),rdplus(:,:,:)
      double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, allocatable :: rdmsq3(:,:,:)
  
      double precision rlams(ntmax), whts(ntmax)

      double precision, allocatable :: rlsc(:,:,:)
      integer nfourier(ntmax), nphysical(ntmax)
      integer nexptot, nexptotp
      double complex, allocatable :: xshift(:,:)
      double complex, allocatable :: yshift(:,:)
      double precision, allocatable :: zshift(:,:)

      double complex fexpe(50000), fexpo(50000), fexpback(100000)
      double complex, allocatable :: mexp(:,:,:)
      double complex, allocatable :: mexpf1(:),mexpf2(:)
      double complex, allocatable :: mexpp1(:),mexpp2(:),mexppall(:,:)

      double complex, allocatable :: tmp(:,:)

      double precision sourcetmp(3)
      double complex chargetmp

      integer ix,iy,iz,ictr
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype
      double precision wlege(40000)
      integer nterms_eval(4,0:nlevels)

      integer mnlist1, mnlist2,mnlist3,mnlist4,mnbors
      double complex eye, ztmp
      double precision alphaj
      integer ctr,nn,iptr1,iptr2
      double precision, allocatable :: rscpow(:)
      double precision pi,errtmp
      double complex ima
      data ima/(0.0d0,1.0d0)/


c     Initialize routines for plane wave mp loc translation
 
      if(isep.eq.1) then
         if(iprec.le.1) nlams = 12
         if(iprec.eq.2) nlams = 20
         if(iprec.eq.3) nlams = 29
         if(iprec.eq.4) nlams = 37

         nlams = nlams + nadd
      endif
      if(isep.eq.2) then
         if(iprec.le.1) nlams = 9
         if(iprec.eq.2) nlams = 15
         if(iprec.eq.3) nlams = 22
         if(iprec.eq.4) nlams = 29
      endif

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
      enddo
      allocate(rscpow(0:nmax))
      allocate(carray(4*nmax+1,4*nmax+1))
      allocate(dc(0:4*nmax,0:4*nmax))
      allocate(rdplus(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdminus(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdsq3(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdmsq3(0:nmax,0:nmax,-nmax:nmax))
      allocate(rlsc(0:nmax,0:nmax,nlams))



c     generate rotation matrices and carray
      call rotgen(nterms,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)


c     generate rlams and weights (these are the nodes
c     and weights for the lambda integral)

      if(isep.eq.1) call vwts(rlams,whts,nlams)
      if(isep.eq.2) call lwtsexp3sep2(nlams,rlams,whts,errtmp)


c     generate the number of fourier modes required to represent the
c     moment function in fourier space

      if(isep.eq.1) call numthetahalf(nfourier,nlams)
      if(isep.eq.2) call numthetahalf_isep2(nfourier,nlams)
 
c     generate the number of fourier modes in physical space
c     required for the exponential representation
      if(isep.eq.1) call numthetafour(nphysical,nlams)
      if(isep.eq.2) call numthetasix(nphysical,nlams)

c     Generate powers of lambda for the exponential basis
      call rlscini(rlsc,nlams,rlams,nmax)

c     Compute total number of plane waves
      nexptotp = 0
      nexptot = 0
      nthmax = 0
      nphmax = 0
      do i=1,nlams
         nexptot = nexptot + nfourier(i)
         nexptotp = nexptotp + nphysical(i)
         if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
         if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
      enddo
      allocate(tmp(0:nmax,-nmax:nmax))

      allocate(xshift(-5:5,nexptotp))
      allocate(yshift(-5:5,nexptotp))
      allocate(zshift(5,nexptotp))

      allocate(mexpf1(nexptot),mexpf2(nexptot),mexpp1(nexptotp))
      allocate(mexpp2(nexptotp),mexppall(nexptotp,16))

      allocate(mexp(nexptotp,nboxes,6))

c     Precompute table for shifting exponential coefficients in 
c     physical domain
      call mkexps(rlams,nlams,nphysical,nexptotp,xshift,yshift,zshift)

c     Precompute table of exponentials for mapping from
c     fourier to physical domain
      call mkfexp(nlams,nfourier,nphysical,fexpe,fexpo,fexpback)
      
c
cc    compute array of factorials

     
      nmax2 = 2*nmax
      allocate(fact(0:nmax2),cs(0:nmax,-nmax:nmax))
      
      d = 1
      fact(0) = d
      do i=1,nmax2

      d=d*sqrt(i+0.0d0)
      fact(i) = d

      enddo

      cs(0,0) = 1.0d0
      do l=1,nmax
      do m=0,l

      cs(l,m) = ((-1)**l)/(fact(l-m)*fact(l+m))
      cs(l,-m) = cs(l,m)

      enddo
      enddo

      call prin2('end of generating plane wave info*',i,0)



      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, 
c     and other things if ifprint=2.
c       
        ifprint=1
        ifhesstarg = 0
c
c
c     ... set the expansion coefficients to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
      do i=1,nexpc
         do j = 0,ntj
           do k=-ntj,ntj
              tsort(j,k,i)=0
           enddo
         enddo
      enddo
C$OMP END PARALLEL DO

c       
        do i=1,10
          timeinfo(i)=0
        enddo

c
c       ... set all multipole and local expansions to zero
c

      do ilev = 0,nlevels

         nn = nterms(ilev)
         istart = laddr(1,ilev)
         iend = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT (SHARED)     
C$OMP$PRIVATE(ibox,iptr1,iptr2)
         do ibox = istart,iend
             iptr1 = iaddr(1,ibox)
             iptr2 = iaddr(2,ibox)
            call l3dzero(rmlexp(iptr1),nn)
            call l3dzero(rmlexp(iptr2),nn)
         enddo
C$OMP END PARALLEL DO         
       enddo

c    initialize legendre function evaluation routines
      nlege = 100
      lw7 = 40000
      call ylgndrfwini(nlege,wlege,lw7,lused7)

c    initialize nterms_eval
      do ilev=0,nlevels
         do itype = 1,4
            call l3dterms_eval(itype,epsfmm,
     1           nterms_eval(itype,ilev),ier)
         enddo
      enddo
       
c
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
        time1=second()
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions



      do ilev=2,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nchild)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(ipointer(3)+ibox-1)
            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)
            npts = iend-istart+1

            if(nchild.eq.0.and.npts.gt.0) then
               if(ifcharge.eq.1) then
                  call l3dformmp_add_trunc(ier,scales(ilev),
     1            sourcesort(1,istart),chargesort(istart),npts,
     2            centers(1,ibox),nterms(ilev), nterms_eval(1,ilev),
     3            rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
               if(ifdipole.eq.1) then
                  call l3dformmp_dp_add_trunc(ier,scales(ilev),
     1            sourcesort(1,istart),dipstrsort(istart),
     2            dipvecsort(1,istart),          
     3            npts,centers(1,ibox),
     4            nterms(ilev),nterms_eval(1,ilev),
     5            rmlexp(iaddr(1,ibox)),wlege,nlege)
               endif

            endif
         enddo
C$OMP END PARALLEL DO      
      enddo



      time2=second()
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1



      if(ifprint.ge.1)
     $   call prinf('=== STEP 2 (form lo) ===*',i,0)
      time1=second()
C$    time1=omp_get_wtime()

      do ilev=2,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nlist4 = itree(ipointer(26)+ibox-1)
            do i=1,nlist4
               jbox = itree(ipointer(27)+(ibox-1)*mnlist4+i-1)

c              Form local expansion for all boxes in list3
c              of the current box


               istart = itree(ipointer(10)+jbox-1)
               iend = itree(ipointer(11)+jbox-1)
               npts = iend-istart+1
               if(ifcharge.eq.1.and.npts.gt.0) then
                  call l3dformta_add_trunc(ier,scales(ilev),
     1            sourcesort(1,istart),chargesort(istart),npts,
     2            centers(1,ibox),nterms(ilev),nterms_eval(1,ilev),
     3            rmlexp(iaddr(2,ibox)),wlege,nlege)
               endif
               if(ifdipole.eq.1.and.npts.gt.0) then
                  call l3dformta_dp_add_trunc(ier,scales(ilev),
     1            sourcesort(1,istart),dipstrsort(istart),
     2            dipvecsort(1,istart),npts,
     2            centers(1,ibox),nterms(ilev),nterms_eval(1,ilev),
     3            rmlexp(iaddr(2,ibox)),wlege,nlege)
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo
      time2=second()
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1


c       
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 3 (merge mp) ====*',i,0)
      time1=second()
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,0,-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  istart = itree(ipointer(10)+jbox-1)
                  iend = itree(ipointer(11)+jbox-1)
                  npts = iend-istart+1
                  if(npts.gt.0) then
                     call l3dmpmpquadu_add(scales(ilev+1),
     1               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),scales(ilev),centers(1,ibox),
     3               rmlexp(iaddr(1,ibox)),nterms(ilev),nterms(ilev),
     4               ier)
                  endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      time2=second()
C$    time2=omp_get_wtime()
      timeinfo(3)=time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (mp to loc) ===*',i,0)
c      ... step 4, convert multipole expansions into local
c       expansions

      time1 = second()
C$        time1=omp_get_wtime()

c
cc     zero out mexp
c 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k)
      do i=1,nboxes
         do j=1,nexptotp
            do k=1,6
               mexp(j,i,k) = 0.0d0
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO      


      do ilev=2,nlevels

         rscpow(0) = 1
         rtmp = scales(ilev)**2
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,tmp,mexpf1,mexpf2,mptemp,ictr)
C$OMP$PRIVATE(ii,jj)
         do ibox=laddr(1,ilev),laddr(2,ilev)

            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)

            npts = iend-istart+1

            if(npts.gt.0) then
c            rescale the multipole expansion

                ictr = 0
                do ii=-nterms(ilev),nterms(ilev)
                   do jj=0,nterms(ilev)
                      tmp(jj,ii) = rmlexp(iaddr(1,ibox)+ictr)+
     1                    ima*rmlexp(iaddr(1,ibox)+ictr+1)
                      tmp(jj,ii) = tmp(jj,ii)/scales(ilev)**(2*jj+1)
                      ictr = ictr + 2
                   enddo
                enddo
c
cc                process up down for current box
c
                call mpoletoexp(tmp,nterms(ilev),nlams,nfourier,
     1              nexptot,mexpf1,mexpf2,rlsc)

                call ftophys(mexpf1,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,ibox,1),fexpe,fexpo)

                call ftophys(mexpf2,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,ibox,2),fexpe,fexpo)


c
cc                process north-south for current box
c
                call rotztoy(nterms(ilev),tmp,mptemp,rdminus)
                call mpoletoexp(mptemp,nterms(ilev),nlams,nfourier,
     1              nexptot,mexpf1,mexpf2,rlsc)

                call ftophys(mexpf1,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,ibox,3),fexpe,fexpo)

                call ftophys(mexpf2,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,ibox,4),fexpe,fexpo)

c
cc                process east-west for current box

                call rotztox(nterms(ilev),tmp,mptemp,rdplus)
                call mpoletoexp(mptemp,nterms(ilev),nlams,nfourier,
     1              nexptot,mexpf1,mexpf2,rlsc)

                call ftophys(mexpf1,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,ibox,5),fexpe,fexpo)


                call ftophys(mexpf2,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,ibox,6),fexpe,fexpo)

            endif

         enddo
C$OMP END PARALLEL DO         
c
c
cc         loop over parent boxes and ship plane wave
c          expansions to the first child of parent 
c          boxes. 
c          The codes are now written from a gathering perspective
c
c          so the first child of the parent is the one
c          recieving all the local expansions
c          coming from all the lists
c
c          
c
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nchild)
C$OMP$PRIVATE(mexpf1,mexpf2,mexpp1,mexpp2,mexppall)
C$OMP$PRIVATE(nuall,uall,ndall,dall,nnall,nall,nsall,sall)
C$OMP$PRIVATE(neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678)
C$OMP$PRIVATE(nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468)
C$OMP$PRIVATE(nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,e57)
C$OMP$PRIVATE(nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7)
C$OMP$PRIVATE(nw2,w2,nw4,w4,nw6,w6,nw8,w8)
         do ibox = laddr(1,ilev-1),laddr(2,ilev-1)
            
            istart = itree(ipointer(12)+ibox-1)
            iend = itree(ipointer(13)+ibox-1)
            npts = iend-istart+1
 
            istart = itree(ipointer(14)+ibox-1)
            iend = itree(ipointer(17)+ibox-1)
            npts = npts + iend-istart+1

            if(ifpot.eq.1.or.iffld.eq.1) then
               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = npts + iend-istart+1
            endif


            nchild = itree(ipointer(3)+ibox-1)

            if(npts.gt.0.and.nchild.gt.0) then

               call getpwlistall(ibox,boxsize(ilev),nboxes,
     1         itree(ipointer(18)+ibox-1),itree(ipointer(19)+
     2         mnbors*(ibox-1)),nchild,itree(ipointer(4)),centers,
     3         isep,nuall,uall,ndall,dall,nnall,nall,nsall,sall,neall,
     4         eall,nwall,wall,nu1234,u1234,nd5678,d5678,nn1256,n1256,
     5         ns3478,s3478,ne1357,e1357,nw2468,w2468,nn12,n12,nn56,n56,
     6         ns34,s34,ns78,s78,ne13,e13,ne57,e57,nw24,w24,nw68,w68,
     7         ne1,e1,ne3,e3,ne5,e5,ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)

               call processudexp(ibox,ilev,nboxes,centers,
     1         itree(ipointer(4)),scales(ilev),nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nuall,uall,nu1234,u1234,ndall,dall,nd5678,d5678,
     5         mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1),mexppall(1,2),
     6         mexppall(1,3),mexppall(1,4),xshift,yshift,zshift,
     7         fexpback,rlsc,rscpow)

               call processnsexp(ibox,ilev,nboxes,centers,
     1         itree(ipointer(4)),scales(ilev),nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nnall,nall,nn1256,n1256,nn12,n12,nn56,n56,nsall,sall,
     5         ns3478,s3478,ns34,s34,ns78,s78,
     6         mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1),mexppall(1,2),
     7         mexppall(1,3),mexppall(1,4),mexppall(1,5),mexppall(1,6),
     8         mexppall(1,7),mexppall(1,8),rdplus,xshift,yshift,zshift,
     9         fexpback,rlsc,rscpow)
               
               call processewexp(ibox,ilev,nboxes,centers,
     1         itree(ipointer(4)),scales(ilev),nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         neall,eall,ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,
     5         ne3,e3,ne5,e5,ne7,e7,nwall,wall,
     5         nw2468,w2468,nw24,w24,nw68,w68,
     5         nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6         mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1),mexppall(1,2),
     7         mexppall(1,3),mexppall(1,4),mexppall(1,5),mexppall(1,6),
     8         mexppall(1,7),mexppall(1,8),mexppall(1,9),
     9         mexppall(1,10),mexppall(1,11),mexppall(1,12),
     9         mexppall(1,13),mexppall(1,14),mexppall(1,15),
     9         mexppall(1,16),rdminus,xshift,yshift,zshift,
     9         fexpback,rlsc,rscpow)

            endif

         enddo
C$OMP END PARALLEL DO         
      enddo
      time2 = second()
C$        time2=omp_get_wtime()
      timeinfo(4) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (split loc) ===*',i,0)

      time1 = second()
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,mptemp)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  call l3dzero(mptemp,nterms(ilev+1))

                  call l3dloclocquadu_add(scales(ilev),centers(1,ibox),
     1            rmlexp(iaddr(2,ibox)),nterms(ilev),
     2            scales(ilev+1),centers(1,jbox),mptemp,nterms(ilev+1),
     3            nterms(ilev+1),ier)

                  call l3dadd(mptemp,rmlexp(iaddr(2,jbox)),
     1            nterms(ilev+1))
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo
      time2 = second()
C$        time2=omp_get_wtime()
      timeinfo(5) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== step 6 (mp eval) ===*',i,0)
      time1 = second()
C$        time1=omp_get_wtime()

      do ilev=1,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,j,i,jbox)
C$OMP$PRIVATE(mptemp,pottmp,fldtmp)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nlist3 = itree(ipointer(24)+ibox-1)

c
cc            shift multipole expansion
c             to local expansion at expansion
c             centers (Note: this part is
c             not relevant for particle codes.
c             It is relevant only for QBX codes)
c 

            istart = itree(ipointer(16)+ibox-1)
            iend = itree(ipointer(17)+ibox-1)

            do j=istart,iend
               do i=1,nlist3
                  jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
c
cc                  shift multipole expansion directly from box
c                   center to expansion center
c  
                    call l3dmplocquadu_add(scales(ilev+1),
     $                centers(1,jbox),rmlexp(iaddr(1,jbox)),
     $                nterms(ilev+1),scjsort(j),expcsort(1,j),mptemp,
     $                ntj,ntj,ier)

                    call l3dadd(mptemp,tsort(0,-ntj,j),ntj)
               enddo
            enddo
           
c
cc           evaluate multipole expansion at source locations
c
            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)

            do j=istart,iend
                do i=1,nlist3
                   jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
                   pottmp = 0
                   fldtmp(1) = 0
                   fldtmp(2) = 0
                   fldtmp(3) = 0

                   call l3dmpevalall_trunc(scales(ilev+1),
     $               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     $               nterms(ilev+1),nterms_eval(1,ilev+1),
     $               sourcesort(1,j),npts,ifpot,
     $               pottmp,iffld,fldtmp,wlege,nlege,ier)

                   if(ifpot.eq.1) pot(j) = pot(j)+pottmp
                   if(iffld.eq.1) then
                      fld(1,j) = fld(1,j) + fldtmp(1)
                      fld(2,j) = fld(2,j) + fldtmp(2)
                      fld(3,j) = fld(3,j) + fldtmp(3)
                   endif
                enddo
            enddo
c
cc           evaluate multipole expansion at target locations
c

            istart = itree(ipointer(12)+ibox-1)
            iend = itree(ipointer(13)+ibox-1)

            npts = 1
            do j=istart,iend
                if(flagsort(j).eq.-1) then
                   do i=1,nlist3
                      jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
                      pottmp = 0
                      fldtmp(1) = 0
                      fldtmp(2) = 0
                      fldtmp(3) = 0

                      call l3dmpevalall_trunc(scales(ilev+1),
     $                  centers(1,jbox),rmlexp(iaddr(1,jbox)),
     $                  nterms(ilev+1),nterms_eval(1,ilev+1),
     $                  targetsort(1,j),npts,ifpottarg,
     $                  pottmp,iffldtarg,fldtmp,wlege,nlege,ier)

                      if(ifpottarg.eq.1) pottarg(j) = pottarg(j)+pottmp
                      if(iffldtarg.eq.1) then
                         fldtarg(1,j) = fldtarg(1,j) + fldtmp(1)
                         fldtarg(2,j) = fldtarg(2,j) + fldtmp(2)
                         fldtarg(3,j) = fldtarg(3,j) + fldtmp(3)
                      endif
                   enddo
                endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      time2 = second()
C$        time2=omp_get_wtime()
      timeinfo(6) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== step 7 (eval lo) ===*',i,0)

c     ... step 7, evaluate all local expansions
      time1 = second()
C$        time1=omp_get_wtime()

      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,i,istart,iend,mptemp,npts,pottmp,fldtmp)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild = itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then
c
cc            shift local expansion
c             to local expansion at expansion
c             centers (Note: this part is
c             not relevant for particle codes.
c             It is relevant only for QBX codes)
c 
               istart = itree(ipointer(16)+ibox-1)
               iend = itree(ipointer(17)+ibox-1)
               do i=istart,iend
             
                   call l3dzero(mptemp,ntj)
  
                    call d3tataf(scales(ilev),centers(1,ibox),
     1              rmlexp(iaddr(2,ibox)),nterms(ilev),scjsort(i),
     2              expcsort(1,i),mptemp,ntj,cs,nmax,wlege,nlege)

                    call l3dadd(mptemp,tsort(0,-ntj,i),ntj)
               enddo

c
cc              evaluate local expansion at source locations
c

               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = 1
               do i=istart,iend
                  pottmp = 0.0d0
                  fldtmp(1) = 0.0d0
                  fldtmp(2) = 0.0d0
                  fldtmp(3) = 0.0d0
                  call l3dtaevalall_trunc(scales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),nterms_eval(1,ilev),
     3                sourcesort(1,i),          
     3                npts,ifpot,pottmp,iffld,fldtmp,
     3                wlege,nlege,ier)

                  if(ifpot.eq.1) pot(i) = pot(i)+pottmp
                  if(iffld.eq.1) then 
                     fld(1,i) = fld(1,i)+fldtmp(1)
                     fld(2,i) = fld(2,i)+fldtmp(2)
                     fld(3,i) = fld(3,i)+fldtmp(3)
                  endif
              enddo
c
cc              evaluate local expansion at target locations
c

               istart = itree(ipointer(12)+ibox-1)
               iend = itree(ipointer(13)+ibox-1)
               npts = 1
               do i=istart,iend
                  if(flagsort(i).eq.-1) then
                     pottmp = 0.0d0
                     fldtmp(1) = 0.0d0
                     fldtmp(2) = 0.0d0
                     fldtmp(3) = 0.0d0
                     call l3dtaevalall_trunc(scales(ilev),
     1                   centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                   nterms(ilev),nterms_eval(1,ilev),
     3                   targetsort(1,i),          
     3                   npts,ifpottarg,pottmp,iffldtarg,fldtmp,
     3                   wlege,nlege,ier)

                     if(ifpottarg.eq.1) pottarg(i) = pottarg(i)+pottmp
                     if(iffldtarg.eq.1) then 
                        fldtarg(1,i) = fldtarg(1,i)+fldtmp(1)
                        fldtarg(2,i) = fldtarg(2,i)+fldtmp(2)
                        fldtarg(3,i) = fldtarg(3,i)+fldtmp(3)
                     endif
                 endif
              enddo
            endif
         enddo
C$OMP END PARALLEL DO      
      enddo

      time2 = second()
C$        time2=omp_get_wtime()
      timeinfo(7) = time2 - time1

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
      time1=second()
C$        time1=omp_get_wtime()

      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,istarte,iende,nlist1,i,jbox)
C$OMP$PRIVATE(jstart,jend,istarts,iends)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istarte = itree(ipointer(16)+ibox-1)
            iende = itree(ipointer(17)+ibox-1)

            istarts = itree(ipointer(10)+ibox-1)
            iends = itree(ipointer(11)+ibox-1)

            istartt = itree(ipointer(12)+ibox-1)
            iendt = itree(ipointer(13)+ibox-1)

            nlist1 = itree(ipointer(20)+ibox-1)
            do i =1,nlist1
               jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)

               jstart = itree(ipointer(10)+jbox-1)
               jend = itree(ipointer(11)+jbox-1)

c
cc            directly form local expansion for list1
c             sources at expansion centers
c             (Note: this part is
c             not relevant for particle codes.
c             It is relevant only for QBX codes)
c 
               call lfmm3dexpc_direct(jstart,jend,
     1             istarte,iende,scjsort,sourcesort, 
     2             ifcharge,chargesort,ifdipole,dipstrsort,
     3             dipvecsort,expcsort,tsort,ntj,wlege,nlege)
c
cc               directly evaluate potential at
c                sources due to sources in the near-field
c
               call lfmm3dsrc_direct(jstart,jend,
     1          istarts,iends,sourcesort,ifcharge,chargesort,
     2          ifdipole,dipstrsort,dipvecsort,ifpot,
     3          pot,iffld,fld)

c
cc                directly evaluate potential at targets
c                 due to sources in the near-field

               call lfmm3dtarg_direct(jstart,jend,
     1          istartt,iendt,sourcesort,ifcharge,chargesort,
     2          ifdipole,dipstrsort,dipvecsort,targetsort,ifpottarg,
     3          pottarg,iffldtarg,fldtarg,flagsort)

            enddo   
         enddo
C$OMP END PARALLEL DO      
      enddo

      time2 = second()
C$        time2=omp_get_wtime()

      timeinfo(8) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,8)
      d = 0
      do i = 1,8
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

      return
      end
c------------------------------------------------
      subroutine lfmm3dexpc_direct(istart,iend,jstart,jend,
     $     scj,source,ifcharge,charge,ifdipole,dipstr,
     $     dipvec,expc,texps,ntj,wlege,nlege)
c--------------------------------------------------------------------
c     This subroutine adds the local expansions due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the expansion center array to the existing 
c     local expansions at the corresponding expansion centers.
c
c     INPUT arguments
c-------------------------------------------------------------------
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in the expansion center array at 
c                  which we  wish to compute the expansions
c 
c     jend         in:Integer
c                  Last index in expansion center array at 
c                  which we wish to compute the expansions
c 
c     scjsort      in: double precision(*)
c                  Scale of expansions formed at the expansion centers
c
c     source       in: double precision(3,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: double complex
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: double complex(ns)
c                   dip strengths at the source locations
c
c     dipvec      in: double precision(3,ns)
c                 Dipole orientation vector at the source locations
c
c     expc        in: double precision(3,nexpc)
c                 Expansion center locations
c
c     ntj         in: Integer
c                 Number of terms in expansion
c
c     wlege       in: double precision(0:nlege,0:nlege)
c                 precomputed array of recurrence relation
c                 coeffs for Ynm calculation.
c
c    nlege        in: integer
c                 dimension parameter for wlege
c------------------------------------------------------------
c     OUTPUT
c
c   Updated expansions at the targets
c   texps       out: double complex(0:ntj,-ntj:ntj,expc) 
c                 coeffs for local expansions
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j, nlege
        integer ifcharge,ifdipole,ier
        double precision source(3,*)
        double precision scj(*)
        double precision wlege(*)
        double complex charge(*),dipstr(*)
        double precision dipvec(3,*)
        double precision expc(3,*)

        integer nlevels,ntj
c
        double complex texps(0:ntj,-ntj:ntj,*)
        
c
        ns = iend - istart + 1

        
        do j=jstart,jend
           if(ifcharge.eq.1) then

              call l3dformta_add_trunc(ier,scj(j),
     1        source(1,istart),charge(istart),ns,expc(1,j),
     2        ntj,ntj,texps(0,-ntj,j),wlege,nlege)

           endif

           if(ifdipole.eq.1) then
               call l3dformta_dp_add_trunc(ier,scj(j),
     1         source(1,istart),dipstr(istart),dipvec(1,istart),
     2         ns,expc(1,j),ntj,ntj,texps(0,-ntj,j),wlege,nlege)
           endif        
        enddo
c
        return
        end
c------------------------------------------------------------------     
      subroutine lfmm3dsrc_direct(istart,iend,jstart,jend,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to sources
c     istart to iend in the source array at the sources
c     jstart to jend in the source array to the computed
c     potentials and field strengths
c
c     INPUT arguments
c-------------------------------------------------------------------
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to update the potential and gradients
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to update the potential and gradients
c
c     source       in: double precision(3,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: double complex
c                  Charge at the source locations
c
c     ifdipole    in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr      in: double complex(ns)
c                 dipole strengths at the source locations
c
c     ifpot    in: Integer
c              Flag for computing the potential. The potential
c              will be updated if ifpot = 1
c
c     ifgrad  in: Integer
c             Flag for computing the gradient.
c             The gradient will be updated if ifgrad == 1
c
c------------------------------------------------------------
c     OUTPUT
c
c   Updated potential and field strengths at the sources
c   pot :     out: double complex(*)
c             potential at the targets
c
c   grad:     out: double complex(3,*)
c             field strength at the targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole

        double precision source(3,*)
        double complex charge(*),dipstr(*)
        double precision dipvec(3,*)

        integer ifpot,ifgrad,ifhess
c
        double complex pot(*)
        double complex grad(3,*)

        double complex pottmp,gradtmp(3)

        double complex hesstmp(3)
c
        ifhess = 0
        ns = iend - istart + 1
        do j=jstart,jend
           do i=istart,iend
              if((source(1,j)-source(1,i))**2 +
     1           (source(2,j)-source(2,i))**2 + 
     2           (source(3,j)-source(3,i))**2.gt.1.0d-28) then
                 if(ifcharge.eq.1) then
                    call lpotfld3d(ifgrad,source(1,i),
     1              charge(i),source(1,j),pottmp,gradtmp)
                    if(ifpot.eq.1) pot(j) = pot(j)+
     1                                 pottmp
                    if(ifgrad.eq.1) then
                       grad(1,j) = grad(1,j)+gradtmp(1)
                       grad(2,j) = grad(2,j)+gradtmp(2)
                       grad(3,j) = grad(3,j)+gradtmp(3)
                    endif
                 endif
                 if(ifdipole.eq.1) then
                    call lpotfld3d_dp(ifgrad,source(1,i),
     1              dipstr(i),dipvec(1,i),source(1,j),pottmp,gradtmp)
                    if(ifpot.eq.1) pot(j) = pot(j)+ pottmp
                    if(ifgrad.eq.1) then
                       grad(1,j) = grad(1,j)+gradtmp(1)
                       grad(2,j) = grad(2,j)+gradtmp(2)
                       grad(3,j) = grad(3,j)+gradtmp(3)
                    endif
                 endif
              endif
           enddo
        enddo
c
        return
        end
c------------------------------------------------------------------    
 
      subroutine lfmm3dtarg_direct(istart,iend,jstart,jend,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     targ,ifpottarg,pottarg,ifgradtarg,gradtarg,flagsort)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the computed potential
c     and field strengths
c
c     INPUT arguments
c-------------------------------------------------------------------
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to update the potential and gradients
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to update the potential and gradients
c
c     source       in: double precision(3,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: double complex
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: double complex(ns)
c                 dipole strengths at the source locations
c
c     targ        in: double precision(3,nt)
c                 target locations
c
c     ifpottarg    in: Integer
c                  Flag for computing the potential. The potential
c                  will be updated if ifpottarg = 1
c
c     ifgradtarg  in: Integer
c                  Flag for computing the gradient.
c                  The gradient will be updated if
c                  ifgradtarg == 1
c
c     flagsort     in: Integer
c                  Flag for computing velocities and gradients
c                  at targets using point fmm. The potential/
c                  gradient at target i will be evaluated if 
c                  flagsort(i) = -1
c
c
c------------------------------------------------------------
c     OUTPUT
c
c   Updated potential and field at the targets
c   pottarg     out: double complex(*)
c               potential at the targets
c   gradtarg    out: double complex (3,*)
c               field at the targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole

        integer flagsort(1)

        double precision source(3,*)
        double complex charge(*),dipstr(*)
        double precision dipvec(3,*)

        integer ifpottarg,ifgradtarg,ifhesstarg
        double precision targ(3,*)
c
        double complex pottarg(*)
        double complex gradtarg(3,*)

        double complex pottmp,gradtmp(3)

        double complex hesstmp(3)
c
        ifhesstarg = 0
        ns = iend - istart + 1
        do j=jstart,jend
           if(flagsort(j).eq.-1) then
              do i=istart,iend
                 if((targ(1,j)-source(1,i))**2 +
     1              (targ(2,j)-source(2,i))**2 + 
     2              (targ(3,j)-source(3,i))**2.gt.1.0d-28) then
                    if(ifcharge.eq.1) then
                       call lpotfld3d(ifgradtarg,source(1,i),
     1                 charge(i),targ(1,j),pottmp,gradtmp)
                       if(ifpottarg.eq.1) pottarg(j) = pottarg(j)+
     1                                    pottmp
                       if(ifgradtarg.eq.1) then
                          gradtarg(1,j) = gradtarg(1,j)+gradtmp(1)
                          gradtarg(2,j) = gradtarg(2,j)+gradtmp(2)
                          gradtarg(3,j) = gradtarg(3,j)+gradtmp(3)
                       endif
                    endif
                    if(ifdipole.eq.1) then
                       call lpotfld3d_dp(ifgradtarg,source(1,i),
     1                 dipstr(i),dipvec(1,i),targ(1,j),pottmp,gradtmp)
                       if(ifpottarg.eq.1) pottarg(j) = pottarg(j)+
     1                                    pottmp
                       if(ifgradtarg.eq.1) then
                          gradtarg(1,j) = gradtarg(1,j)+gradtmp(1)
                          gradtarg(2,j) = gradtarg(2,j)+gradtmp(2)
                          gradtarg(3,j) = gradtarg(3,j)+gradtmp(3)
                       endif
                    endif
                 endif
              enddo
           endif
        enddo
c
        return
        end
c------------------------------------------------------------------    
 
