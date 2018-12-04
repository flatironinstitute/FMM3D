cc Copyright (C) 2017-2018: Leslie Greengard and
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


        subroutine hfmm3dpartstost(ier,iprec,zk,nsource,source,ifcharge,
     $    charge,ifdipole,dipstr,dipvec,ifpot,pot,iffld,fld,ntarg,
     $    targ,ifpottarg,pottarg,iffldtarg,fldtarg)
c       
c       
c       Generalized helmholtz FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets 
c
c       We use exp(ikr)/r for the Green's function., without
c       the 1/(4\pi ) scaling.
c
c       this is primarily a memory management code. 
c       The actual work is carried out in subroutine hfmm3dparttargmain
c   
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c
c   zk:    double complex: helmholtz parameter                
c
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
c   dipstr   in: double complex (nsource)
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
c                gradient will be evaluated at the targs if
c                iffldtarg = 1
c
c   ntarg  in: integer  
c                 number of targs 
c
c   targ  in: double precision (3,ntarg)
c               targ(k,j) is the kth component of the jth
c               targ location
c
c   ifpottarg   in: integer
c              flag for evaluating potential at the targs
c              The potential will be evaluated if ifpottarg = 1
c
c   iffldtarg   in:integer
c                flag for evaluating the gradient at the targs
c                gradient will be evaluated at the targs if
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
c   pottarg:    out: double complex(ntarg) 
c               potential at the targ locations
c
c   fldtarg:   out: double complex(3,ntarg)
c               gradient at the targ locations
     
c------------------------------------------------------------------

      implicit none

      double complex zk

      integer ier,iprec,ifcharge,ifdipole
      integer ifpot,iffld,ifpottarg,iffldtarg

      integer nsource,ntarg

      double precision source(3,*),targ(3,*)
      double complex charge(*)

      double complex dipstr(*)
      double precision dipvec(3,*)

      double complex pot(*),fld(3,*),pottarg(3,*),fldtarg(3,*)

c       Tree variables
      integer ltree,mhung,idivflag,ndiv,isep,nboxes,nbmax,nlevels
      integer nlmax
      integer mnbors,mnlist1,mnlist2,mnlist3,mnlist4
      integer ipointer(32)
      integer, allocatable :: itree(:)
      double precision, allocatable :: treecenters(:,:),boxsize(:)

c
cc      temporary sorted arrays
c
      double precision, allocatable :: sourcesort(:,:),targsort(:,:)
      double precision, allocatable :: radsrc(:)
      double complex, allocatable :: chargesort(:),dipstrsort(:)
      double precision, allocatable :: dipvecsort(:,:)

      integer, allocatable :: flagsort(:)

      double complex, allocatable :: potsort(:),fldsort(:,:)
      double complex, allocatable :: pottargsort(:),fldtargsort(:,:)

c
cc       temporary fmm arrays
c
      double precision epsfmm
      integer, allocatable :: nterms(:),iaddr(:,:)
      double precision, allocatable :: scales(:)
      double precision, allocatable :: rmlexp(:)

      integer lmptemp,nmax,lmptot
      double complex, allocatable :: mptemp(:),mptemp2(:)

c
cc       temporary variables not used in particle code
c
      double precision expc(3),texpssort(100),scjsort,radexp
      double precision expcsort(3),radssort
      integer ntj,nexpc,nadd

c
cc        other temporary variables
c
       integer i,iert,ifprint,ilev
       double precision time1,time2,omp_get_wtime,second

c
cc       figure out tree structure
c
   
c
cc        set criterion for box subdivision
c
       if(iprec.eq.1) ndiv = 100
       if(iprec.eq.2) ndiv = 100
       if(iprec.eq.3) ndiv = 100
       if(iprec.eq.4) ndiv = 100

c
cc         set tree flags
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

       allocate(flagsort(ntarg))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
       do i=1,ntarg
           flagsort(i) = -1
       enddo
C$OMP END PARALLEL DO  

       radexp = 0

c
cc      memory management code for constructing level restricted tree
        iert = 0


        call mklraptreemem(iert,source,nsource,radsrc,targ,ntarg,
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
        call mklraptree(source,nsource,radsrc,targ,ntarg,expc,
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
      ifprint=1
c
c     set fmm tolerance based on iprec flag.
c       
      if( iprec .eq. -2 ) epsfmm=.5d-0 
      if( iprec .eq. -1 ) epsfmm=.5d-1
      if( iprec .eq. 0 ) epsfmm=.5d-2
      if( iprec .eq. 1 ) epsfmm=.5d-3
      if( iprec .eq. 2 ) epsfmm=.5d-6
      if( iprec .eq. 3 ) epsfmm=.5d-9
      if( iprec .eq. 4 ) epsfmm=.5d-12
      if( iprec .eq. 5 ) epsfmm=.5d-15
      if( iprec .eq. 6 ) epsfmm=0
c      
      if(ifprint .eq. 1) call prin2('epsfmm=*',epsfmm,1)

c     Allocate sorted source and target arrays      

      allocate(sourcesort(3,nsource))
      allocate(targsort(3,ntarg))
      if(ifcharge.eq.1) allocate(chargesort(nsource))
      if(ifdipole.eq.1) then
         allocate(dipstrsort(nsource),dipvecsort(3,nsource))
      endif

      
      if(ifpot.eq.1) allocate(potsort(nsource))
      if(iffld.eq.1) allocate(fldsort(3,nsource))
      
      if(ifpottarg.eq.1) allocate(pottargsort(ntarg))
      if(iffldtarg.eq.1) allocate(fldtargsort(3,ntarg))


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
cc       initialize potential and field at targ
c        locations
c

      if(ifpottarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
         pottargsort(i) = 0.0d0
      enddo
C$OMP END PARALLEL DO
      endif

      if(iffldtarg.eq.1) then 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
         fldtargsort(1,i) = 0.0d0
         fldtargsort(2,i) = 0.0d0
         fldtargsort(3,i) = 0.0d0
      enddo
C$OMP END PARALLEL DO
      endif


c     Compute length of expansions at each level      
      nmax = 0
      do i=0,nlevels
         call h3dterms(boxsize(i),zk,epsfmm,nterms(i),ier)
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
cc       reorder sources
c
      call dreorderf(3,nsource,source,sourcesort,itree(ipointer(5)))
      if(ifcharge.eq.1) call dreorderf(2,nsource,charge,chargesort,
     1                     itree(ipointer(5)))

      if(ifdipole.eq.1) then
         call dreorderf(2,nsource,dipstr,dipstrsort,itree(ipointer(5)))
         call dreorderf(3,nsource,dipvec,dipvecsort,itree(ipointer(5)))
      endif

c
cc      reorder targs
c
      call dreorderf(3,ntarg,targ,targsort,itree(ipointer(6)))
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
c
      time1=second()
C$      time1=omp_get_wtime()
      call hfmm3dmain(ier,iprec,zk,
     $   nsource,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipstrsort,dipvecsort,
     $   ntarg,targsort,nexpc,expcsort,radssort,
     $   epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $   itree,ltree,ipointer,isep,ndiv,nlevels,
     $   nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $   scales,treecenters,itree(ipointer(1)),nterms,
     $   ifpot,potsort,iffld,fldsort,ifpottarg,pottargsort,
     $   iffldtarg,fldtargsort,flagsort,ntj,texpssort,scjsort,nadd)

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

      if(ifpottarg.eq.1) call dreorderi(2,ntarg,pottargsort,pottarg,
     1                 itree(ipointer(6)))
      if(iffldtarg.eq.1) call dreorderi(6,ntarg,fldtargsort,fldtarg,
     1                 itree(ipointer(6)))

   

      return
      end
c
c
c
c
      subroutine hfmm3dmain(ier,iprec,zk,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ntarg,targsort,nexpc,expcsort,radssort,
     $     epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $     itree,ltree,ipointer,isep,ndiv,nlevels, 
     $     nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $     scales,centers,laddr,nterms,
     $     ifpot,pot,iffld,fld,ifpottarg,pottarg,
     $     iffldtarg,fldtarg,flagsort,ntj,jsort,scjsort,nadd)
      implicit none

      integer ier,iprec
      integer nsource,ntarg,nexpc
      integer ndiv,nlevels

      integer ifcharge,ifdipole
      integer ifpot,iffld
      integer ifpottarg,iffldtarg
      double precision epsfmm

      double complex zk,zk2

      double precision sourcesort(3,nsource)

      double complex chargesort(nsource)
      double complex dipstrsort(nsource)
      double precision dipvecsort(3,nsource)

      double precision targsort(3,ntarg)
      integer flagsort(ntarg)

      double complex pot(*),fld(3,*)
      double complex pottarg(*),fldtarg(3,*)

      integer ntj
      double precision expcsort(3,nexpc)
      double complex jsort(0:ntj,-ntj:ntj,nexpc)
      integer nadd


      integer iaddr(2,nboxes), lmptot, lmptemp
      double precision rmlexp(lmptot)
      double precision mptemp(lmptemp)
      double precision mptemp2(lmptemp)
       
      double precision timeinfo(10)
      double precision centers(3,nboxes)
c
cc      tree variables
c
      integer isep, ltree
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer ipointer(32)
      integer itree(ltree)
      integer nboxes
      double precision scales(0:nlevels)
      double precision boxsize(0:nlevels)
c
cc      pw stuff
c

      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer uall(200),dall(200),nall(120),sall(120),eall(72),wall(72)
      integer u1234(36),d5678(36),n1256(24),s3478(24)
      integer e1357(16),w2468(16),n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(20),e3(5),e5(5),e7(5),w2(5),w4(5),w6(5),w8(5)

      integer ntmax, nexpmax, nlams, nmax, nthmax, nphmax
      parameter (ntmax = 1000)
      double precision, allocatable :: carray(:,:), dc(:,:)
      double precision, allocatable :: rdplus(:,:,:)
      double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, allocatable :: rdmsq3(:,:,:)
      double complex, allocatable :: rdminus2(:,:,:),zeyep(:)
      double complex, allocatable :: rdplus2(:,:,:)
      double precision, allocatable :: zmone(:)
      integer nn,nnn
  
      double complex rlams(ntmax), whts(ntmax)

      double complex, allocatable :: rlsc(:,:,:)
      integer nfourier(ntmax), nphysical(ntmax)
      integer nexptot, nexptotp
      double complex, allocatable :: xshift(:,:),yshift(:,:),zshift(:,:)

      double complex fexp(100000), fexpback(100000)

      double complex, allocatable :: mexp(:,:,:)
      double complex, allocatable :: tmp(:,:)
      double complex, allocatable :: mexpf1(:),mexpf2(:)
      double complex, allocatable :: mexpp1(:),mexpp2(:),mexppall(:,:)

      double precision scjsort(nexpc),radssort(nexpc)

c     temp variables
      integer i,j,k,l,ii,jj,kk,ll
      integer ibox,jbox,ilev,npts
      integer nchild,nlist1,nlist2,nlist3,nlist4

      integer istart,iend,istartt,iendt,istarte,iende
      integer istarts,iends
      integer jstart,jend

      integer ifprint

      integer ifhesstarg
      double precision d,time1,time2,omp_get_wtime
      double complex pottmp,fldtmp(3),hesstmp(3)

      double precision sourcetmp(3)
      double complex chargetmp

      integer ix,iy,iz
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype
      double precision wlege(40000)
      integer nterms_eval(4,0:nlevels)

      integer mnbors,mnlist1, mnlist2,mnlist3,mnlist4
      double complex eye, ztmp,zmult
      double precision alphaj
      integer ctr,ifinit2
      double precision, allocatable :: xnodes(:),wts(:)
      double precision radius
      integer nquad2
      integer maX_nodes
      double precision pi
      
      integer istart0,istart1,istartm1,nprin
      double precision rtmp1,rtmp2,rtmp3,rtmp4
      double complex ima
      data ima/(0.0d0,1.0d0)/

      integer nlfbox


      pi = 4.0d0*atan(1.0d0)

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
      enddo

      call prinf('nmax=*',nmax,1)
      call prinf('nterms=*',nterms,nlevels+1)
      call prinf('nlevels=*',nlevels,1)


      allocate(zeyep(-nmax:nmax),zmone(0:2*nmax))
      
      zeyep(0) = 1
      zmult = -ima
      do i=1,nmax
         zeyep(i) = zeyep(i-1)*zmult
         zeyep(-i) = zeyep(-i+1)/zmult
      enddo


      zmone(0) = 1
      do i=1,2*nmax
         zmone(i) = -zmone(i-1)
      enddo



c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
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
              jsort(j,k,i)=0
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
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            call h3dzero(rmlexp(iaddr(1,ibox)),nterms(ilev))
            call h3dzero(rmlexp(iaddr(2,ibox)),nterms(ilev))
         enddo
C$OMP END PARALLEL DO          
       enddo


c
ccc       set scjsort
c
      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(ipointer(3)+ibox-1)
            if(nchild.gt.0) then
               istart = itree(ipointer(16)+ibox-1)
               iend = itree(ipointer(17)+ibox-1)
               do i=istart,iend
                  scjsort(i) = scales(ilev)
                  radssort(i) = min(radssort(i),boxsize(ilev)/32*
     1                            sqrt(3.0d0))
               enddo
            endif
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
            call h3dterms_eval(itype,boxsize(ilev),zk,epsfmm,
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
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,npts,istart,iend,nchild)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            call h3dzero(rmlexp(iaddr(1,ibox)),nterms(ilev))

            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)
            npts = iend-istart+1

            nchild = itree(ipointer(3)+ibox-1)

            if(npts.gt.0.and.nchild.eq.0) then
               if(ifcharge.eq.1) then
                  call h3dformmp_add_trunc(ier,zk,scales(ilev),
     1            sourcesort(1,istart),chargesort(istart),npts,
     2            centers(1,ibox),nterms(ilev), nterms(ilev),
     3            rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
               if(ifdipole.eq.1) then
                  call h3dformmp_dp_add_trunc(ier,zk,scales(ilev),
     1            sourcesort(1,istart),dipstrsort(istart),
     2            dipvecsort(1,istart),          
     3            npts,centers(1,ibox),
     4            nterms(ilev),nterms(ilev),
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


      if(ifcharge.eq.1) then
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
               if(npts.gt.0) then
                  call h3dformta_add_trunc(ier,zk,scales(ilev),
     1             sourcesort(1,istart),chargesort(istart),npts,
     2             centers(1,ibox),nterms(ilev),nterms(ilev),
     3             rmlexp(iaddr(2,ibox)),wlege,nlege)
               endif
            enddo
         enddo
C$OMP END PARALLEL DO
      enddo
      endif


      if(ifdipole.eq.1) then
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
               if(npts.gt.0) then
                   call h3dformta_dp_add_trunc(ier,zk,scales(ilev),
     1              sourcesort(1,istart),dipstrsort(istart),
     2              dipvecsort(1,istart),npts,centers(1,ibox),
     3              nterms(ilev),nterms(ilev),
     4              rmlexp(iaddr(2,ibox)),wlege,nlege)
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      endif
      time2=second()
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

c       
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 3 (merge mp) ====*',i,0)
      time1=second()
C$    time1=omp_get_wtime()
c

      max_nodes = 10000
      allocate(xnodes(max_nodes))
      allocate(wts(max_nodes))

      do ilev=nlevels-1,0,-1
         nquad2 = nterms(ilev)*2.5
         nquad2 = max(6,nquad2)
         ifinit2 = 1
         call legewhts(nquad2,xnodes,wts,ifinit2)
         radius = boxsize(ilev)/2*sqrt(3.0d0)

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
                     call h3dmpmpquadu_add(zk,scales(ilev+1),
     1               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),scales(ilev),centers(1,ibox),
     3               rmlexp(iaddr(1,ibox)),nterms(ilev),nterms(ilev),
     4               radius,xnodes,wts,nquad2,ier)
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
c      ... step 3, convert multipole expansions into local
c       expansions

      time1 = second()
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels

c
cc       load the necessary quadrature for plane waves
c
      
         zk2 = zk*boxsize(ilev)
         if(real(zk2).le.pi.and.imag(zk2).le.0.02d0) then
            call lreadall(iprec,zk2,nlams,rlams,whts,nfourier,
     1           nphysical,ntmax,ier)


            nphmax = 0
            nthmax = 0
            nexptotp = 0
            nexptot = 0
            do i=1,nlams
               nexptotp = nexptotp + nphysical(i)
               nexptot = nexptot + 2*nfourier(i)+1
               if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
               if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
            enddo

            allocate(xshift(-5:5,nexptotp))
            allocate(yshift(-5:5,nexptotp))
            allocate(zshift(5,nexptotp))
            allocate(rlsc(0:nterms(ilev),0:nterms(ilev),nlams))
            allocate(tmp(0:nterms(ilev),-nterms(ilev):nterms(ilev)))
 
            allocate(mexpf1(nexptot),mexpf2(nexptot),mexpp1(nexptotp))
            allocate(mexpp2(nexptotp),mexppall(nexptotp,16))


c
cc      NOTE: there can be some memory savings here
c
            allocate(mexp(nexptotp,nboxes,6))

            nn = nterms(ilev)
            allocate(carray(4*nn+1,4*nn+1))
            allocate(dc(0:4*nn,0:4*nn))
            allocate(rdplus(0:nn,0:nn,-nn:nn))
            allocate(rdminus(0:nn,0:nn,-nn:nn))
            allocate(rdsq3(0:nn,0:nn,-nn:nn))
            allocate(rdmsq3(0:nn,0:nn,-nn:nn))

c     generate rotation matrices and carray
            call rotgen(nn,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)


            call rlscini(rlsc,nlams,rlams,zk2,nterms(ilev))
            call mkexps(rlams,nlams,nphysical,nexptotp,zk2,xshift,
     1           yshift,zshift)
            call mkfexp(nlams,nfourier,nphysical,fexp,fexpback)
c
cc      zero out mexp
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

c
cc         create multipole to plane wave expansion for
c          all boxes at this level
c
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,tmp,mexpf1,mexpf2,mptemp,ctr)
C$OMP$PRIVATE(ii,jj)
            do ibox = laddr(1,ilev),laddr(2,ilev)
               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = iend - istart+1
               if(npts.gt.0) then

c           rescale multipole expansion
                  ctr = 0
                  do ii = -nterms(ilev),nterms(ilev)
                     do jj = 0,nterms(ilev)

                        tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,ibox)+ctr),
     1                         rmlexp(iaddr(1,ibox)+ctr+1))
                       tmp(jj,ii) = tmp(jj,ii)*scales(ilev)**(jj)
                       ctr = ctr + 2

                     enddo
                  enddo

                  call mpoletoexp(tmp,nterms(ilev),
     1                  nlams,nfourier,nexptot,mexpf1,mexpf2,rlsc) 

                  call ftophys(mexpf1,nlams,nfourier,nphysical,
     1                 mexp(1,ibox,1),fexp)           

                  call ftophys(mexpf2,nlams,nfourier,nphysical,
     1                 mexp(1,ibox,2),fexp)


c             form mexpnorth, mexpsouth for current box

c             Rotate mpole for computing mexpnorth and
c             mexpsouth
                  call rotztoy(nterms(ilev),tmp,
     1                           mptemp,rdminus)

                  call mpoletoexp(mptemp,nterms(ilev),nlams,
     1                  nfourier,nexptot,mexpf1,mexpf2,rlsc)

                  call ftophys(mexpf1,nlams,nfourier,
     1                 nphysical,mexp(1,ibox,3),fexp)           

                  call ftophys(mexpf2,nlams,nfourier,
     1                 nphysical,mexp(1,ibox,4),fexp)   


c             Rotate mpole for computing mexpeast, mexpwest
                  call rotztox(nterms(ilev),tmp,
     1                              mptemp,rdplus)
                  call mpoletoexp(mptemp,nterms(ilev),nlams,
     1                  nfourier,nexptot,mexpf1,mexpf2,rlsc)

                  call ftophys(mexpf1,nlams,nfourier,
     1                 nphysical,mexp(1,ibox,5),fexp)

                  call ftophys(mexpf2,nlams,nfourier,
     1                 nphysical,mexp(1,ibox,6),fexp)           

               endif
            enddo
C$OMP END PARALLEL DO         


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

               nchild = itree(ipointer(3)+ibox-1)

               if(ifpot.eq.1.or.iffld.eq.1) then
                  istart = itree(ipointer(10)+ibox-1)
                  iend = itree(ipointer(11)+ibox-1)
                  npts = npts + iend-istart+1
               endif


               if(npts.gt.0.and.nchild.gt.0) then

              
                  call getpwlistall(ibox,boxsize(ilev),nboxes,
     1            itree(ipointer(18)+ibox-1),itree(ipointer(19)+
     2            mnbors*(ibox-1)),nchild,itree(ipointer(4)),centers,
     3            isep,nuall,uall,ndall,dall,nnall,nall,nsall,sall,
     4            neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678,
     5            nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468,
     6            nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,
     7            e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7,
     8            nw2,w2,nw4,w4,nw6,w6,nw8,w8)


                  call processudexp(zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(4)),scales(ilev),nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            nuall,uall,nu1234,u1234,ndall,dall,nd5678,d5678,
     5            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1),
     6            mexppall(1,2),mexppall(1,3),mexppall(1,4),xshift,
     7            yshift,zshift,fexpback,rlsc)


                  call processnsexp(zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(4)),scales(ilev),nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            nnall,nall,nn1256,n1256,nn12,n12,nn56,n56,nsall,sall,
     5            ns3478,s3478,ns34,s34,ns78,s78,
     6            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1),
     7            mexppall(1,2),mexppall(1,3),mexppall(1,4),
     8            mexppall(1,5),mexppall(1,6),mexppall(1,7),
     9            mexppall(1,8),rdplus,xshift,yshift,zshift,
     9            fexpback,rlsc)

                  call processewexp(zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(4)),scales(ilev),nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            neall,eall,ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,
     5            ne3,e3,ne5,e5,ne7,e7,nwall,wall,
     5            nw2468,w2468,nw24,w24,nw68,w68,
     5            nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1),
     7            mexppall(1,2),mexppall(1,3),mexppall(1,4),
     8            mexppall(1,5),mexppall(1,6),
     8            mexppall(1,7),mexppall(1,8),mexppall(1,9),
     9            mexppall(1,10),mexppall(1,11),mexppall(1,12),
     9            mexppall(1,13),mexppall(1,14),mexppall(1,15),
     9            mexppall(1,16),rdminus,xshift,yshift,zshift,
     9            fexpback,rlsc)
               endif
            enddo
C$OMP END PARALLEL DO        

            deallocate(xshift,yshift,zshift,rlsc,tmp)
            deallocate(carray,dc,rdplus,rdminus,rdsq3,rdmsq3)

            deallocate(mexpf1,mexpf2,mexpp1,mexpp2,mexppall,mexp)

         endif

         if(real(zk2).ge.pi.or.imag(zk2).ge.0.02d0) then

            nquad2 = nterms(ilev)*1.2
            nquad2 = max(6,nquad2)
            ifinit2 = 1
            ier = 0

            call legewhts(nquad2,xnodes,wts,ifinit2)

            radius = boxsize(ilev)/2*sqrt(3.0d0)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nlist2,i,jbox)
            do ibox = laddr(1,ilev),laddr(2,ilev)
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


               nlist2 = itree(ipointer(22)+ibox-1)
               if(npts.gt.0) then
                  do i =1,nlist2
                     jbox = itree(ipointer(23)+mnlist2*(ibox-1)+i-1)

                     istart = itree(ipointer(10)+jbox-1)
                     iend = itree(ipointer(11)+jbox-1)
                     npts = iend-istart+1

                     if(npts.gt.0) then
                        call h3dmplocquadu_add_trunc(zk,scales(ilev),
     1                  centers(1,jbox),
     1                  rmlexp(iaddr(1,jbox)),nterms(ilev),nterms(ilev),
     2                  scales(ilev),centers(1,ibox),
     2                  rmlexp(iaddr(2,ibox)),nterms(ilev),
     3                  nterms(ilev),radius,xnodes,wts,nquad2,ier)
                     endif
                  enddo
               endif
           enddo
C$OMP END PARALLEL DO        
         endif
      enddo
      time2 = second()
C$        time2=omp_get_wtime()
      timeinfo(4) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (split loc) ===*',i,0)

      time1 = second()
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels-1

        nquad2 = nterms(ilev)*2
        nquad2 = max(6,nquad2)
        ifinit2 = 1
        call legewhts(nquad2,xnodes,wts,ifinit2)
        radius = boxsize(ilev+1)/2*sqrt(3.0d0)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)

            npts = 0
            istart = itree(ipointer(12)+ibox-1)
            iend = itree(ipointer(13)+ibox-1)

            npts = npts + iend-istart+1

            istart = itree(ipointer(14)+ibox-1)
            iend = itree(ipointer(17)+ibox-1)
            npts = npts + iend-istart+1

            if(ifpot.eq.1.or.iffld.eq.1) then
               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = npts + iend-istart+1
            endif

            if(npts.gt.0) then
               do i=1,8
                  jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
                  if(jbox.gt.0) then
                     call h3dloclocquadu_add(zk,scales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),scales(ilev+1),centers(1,jbox),
     3                rmlexp(iaddr(2,jbox)),nterms(ilev+1),
     4                nterms(ilev+1),radius,xnodes,wts,nquad2,ier)
                  endif
               enddo
            endif
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

c
cc       shift mutlipole expansions to expansion center
c        (Note: this part is not relevant for particle codes.
c         It is relevant only for QBX codes)


      nquad2 = 2*ntj
      call legewhts(nquad2,xnodes,wts,ifinit2)
      do ilev=1,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,j,i,jbox)
C$OMP$PRIVATE(mptemp)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nlist3 = itree(ipointer(24)+ibox-1)

            istart = itree(ipointer(16)+ibox-1)
            iend = itree(ipointer(17)+ibox-1)

            do j=istart,iend
               do i=1,nlist3
                  jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
c
cc                  shift multipole expansion directly from box
c                   center to expansion center
                     call h3dmplocquadu_add_trunc(zk,scales(ilev+1),
     1               centers(1,jbox),
     1               rmlexp(iaddr(1,jbox)),nterms(ilev+1),
     2               nterms(ilev+1),scjsort(j),expcsort(1,j),
     2               jsort(0,-ntj,j),ntj,
     3               ntj,radssort(j),xnodes,wts,nquad2,ier)
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO  
      enddo

c
cc       evaluate multipole expansions at source locations
c

      if(ifpot.eq.1.or.iffld.eq.1) then
      do ilev=1,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,j,i,jbox)
C$OMP$PRIVATE(pottmp,fldtmp)
C$OMP$SCHEDULE(DYNAMIC)
      
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nlist3 = itree(ipointer(24)+ibox-1)
            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)

            npts = 1
            do j=istart,iend
               do i=1,nlist3
                  jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
                  pottmp = 0
                  fldtmp(1) = 0
                  fldtmp(2) = 0
                  fldtmp(3) = 0
                  call h3dmpevalall_trunc(zk,scales(ilev+1),
     $               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     $               nterms(ilev+1),nterms(ilev+1),
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
         enddo
C$OMP END PARALLEL DO
      enddo
      endif

c
cc        evaluate multipole expansions at target locations
c
      if(ifpottarg.eq.1.or.iffldtarg.eq.1) then
      do ilev=1,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,j,i,jbox)
C$OMP$PRIVATE(pottmp,fldtmp)
C$OMP$SCHEDULE(DYNAMIC)
      
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nlist3 = itree(ipointer(24)+ibox-1)
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
                     call h3dmpevalall_trunc(zk,scales(ilev+1),
     $                  centers(1,jbox),rmlexp(iaddr(1,jbox)),
     $                  nterms(ilev+1),nterms(ilev+1),
     $                  targsort(1,j),npts,ifpottarg,
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
      endif

      time2 = second()
C$        time2=omp_get_wtime()
      timeinfo(6) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== step 7 (eval lo) ===*',i,0)

c     ... step 7, evaluate all local expansions
c

      nquad2 = 2*ntj
      call legewhts(nquad2,xnodes,wts,ifinit2)
      time1 = second()
C$        time1=omp_get_wtime()
C

c
cc       shift local expansion to local epxanion at expansion centers
c        (note: this part is not relevant for particle codes.
c        it is relevant only for qbx codes)

      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,mptemp,i)
C$OMP$SCHEDULE(DYNAMIC)      
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
               istart = itree(ipointer(16)+ibox-1)
               iend = itree(ipointer(17)+ibox-1)
               do i=istart,iend

                  call h3dzero(mptemp,ntj)
                  call h3dloclocquadu_add(zk,scales(ilev),
     1             centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2             nterms(ilev),scales(ilev),expcsort(1,i),
     3             mptemp,ntj,ntj,radssort(i),xnodes,wts,nquad2,
     3             ier)

                  call h3dadd(mptemp,jsort(0,-ntj,i),ntj)
               enddo
            endif
         enddo
C$OMP END PARALLEL DO
      enddo

c
cc        evaluate local expansion at source locations
c
      if(ifpot.eq.1.or.iffld.eq.1) then
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,pottmp,fldtmp,i)
C$OMP$SCHEDULE(DYNAMIC)      
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               do i=istart,iend
                  pottmp = 0.0d0
                  fldtmp(1) = 0.0d0
                  fldtmp(2) = 0.0d0
                  fldtmp(3) = 0.0d0

                  call h3dtaevalall_trunc(zk,scales(ilev),
     1            centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2            nterms(ilev),nterms(ilev),sourcesort(1,i),1,
     3            ifpot,pottmp,iffld,fldtmp,
     4            wlege,nlege,ier)
                  if(ifpot.eq.1) pot(i) = pot(i)+pottmp
                  if(iffld.eq.1) then 
                     fld(1,i) = fld(1,i)+fldtmp(1)
                     fld(2,i) = fld(2,i)+fldtmp(2)
                     fld(3,i) = fld(3,i)+fldtmp(3)
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo
      endif

    
c
cc         evaluate local expansion at targret locations
c
      if(ifpottarg.eq.1.or.iffldtarg.eq.1) then
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,pottmp,fldtmp,i)
C$OMP$SCHEDULE(DYNAMIC)      
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
               istart = itree(ipointer(12)+ibox-1)
               iend = itree(ipointer(13)+ibox-1)
               do i=istart,iend
                  if(flagsort(i).eq.-1) then
                     pottmp = 0.0d0
                     fldtmp(1) = 0.0d0
                     fldtmp(2) = 0.0d0
                     fldtmp(3) = 0.0d0

                     call h3dtaevalall_trunc(zk,scales(ilev),
     1               centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2               nterms(ilev),nterms(ilev),targsort(1,i),1,
     3               ifpottarg,pottmp,iffldtarg,fldtmp,
     4               wlege,nlege,ier)
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
      endif

      time2 = second()
C$        time2=omp_get_wtime()
      timeinfo(7) = time2 - time1

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
      time1=second()
C$        time1=omp_get_wtime()

c
cc       directly form local expansions for list1 sources
c        at expansion centers. 
c        (note: this part is not relevant for particle codes.
c         It is relevant only for qbx codes)


      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarte,iende,nlist1,i,jbox)
C$OMP$PRIVATE(jstart,jend)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istarte = itree(ipointer(16)+ibox-1)
            iende = itree(ipointer(17)+ibox-1)

            nlist1 = itree(ipointer(20)+ibox-1)
   
            do i =1,nlist1
               jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)


               jstart = itree(ipointer(10)+jbox-1)
               jend = itree(ipointer(11)+jbox-1)

               call hfmm3dexpc_direct(zk,jstart,jend,istarte,
     1         iende,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,dipvecsort,expcsort,jsort,scjsort,ntj,
     3         wlege,nlege)
            enddo
         enddo
C$OMP END PARALLEL DO
      enddo

c
cc        directly evaluate potential at sources 
c         due to sources in list1

      if(ifpot.eq.1.or.iffld.eq.1) then
      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,nlist1,i,jbox)
C$OMP$PRIVATE(jstart,jend)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istarts = itree(ipointer(10)+ibox-1)
            iends = itree(ipointer(11)+ibox-1)

            nlist1 = itree(ipointer(20)+ibox-1)
   
            do i =1,nlist1
               jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
               jstart = itree(ipointer(10)+jbox-1)
               jend = itree(ipointer(11)+jbox-1)

               call hfmm3dsrc_direct(zk,jstart,jend,
     1         istarts,iends,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,dipvecsort,ifpot,pot,iffld,fld)
            enddo   
         enddo
C$OMP END PARALLEL DO     
      enddo
      endif
 

c
cc        directly evaluate potential at targets 
c         due to sources in list1

      if(ifpottarg.eq.1.or.iffldtarg.eq.1) then
      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,nlist1,i,jbox)
C$OMP$PRIVATE(jstart,jend)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istartt = itree(ipointer(12)+ibox-1)
            iendt = itree(ipointer(13)+ibox-1)

            nlist1 = itree(ipointer(20)+ibox-1)
   
            do i =1,nlist1
               jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)


               jstart = itree(ipointer(10)+jbox-1)
               jend = itree(ipointer(11)+jbox-1)

               call hfmm3dtarg_direct(zk,jstart,jend,
     1         istartt,iendt,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,dipvecsort,targsort,ifpottarg,pottarg,
     3         iffldtarg,fldtarg,flagsort)
            enddo   
         enddo
C$OMP END PARALLEL DO     
      enddo
      endif
 
      time2 = second()
C$        time2=omp_get_wtime()
      timeinfo(8) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
      d = 0
      do i = 1,8
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

      return
      end
c------------------------------------------------------------------
      subroutine hfmm3dexpc_direct(zk,istart,iend,jstart,
     $     jend,source,ifcharge,charge,ifdipole,dipstr,
     $     dipvec,targ,texps,scj,ntj,wlege,nlege)
c---------------------------------------------------------------
c     This subroutine adds the local expansions due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the existing local
c     expansions
c
c     INPUT arguments
c-------------------------------------------------------------------
c     zk           in: double complex
c                  helmholtz parameter
c
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
c                  wish to compute the expansions
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to compute the expansions
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
c     targ        in: double precision(3,nexpc)
c                 Expansion center locations
c
c     scj         in: double precision(nexpc)
c                 scaling parameters for local expansions
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
c   texps : coeffs for local expansions
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j, nlege
        integer ifcharge,ifdipole,ier
        double complex zk
        double precision source(3,*)
        double precision wlege(0:nlege,0:nlege)
        double complex charge(*),dipstr(*)
        double precision dipvec(3,*)
        double precision targ(3,*),scj(*)

        integer nlevels,ntj
c
        double complex texps(0:ntj,-ntj:ntj,*)
        
c
        ns = iend - istart + 1
        do j=jstart,jend
           if(ifcharge.eq.1) then
              call h3dformta_add_trunc(ier,zk,scj(j),
     1        source(1,istart),charge(istart),ns,targ(1,j),
     2        ntj,ntj,texps(0,-ntj,j),wlege,nlege)
           endif

           if(ifdipole.eq.1) then
               call h3dformta_dp_add_trunc(ier,zk,scj(j),
     1         source(1,istart),dipstr(istart),dipvec(1,istart),
     2         ns,targ(1,j),ntj,ntj,texps(0,-ntj,j),wlege,nlege)
           endif        
        enddo
c
        return
        end
c------------------------------------------------------------------     
      subroutine hfmm3dsrc_direct(zk,istart,iend,jstart,
     $     jend,source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to sources
c     istart to iend in the source array at the sources
c     jstart to jend in the source array to the potential
c     and gradients
c
c     INPUT arguments
c-------------------------------------------------------------------
c     zk           in: double complex
c                  helmholtz parameter
c
c     istart       in:Integer
c                  Starting index in source array 
c
c     iend         in:Integer
c                  Last index in source array 
c
c     jstart       in: Integer
c                  First index in source array at which we
c                  wish to update the potential and gradients
c 
c     jend         in:Integer
c                  Last index in source array at which we wish
c                  to update the potential and gradients
c
c     source       in: double precision(3,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including contribution due to charges
c                  The charge contribution will be included
c                  if ifcharge == 1
c
c     charge       in: double complex
c                  Charge strength at the source locations
c
c     ifdipole     in: Integer
c                 flag for including contribution due to dipoles
c                 The dipole contributions will be included
c                 if ifdipole == 1
c
c     dipstr        in: double complex(ns)
c                 dipole strengths at the source locations
c
c     dipvec       in: double precision(3,ns)
c                  dipole orientation at source locations
c
c
c     ifpot    in: Integer
c                  Flag for updating the potential. The potential
c                  will be updated if ifpot = 1
c
c     ifgrad  in: Integer
c                  Flag for updating the gradient.
c                  The gradient will be updated if
c                  ifgrad == 1
c
c
c------------------------------------------------------------
c     OUTPUT
c
c   Updated potential and field strengths at the sources
c   pot :      out: double  complex(*)
c              potential at the sources
c
c   grad:      out: double complex(*)
c              gradient at the sources
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole

 
        double complex zk
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
                    pottmp = 0.0d0
                    gradtmp(1) = 0.0d0
                    gradtmp(2) = 0.0d0
                    gradtmp(3) = 0.0d0
                    call hpotfld3d(ifgrad,source(1,i),
     1              charge(i),source(1,j),zk,pottmp,gradtmp)
                    if(ifpot.eq.1) pot(j) = pot(j)+pottmp
                    if(ifgrad.eq.1) then
                       grad(1,j) = grad(1,j)+gradtmp(1)
                       grad(2,j) = grad(2,j)+gradtmp(2)
                       grad(3,j) = grad(3,j)+gradtmp(3)
                    endif
                 endif
                 if(ifdipole.eq.1) then
                    pottmp = 0.0d0
                    gradtmp(1) = 0.0d0
                    gradtmp(2) = 0.0d0
                    gradtmp(3) = 0.0d0
                    call hpotfld3d_dp(ifgrad,source(1,i),
     1              dipstr(i),dipvec(1,i),source(1,j),zk,pottmp,
     2              gradtmp)
                    if(ifpot.eq.1) pot(j) = pot(j)+pottmp
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
c----------------------------------------------------        
      subroutine hfmm3dtarg_direct(zk,istart,iend,jstart,
     $     jend,source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     targ,ifpottarg,pottarg,ifgradtarg,gradtarg,flagsort)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the computed velocities
c     and gradients
c
c     INPUT arguments
c-------------------------------------------------------------------
c     zk           in: double complex
c                  helmholtz parameter
c
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
c------------------------------------------------------------
c     OUTPUT
c
c   Updated velocity and gradients at the targets
c   pottarg : potential at the targets
c   gradtarg: gradient at the targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole

        integer flagsort(1)
 
        double complex zk
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
                       pottmp = 0.0d0
                       gradtmp(1) = 0.0d0
                       gradtmp(2) = 0.0d0
                       gradtmp(3) = 0.0d0
                       call hpotfld3d(ifgradtarg,source(1,i),
     1                 charge(i),targ(1,j),zk,pottmp,gradtmp)
                       if(ifpottarg.eq.1) pottarg(j) = pottarg(j)+
     1                                    pottmp
                       if(ifgradtarg.eq.1) then
                          gradtarg(1,j) = gradtarg(1,j)+gradtmp(1)
                          gradtarg(2,j) = gradtarg(2,j)+gradtmp(2)
                          gradtarg(3,j) = gradtarg(3,j)+gradtmp(3)
                       endif
                    endif
                    if(ifdipole.eq.1) then
                       pottmp = 0.0d0
                       gradtmp(1) = 0.0d0
                       gradtmp(2) = 0.0d0
                       gradtmp(3) = 0.0d0
                       call hpotfld3d_dp(ifgradtarg,source(1,i),
     1                 dipstr(i),dipvec(1,i),targ(1,j),zk,pottmp,
     2                 gradtmp)
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
c----------------------------------------------------        
