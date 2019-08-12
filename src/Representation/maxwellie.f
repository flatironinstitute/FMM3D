c  This file contains the Maxwell IE FMM wrappers

      subroutine guruie_eh(nd, eps, zk, 
     $                 nsource, source,
     $                 ifj, j, ifr, r, ifm, m, ifq, q,
     $                 ifeh, e, h, grade, gradh,
     $                 ntarg, targ, 
     $                 ifehtarg, etarg, htarg, gradetarg, gradhtarg)
c
c        Guru interface for Maxwell's equations, 
c        FMM in R^{3}: evaluate all pairwise particle
c        interactions (ignoring self-interactions) and interactions
c        with targs.
c
c        We use exp(ikr)/r for the Helmholtz Green's function, without the
c        1/(4 \pi) scaling.
c
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:    in: integer
c              number of densities
c   
c   eps:   in: double precision
c              requested precision
c
c   zk:    in: double complex
c              helmholtz parameter                
c
c   nsource in: integer  
c               number of sources
c
c   source  in: double precision (3,nsource)
c               source(k,j) is the kth component of the jth
c               source locations
c
c   ifj     in: integer  
c               fictitious surface electric current computation flag
c               ifj = 1   =>  include surface current contribution
c                                otherwise do not
c 
c   j       in: double complex (nd,3,nsource) 
c               fictitious surface electric current
c
c   ifr     in: integer
c               fictitious electric charge computation flag
c               ifrho = 1   =>  include charge contribution
c                                   otherwise do not
c
c   r       in: double complex (nd,3,nsource) 
c               fictitious electric charge strengths
c
c   ifm     in: integer  
c               fictitious surface magnetic current computation flag
c               ifj = 1   =>  include surface current contribution
c                                otherwise do not
c 
c   m       in: double complex (nd,3,nsource) 
c               fictitious surface magnetic current
c
c   ifq     in: integer
c               fictitious magnetic charge computation flag
c               ifrho = 1   =>  include charge contribution
c                                   otherwise do not
c
c   q       in: double complex (nd,3,nsource) 
c               fictitious magnetic charge strengths
c
c   ifeh    in: integer
c               flag for evaluating potential/gradient at the sources
c               ifeh = 1, e potential is evaluated
c               ifeh = 2, e potential and e gradients are evaluated
c               ifeh = 3, h potential is evaluated
c               ifeh = 4, h potential and h gradients are evaluated
c               ifeh = 5, e potential and h potential are evaluated
c               ifeh = 6, e h potential and e h gradients are evaluated
c
c   ntarg   in: integer  
c              number of targs 
c
c   targ    in: double precision (3,ntarg)
c             targ(k,j) is the kth component of the jth
c             targ location
c
c   ifehtarg  in: integer
c                  flag for evaluating potential/gradient at the targs
c                  ifehtarg = 1, e potential is evaluated
c                  ifehtarg = 2, e potential and e gradients are evaluated
c                  ifehtarg = 3, h potential is evaluated
c                  ifehtarg = 4, h potential and h gradients are evaluated
c                  ifehtarg = 5, e potential and h potential are evaluated
c                  ifehtarg = 6, e h potential and e h gradients are evaluated
c
c
c   OUTPUT parameters:
c
c   e:    out: double complex(nd,3,nsource) 
c              e potential at the source locations
c
c   h:    out: double complex(nd,3,nsource) 
c              h potential at the source locations
c
c   grade:   out: double complex(nd,3,3,nsource)
c                 e gradient at the source locations
c
c   gradh:   out: double complex(nd,3,3,nsource)
c                 h gradient at the source locations
c
c   etarg:    out: double complex(nd,3,ntarg) 
c                  e potential at the target locations
c
c   htarg:    out: double complex(nd,3,ntarg) 
c                  h potential at the target locations
c
c   gradetarg:   out: double complex(nd,3,3,ntarg)
c                     e gradient at the target locations
c
c   gradhtarg:   out: double complex(nd,3,3,ntarg)
c                     h gradient at the target locations
c
c------------------------------------------------------------------
      implicit none
      integer nd
      double precision eps
      double complex zk
      integer nsource, ifj, ifr, ifm, ifq, ifeh, ntarg, ifehtarg
      double precision source(3, *), targ(3, *)
      double complex j(nd, 3, *), r(nd, *)
      double complex m(nd, 3, *), q(nd, *)
      double complex e(nd, 3, *), h(nd, 3, *)
      double complex grade(nd, 3, 3, *), gradh(nd, 3, 3, *)
      double complex etarg(nd, 3, *), htarg(nd, 3, *)
      double complex gradetarg(nd, 3, 3, *), gradhtarg(nd, 3, 3, *)

      !work variables
      double complex ima,coef
      integer ndi,dimi,dimj,srci,cnt,trgi
      integer ifpgh, ifpghtarg
      integer ifcharge, ifdipole
      !work array
      double complex, allocatable :: charge(:,:), cpot(:,:)
      double complex, allocatable :: dipvec(:, :, :)
      double complex, allocatable :: cgrad(:,:,:)
      double complex, allocatable :: cgradtarg(:,:,:)
      double complex, allocatable :: cpottarg(:,:)
      double complex hess(nd*3,6)
      double complex hesstarg(nd*3,6)
      !allocate memory space for work array
      !may save the memory space by checking if flags ahead
      allocate(charge(nd*3,nsource))
      allocate(cpot(nd*3,nsource))
      allocate(dipvec(nd*3,3,nsource))
      allocate(cgrad(nd*3,3,nsource))
      allocate(cpottarg(nd*3,ntarg))
      allocate(cgradtarg(nd*3,3,ntarg))

      !coeffient factor needed
      data ima/(0.0d0,1.0d0)/
      coef=ima*zk

      !source to source,target
      ! e = S_kj - \grad phi - \curl S_km
      dipvec=(0d0,0d0)
      charge=(0d0,0d0)
      cpot=(0d0,0d0)
      cgrad=(0d0,0d0)
      cpottarg=(0d0,0d0)
      cgradtarg=(0d0,0d0)
      if((ifeh.eq.1 .or. ifeh.eq.2 .or. ifeh.eq.5 .or. ifeh.eq.6) .or.
     $   (ifehtarg.eq.1 .or. ifehtarg.eq.2 
     $    .or. ifehtarg.eq.5 .or. ifehtarg.eq.6)
     $  ) then
        ifcharge=0
        ifdipole=0
        !pack charge and dipole
        !pack charge
        if(ifj .eq. 1) then
          ifcharge=1
          !fortran is column major storage
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,ndi,cnt)
          do srci=1,nsource
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          charge(cnt,srci) = j(ndi,dimi,srci)*coef
          cnt=cnt+1
          enddo
          enddo
          enddo
C$OMP END PARALLEL DO
        endif
        !pack dipole
        if(ifr .eq. 1) then
          ifdipole=1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,dimj,ndi,cnt)
          do srci=1,nsource
          do dimj=1,3
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          if(dimi .eq. dimj) then
            dipvec(cnt,dimj,srci) = r(ndi, srci)
          endif
          cnt=cnt+1
          enddo
          enddo
          enddo
          enddo
C$OMP END PARALLEL DO
        endif
        !continue pack dipole
        if(ifm .eq. 1) then
          ifdipole=1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,ndi,cnt)
          do srci=1,nsource
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          if(dimi .eq. 1) then
            dipvec(cnt,1,srci) = dipvec(cnt,1,srci) + 0d0
            dipvec(cnt,2,srci) = dipvec(cnt,2,srci) + m(ndi,3,srci)
            dipvec(cnt,3,srci) = dipvec(cnt,3,srci) - m(ndi,2,srci)
          endif
          if(dimi .eq. 2) then
            dipvec(cnt,1,srci) = dipvec(cnt,1,srci) - m(ndi,3,srci)
            dipvec(cnt,2,srci) = dipvec(cnt,2,srci) + 0d0
            dipvec(cnt,3,srci) = dipvec(cnt,3,srci) + m(ndi,1,srci)
          endif
          if(dimi .eq. 3) then
            dipvec(cnt,1,srci) = dipvec(cnt,1,srci) + m(ndi,2,srci)
            dipvec(cnt,2,srci) = dipvec(cnt,2,srci) - m(ndi,1,srci)
            dipvec(cnt,3,srci) = dipvec(cnt,3,srci) + 0d0
          endif
          cnt=cnt+1
          enddo
          enddo 
          enddo
        endif
        ! set hfmm3d flag and call hfmm3d
        ifpgh=0
        ifpghtarg=0
        if(ifeh.eq.1 .or. ifeh.eq.5) then
          ifpgh=1
        endif
        if(ifeh.eq.2 .or. ifeh.eq.6) then
          ifpgh=2
        endif
        if(ifehtarg.eq.1 .or. ifehtarg.eq.5) then
          ifpghtarg=1
        endif
        if(ifehtarg.eq.2 .or. ifehtarg.eq.6) then
          ifpghtarg=2
        endif
        call hfmm3d(nd*3, eps, zk, nsource, source, ifcharge, charge,
     $         ifdipole, dipvec, ifpgh, cpot, cgrad, hess, ntarg, targ,
     $         ifpghtarg, cpottarg, cgradtarg, hesstarg)
        !unpack source output from cpot  and cgrad
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,ndi,cnt)
        do srci=1,nsource
        cnt=1
        do dimi=1,3
        do ndi=1,nd
        if(ifeh.eq.1 .or. ifeh.eq.5 .or. ifeh.eq.2 .or. ifeh.eq.6) then
        e(ndi,dimi,srci) = cpot(cnt,srci)
        endif
        if(ifeh.eq.2 .or. ifeh.eq.6) then
        grade(ndi,dimi,1,srci) = cgrad(cnt,1,srci)
        grade(ndi,dimi,2,srci) = cgrad(cnt,2,srci)
        grade(ndi,dimi,3,srci) = cgrad(cnt,3,srci)
        endif
        cnt=cnt+1
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
        !unpack targ output from cpottarg and cgradtarg
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(trgi,dimi,ndi,cnt)
        do trgi=1,ntarg
        cnt=1
        do dimi=1,3
        do ndi=1,nd
        if(ifehtarg.eq.1 .or. ifehtarg.eq.5 .or.
     $     ifehtarg.eq.2 .or. ifehtarg.eq.6) then
        etarg(ndi,dimi,trgi) = cpottarg(cnt,trgi)
        endif
        if(ifehtarg.eq.2 .or. ifehtarg.eq.6) then
        gradetarg(ndi,dimi,1,trgi) = cgradtarg(cnt,1,trgi)
        gradetarg(ndi,dimi,2,trgi) = cgradtarg(cnt,2,trgi)
        gradetarg(ndi,dimi,3,trgi) = cgradtarg(cnt,3,trgi)
        endif
        cnt=cnt+1
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
      endif

      !source to source,target
      ! h = \curl S_kj + S_km -\grad phi_m
      dipvec=(0d0,0d0)
      charge=(0d0,0d0)
      cpot=(0d0,0d0)
      cgrad=(0d0,0d0)
      cpottarg=(0d0,0d0)
      cgradtarg=(0d0,0d0)
      if((ifeh.eq.3 .or. ifeh.eq.4 .or. ifeh.eq.5 .or. ifeh.eq.6) .or. 
     $   (ifehtarg.eq.3 .or. ifehtarg.eq.4 .or.
     $    ifehtarg.eq.5 .or. ifehtarg.eq.6)
     $  ) then
        ifcharge=0
        ifdipole=0
        !pack charge and dipole
        !pack dipole
        if(ifj .eq. 1) then
          ifdipole=1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,ndi,cnt,dimi)
          do srci=1,nsource
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          if(dimi .eq. 1) then
            dipvec(3*(ndi-1)+1,1,srci) = 0d0
            dipvec(3*(ndi-1)+1,2,srci) = -j(ndi,3,srci)
            dipvec(3*(ndi-1)+1,3,srci) = j(ndi,2,srci)
          endif
          if(dimi .eq. 2) then
            dipvec(3*(ndi-1)+2,1,srci) = j(ndi,3,srci)
            dipvec(3*(ndi-1)+2,2,srci) = 0d0
            dipvec(3*(ndi-1)+2,3,srci) = -j(ndi,1,srci)
          endif
          if(dimi .eq. 3) then
            dipvec(3*(ndi-1)+3,1,srci) = -j(ndi,2,srci)
            dipvec(3*(ndi-1)+3,2,srci) = j(ndi,1,srci)
            dipvec(3*(ndi-1)+3,3,srci) = 0d0
          endif
          cnt=cnt+1
          enddo
          enddo
          enddo
C$OMP END PARALLEL DO
        endif
        !continue pack dipole
        if(ifq .eq. 1) then
          ifdipole=1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,dimj,ndi,cnt)
          do srci=1,nsource
          do dimj=1,3
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          if(dimi .eq. dimj) then
            dipvec(cnt,dimj,srci) = dipvec(cnt,dimj,srci)+q(ndi,srci)
          endif
          cnt=cnt+1
          enddo
          enddo
          enddo
          enddo
C$OMP END PARALLEL DO
        endif
        !pack charge
        if(ifm .eq. 1) then
          ifcharge=1
          do srci=1,nsource
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          charge(cnt,srci) = m(ndi,dimi,srci)*coef
          cnt=cnt+1
          enddo
          enddo
          enddo
        endif
        !set hfmm3d flags and call hfmm3d
        ifpgh=0
        ifpghtarg=0
        if(ifeh.eq.3 .or. ifeh.eq.5) then
          ifpgh=1
        endif
        if(ifeh.eq.4 .or. ifeh.eq.6) then
          ifpgh=2
        endif
        if(ifehtarg.eq.3 .or. ifehtarg.eq.5) then
          ifpghtarg=1
        endif
        if(ifehtarg.eq.4 .or. ifehtarg.eq.6) then
          ifpghtarg=2
        endif
        call hfmm3d(nd*3, eps, zk, nsource, source, ifcharge, charge,
     $         ifdipole, dipvec, ifpgh, cpot, cgrad, hess, ntarg, targ,
     $         ifpghtarg, cpottarg, cgradtarg, hesstarg)
        !unpack source output from cpot and cgrad
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,ndi,cnt)
        do srci=1,nsource
        cnt=1
        do dimi=1,3
        do ndi=1,nd
        if(ifeh.eq.3 .or. ifeh.eq.5 .or. ifeh.eq.4 .or. ifeh.eq.6) then
        h(ndi,dimi,srci) = cpot(cnt,srci)
        endif
        if(ifeh.eq.4 .or. ifeh.eq.6) then
        gradh(ndi,dimi,1,srci) = cgrad(cnt,1,srci)
        gradh(ndi,dimi,2,srci) = cgrad(cnt,2,srci)
        gradh(ndi,dimi,3,srci) = cgrad(cnt,3,srci)
        endif
        cnt=cnt+1
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
        !unpack targ output from cpottarg and cgradtarg
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(trgi,dimi,ndi,cnt)
        do trgi=1,ntarg
        cnt=1
        do dimi=1,3
        do ndi=1,nd
        if(ifehtarg.eq.3 .or. ifehtarg.eq.5 .or.
     $     ifehtarg.eq.4 .or. ifehtarg.eq.6) then
        htarg(ndi,dimi,trgi) = cpottarg(cnt,trgi)
        endif
        if(ifehtarg.eq.4 .or. ifehtarg.eq.6) then
        gradhtarg(ndi,dimi,1,trgi) = cgradtarg(cnt,1,trgi)
        gradhtarg(ndi,dimi,2,trgi) = cgradtarg(cnt,2,trgi)
        gradhtarg(ndi,dimi,3,trgi) = cgradtarg(cnt,3,trgi)
        endif
        cnt=cnt+1
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
      endif

      return
      end


      subroutine guruie_eh_direct(nd, eps, zk, 
     $                 nsource, source,
     $                 ifj, j, ifr, r, ifm, m, ifq, q,
     $                 ifeh, e, h, grade, gradh,
     $                 ntarg, targ, 
     $                 ifehtarg, etarg, htarg, gradetarg, gradhtarg)
c
c        Guru interface for Maxwell's equations, direct evaluation
c        in R^{3}: evaluate all pairwise particle
c        interactions (ignoring self-interactions) and interactions
c        with targs.
c
c        We use exp(ikr)/r for the Helmholtz Green's function, without the
c        1/(4 \pi) scaling.
c
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:    in: integer
c              number of densities
c   
c   eps:   in: double precision
c              requested precision
c
c   zk:    in: double complex
c              helmholtz parameter                
c
c   nsource in: integer  
c               number of sources
c
c   source  in: double precision (3,nsource)
c               source(k,j) is the kth component of the jth
c               source locations
c
c   ifj     in: integer  
c               fictitious surface electric current computation flag
c               ifj = 1   =>  include surface current contribution
c                                otherwise do not
c 
c   j       in: double complex (nd,3,nsource) 
c               fictitious surface electric current
c
c   ifr     in: integer
c               fictitious electric charge computation flag
c               ifrho = 1   =>  include charge contribution
c                                   otherwise do not
c
c   r       in: double complex (nd,3,nsource) 
c               fictitious electric charge strengths
c
c   ifm     in: integer  
c               fictitious surface magnetic current computation flag
c               ifj = 1   =>  include surface current contribution
c                                otherwise do not
c 
c   m       in: double complex (nd,3,nsource) 
c               fictitious surface magnetic current
c
c   ifq     in: integer
c               fictitious magnetic charge computation flag
c               ifrho = 1   =>  include charge contribution
c                                   otherwise do not
c
c   q       in: double complex (nd,3,nsource) 
c               fictitious magnetic charge strengths
c
c   ifeh    in: integer
c               flag for evaluating potential/gradient at the sources
c               ifeh = 1, e potential is evaluated
c               ifeh = 2, e potential and e gradients are evaluated
c               ifeh = 3, h potential is evaluated
c               ifeh = 4, h potential and h gradients are evaluated
c               ifeh = 5, e potential and h potential are evaluated
c               ifeh = 6, e h potential and e h gradients are evaluated
c
c   ntarg   in: integer  
c              number of targs 
c
c   targ    in: double precision (3,ntarg)
c             targ(k,j) is the kth component of the jth
c             targ location
c
c   ifehtarg  in: integer
c                  flag for evaluating potential/gradient at the targs
c                  ifehtarg = 1, e potential is evaluated
c                  ifehtarg = 2, e potential and e gradients are evaluated
c                  ifehtarg = 3, h potential is evaluated
c                  ifehtarg = 4, h potential and h gradients are evaluated
c                  ifehtarg = 5, e potential and h potential are evaluated
c                  ifehtarg = 6, e h potential and e h gradients are evaluated
c
c
c   OUTPUT parameters:
c
c   e:    out: double complex(nd,3,nsource) 
c              e potential at the source locations
c
c   h:    out: double complex(nd,3,nsource) 
c              h potential at the source locations
c
c   grade:   out: double complex(nd,3,3,nsource)
c                 e gradient at the source locations
c
c   gradh:   out: double complex(nd,3,3,nsource)
c                 h gradient at the source locations
c
c   etarg:    out: double complex(nd,3,ntarg) 
c                  e potential at the target locations
c
c   htarg:    out: double complex(nd,3,ntarg) 
c                  h potential at the target locations
c
c   gradetarg:   out: double complex(nd,3,3,ntarg)
c                     e gradient at the target locations
c
c   gradhtarg:   out: double complex(nd,3,3,ntarg)
c                     h gradient at the target locations
c
c------------------------------------------------------------------
      implicit none
      integer nd
      double precision eps
      double complex zk
      integer nsource, ifj, ifr, ifm, ifq, ifeh, ntarg, ifehtarg
      double precision source(3, *), targ(3, *)
      double complex j(nd, 3, *), r(nd, *)
      double complex m(nd, 3, *), q(nd, *)
      double complex e(nd, 3, *), h(nd, 3, *)
      double complex grade(nd, 3, 3, *), gradh(nd, 3, 3, *)
      double complex etarg(nd, 3, *), htarg(nd, 3, *)
      double complex gradetarg(nd, 3, 3, *), gradhtarg(nd, 3, 3, *)

      !work variables
      double complex ima,coef
      integer ndi,dimi,dimj,srci,cnt,trgi
      integer ifpgh, ifpghtarg
      integer ifcharge, ifdipole
      double precision thresh
      !work array
      double complex, allocatable :: charge(:,:), cpot(:,:)
      double complex, allocatable :: dipvec(:, :, :)
      double complex, allocatable :: cgrad(:,:,:)
      double complex, allocatable :: cgradtarg(:,:,:)
      double complex, allocatable :: cpottarg(:,:)
      double complex hess(nd*3,6)
      double complex hesstarg(nd*3,6)
      !allocate memory space for work array
      !may save the memory space by checking if flags ahead
      allocate(charge(nd*3,nsource))
      allocate(cpot(nd*3,nsource))
      allocate(dipvec(nd*3,3,nsource))
      allocate(cgrad(nd*3,3,nsource))
      allocate(cpottarg(nd*3,ntarg))
      allocate(cgradtarg(nd*3,3,ntarg))

      !coeffient factor needed
      data ima/(0.0d0,1.0d0)/
      coef=ima*zk
      thresh = 2.00d0**(-52)

      !source to source,target
      !S_kJ - \grad phi
      dipvec=(0d0,0d0)
      charge=(0d0,0d0)
      cpot=(0d0,0d0)
      cgrad=(0d0,0d0)
      cpottarg=(0d0,0d0)
      cgradtarg=(0d0,0d0)
      if((ifeh.eq.1 .or. ifeh.eq.2 .or. ifeh.eq.5 .or. ifeh.eq.6) .or.
     $   (ifehtarg.eq.1 .or. ifehtarg.eq.2 
     $    .or. ifehtarg.eq.5 .or. ifehtarg.eq.6)
     $  ) then
        ifcharge=0
        ifdipole=0
        !pack charge and dipole
        !pack charge
        if(ifj .eq. 1) then
          ifcharge=1
          !fortran is column major storage
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,ndi,cnt)
          do srci=1,nsource
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          charge(cnt,srci) = j(ndi,dimi,srci)*coef
          cnt=cnt+1
          enddo
          enddo
          enddo
C$OMP END PARALLEL DO
        endif
        !pack dipole
        if(ifr .eq. 1) then
          ifdipole=1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,dimj,ndi,cnt)
          do srci=1,nsource
          do dimj=1,3
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          if(dimi .eq. dimj) then
            dipvec(cnt,dimj,srci) = r(ndi, srci)
          endif
          cnt=cnt+1
          enddo
          enddo
          enddo
          enddo
C$OMP END PARALLEL DO
        endif
        !continue pack dipole
        if(ifm .eq. 1) then
          ifdipole=1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,ndi,cnt)
          do srci=1,nsource
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          if(dimi .eq. 1) then
            dipvec(cnt,1,srci) = dipvec(cnt,1,srci) + 0d0
            dipvec(cnt,2,srci) = dipvec(cnt,2,srci) + m(ndi,3,srci)
            dipvec(cnt,3,srci) = dipvec(cnt,3,srci) - m(ndi,2,srci)
          endif
          if(dimi .eq. 2) then
            dipvec(cnt,1,srci) = dipvec(cnt,1,srci) - m(ndi,3,srci)
            dipvec(cnt,2,srci) = dipvec(cnt,2,srci) + 0d0
            dipvec(cnt,3,srci) = dipvec(cnt,3,srci) + m(ndi,1,srci)
          endif
          if(dimi .eq. 3) then
            dipvec(cnt,1,srci) = dipvec(cnt,1,srci) + m(ndi,2,srci)
            dipvec(cnt,2,srci) = dipvec(cnt,2,srci) - m(ndi,1,srci)
            dipvec(cnt,3,srci) = dipvec(cnt,3,srci) + 0d0
          endif
          cnt=cnt+1
          enddo
          enddo 
          enddo
        endif
        !call hfmm3d direct
        ifpgh=0
        ifpghtarg=0
        if(ifeh.eq.1 .or. ifeh.eq.5) then
          ifpgh=1
        endif
        if(ifeh.eq.2 .or. ifeh.eq.6) then
          ifpgh=2
        endif
        if(ifehtarg.eq.1 .or. ifehtarg.eq.5) then
          ifpghtarg=1
        endif
        if(ifehtarg.eq.2 .or. ifehtarg.eq.6) then
          ifpghtarg=2
        endif
        !call helmholtz direct evaluations for source to source
        if(ifpgh.eq.1 .and. nsource.gt.0) then
          if(ifcharge.eq.1 .and. ifdipole.eq.1) then
            call h3ddirectcdp(nd*3, zk, source, charge, dipvec, nsource,
     $                        source, nsource, cpot, thresh)
          endif
          if(ifcharge.eq.1 .and. ifdipole.eq.0) then
            call h3ddirectcp(nd*3, zk, source, charge, nsource,
     $                       source, nsource, cpot, thresh)
          endif
          if(ifcharge.eq.0 .and. ifdipole.eq.1) then
            call h3ddirectdp(nd*3, zk, source, dipvec, nsource,
     $                       source, nsource, cpot, thresh)
          endif
        endif
        if(ifpgh.eq.2 .and. nsource.gt.0) then
          if(ifcharge.eq.1 .and. ifdipole.eq.1) then
            call h3ddirectcdg(nd*3, zk, source, charge, dipvec, nsource,
     $                        source, nsource, cpot, cgrad, thresh)
          endif
          if(ifcharge.eq.1 .and. ifdipole.eq.0) then
            call h3ddirectcg(nd*3, zk, source, charge, nsource,
     $                       source, nsource, cpot, cgrad, thresh)
          endif
          if(ifcharge.eq.0 .and. ifdipole.eq.1) then
            call h3ddirectdg(nd*3, zk, source, dipvec, nsource,
     $                       source, nsource, cpot, cgrad, thresh)
          endif
        endif
        !call helmholtz direct evaluations for source to target
        if(ifpghtarg.eq.1 .and. nsource.gt.0 .and. ntarg.gt.0) then
          if(ifcharge.eq.1 .and. ifdipole.eq.1) then
            call h3ddirectcdp(nd*3, zk, source, charge, dipvec, nsource,
     $                        targ, ntarg, cpottarg, thresh)
          endif
          if(ifcharge.eq.1 .and. ifdipole.eq.0) then
            call h3ddirectcp(nd*3, zk, source, charge, nsource,
     $                       targ, ntarg, cpottarg, thresh)
          endif
          if(ifcharge.eq.0 .and. ifdipole.eq.1) then
            call h3ddirectdp(nd*3, zk, source, dipvec, nsource,
     $                       targ, ntarg, cpottarg, thresh)
          endif
        endif
        if(ifpghtarg.eq.2 .and. nsource.gt.0 .and. ntarg.gt.0) then
          if(ifcharge.eq.1 .and. ifdipole.eq.1) then
            call h3ddirectcdg(nd*3, zk, source, charge, dipvec, nsource,
     $                        targ, ntarg, cpottarg, cgradtarg, thresh)
          endif
          if(ifcharge.eq.1 .and. ifdipole.eq.0) then
            call h3ddirectcg(nd*3, zk, source, charge, nsource,
     $                       targ, ntarg, cpottarg, cgradtarg, thresh)
          endif
          if(ifcharge.eq.0 .and. ifdipole.eq.1) then
            call h3ddirectdg(nd*3, zk, source, dipvec, nsource,
     $                       targ, ntarg, cpottarg, cgradtarg, thresh)
          endif
        endif
        !unpack from cpot  and cgrad
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,ndi,cnt)
        do srci=1,nsource
        cnt=1
        do dimi=1,3
        do ndi=1,nd
        if(ifeh.eq.1 .or. ifeh.eq.5 .or. ifeh.eq.2 .or. ifeh.eq.6) then
        e(ndi,dimi,srci) = cpot(cnt,srci)
        endif
        if(ifeh.eq.2 .or. ifeh.eq.6) then
        grade(ndi,dimi,1,srci) = cgrad(cnt,1,srci)
        grade(ndi,dimi,2,srci) = cgrad(cnt,2,srci)
        grade(ndi,dimi,3,srci) = cgrad(cnt,3,srci)
        endif
        cnt=cnt+1
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
        !unpack from cpottarg and cgradtarg
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(trgi,dimi,ndi,cnt)
        do trgi=1,ntarg
        cnt=1
        do dimi=1,3
        do ndi=1,nd
        if(ifehtarg.eq.1 .or. ifehtarg.eq.5 .or.
     $     ifehtarg.eq.2 .or. ifehtarg.eq.6) then
        etarg(ndi,dimi,trgi) = cpottarg(cnt,trgi)
        endif
        if(ifehtarg.eq.2 .or. ifehtarg.eq.6) then
        gradetarg(ndi,dimi,1,trgi) = cgradtarg(cnt,1,trgi)
        gradetarg(ndi,dimi,2,trgi) = cgradtarg(cnt,2,trgi)
        gradetarg(ndi,dimi,3,trgi) = cgradtarg(cnt,3,trgi)
        endif
        cnt=cnt+1
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
      endif

      dipvec=(0d0,0d0)
      charge=(0d0,0d0)
      cpot=(0d0,0d0)
      cgrad=(0d0,0d0)
      cpottarg=(0d0,0d0)
      cgradtarg=(0d0,0d0)
      !curl S_kJ
      if((ifeh.eq.3 .or. ifeh.eq.4 .or. ifeh.eq.5 .or. ifeh.eq.6) .or. 
     $   (ifehtarg.eq.3 .or. ifehtarg.eq.4 .or.
     $    ifehtarg.eq.5 .or. ifehtarg.eq.6)
     $  ) then
        ifcharge=0
        ifdipole=0
        !pack dipole and charge
        if(ifj .eq. 1) then
          ifdipole=1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,ndi,cnt,dimi)
          do srci=1,nsource
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          if(dimi .eq. 1) then
            dipvec(3*(ndi-1)+1,1,srci) = 0d0
            dipvec(3*(ndi-1)+1,2,srci) = -j(ndi,3,srci)
            dipvec(3*(ndi-1)+1,3,srci) = j(ndi,2,srci)
          endif
          if(dimi .eq. 2) then
            dipvec(3*(ndi-1)+2,1,srci) = j(ndi,3,srci)
            dipvec(3*(ndi-1)+2,2,srci) = 0d0
            dipvec(3*(ndi-1)+2,3,srci) = -j(ndi,1,srci)
          endif
          if(dimi .eq. 3) then
            dipvec(3*(ndi-1)+3,1,srci) = -j(ndi,2,srci)
            dipvec(3*(ndi-1)+3,2,srci) = j(ndi,1,srci)
            dipvec(3*(ndi-1)+3,3,srci) = 0d0
          endif
          cnt=cnt+1
          enddo
          enddo
          enddo
C$OMP END PARALLEL DO
        endif
        !continue pack dipole
        if(ifq .eq. 1) then
          ifdipole=1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,dimj,ndi,cnt)
          do srci=1,nsource
          do dimj=1,3
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          if(dimi .eq. dimj) then
            dipvec(cnt,dimj,srci) = dipvec(cnt,dimj,srci)+q(ndi,srci)
          endif
          cnt=cnt+1
          enddo
          enddo
          enddo
          enddo
C$OMP END PARALLEL DO
        endif
        !pack charge
        if(ifm .eq. 1) then
          ifcharge=1
          do srci=1,nsource
          cnt=1
          do dimi=1,3
          do ndi=1,nd
          charge(cnt,srci) = m(ndi,dimi,srci)*coef
          cnt=cnt+1
          enddo
          enddo
          enddo
        endif
        !call hfmm3d
        ifpgh=0
        ifpghtarg=0
        if(ifeh.eq.3 .or. ifeh.eq.5) then
          ifpgh=1
        endif
        if(ifeh.eq.4 .or. ifeh.eq.6) then
          ifpgh=2
        endif
        if(ifehtarg.eq.3 .or. ifehtarg.eq.5) then
          ifpghtarg=1
        endif
        if(ifehtarg.eq.4 .or. ifehtarg.eq.6) then
          ifpghtarg=2
        endif
        !call helmholtz direct evaluations for source to source
        if(ifpgh.eq.1 .and. nsource.gt.0) then
          if(ifcharge.eq.1 .and. ifdipole.eq.1) then
            call h3ddirectcdp(nd*3, zk, source, charge, dipvec, nsource,
     $                        source, nsource, cpot, thresh)
          endif
          if(ifcharge.eq.1 .and. ifdipole.eq.0) then
            call h3ddirectcp(nd*3, zk, source, charge, nsource,
     $                       source, nsource, cpot, thresh)
          endif
          if(ifcharge.eq.0 .and. ifdipole.eq.1) then
            call h3ddirectdp(nd*3, zk, source, dipvec, nsource,
     $                       source, nsource, cpot, thresh)
          endif
        endif
        if(ifpgh.eq.2 .and. nsource.gt.0) then
          if(ifcharge.eq.1 .and. ifdipole.eq.1) then
            call h3ddirectcdg(nd*3, zk, source, charge, dipvec, nsource,
     $                        source, nsource, cpot, cgrad, thresh)
          endif
          if(ifcharge.eq.1 .and. ifdipole.eq.0) then
            call h3ddirectcg(nd*3, zk, source, charge, nsource,
     $                       source, nsource, cpot, cgrad, thresh)
          endif
          if(ifcharge.eq.0 .and. ifdipole.eq.1) then
            call h3ddirectdg(nd*3, zk, source, dipvec, nsource,
     $                       source, nsource, cpot, cgrad, thresh)
          endif
        endif
        !call helmholtz direct evaluations for source to target
        if(ifpghtarg.eq.1 .and. nsource.gt.0 .and. ntarg.gt.0) then
          if(ifcharge.eq.1 .and. ifdipole.eq.1) then
            call h3ddirectcdp(nd*3, zk, source, charge, dipvec, nsource,
     $                        targ, ntarg, cpottarg, thresh)
          endif
          if(ifcharge.eq.1 .and. ifdipole.eq.0) then
            call h3ddirectcp(nd*3, zk, source, charge, nsource,
     $                       targ, ntarg, cpottarg, thresh)
          endif
          if(ifcharge.eq.0 .and. ifdipole.eq.1) then
            call h3ddirectdp(nd*3, zk, source, dipvec, nsource,
     $                       targ, ntarg, cpottarg, thresh)
          endif
        endif
        if(ifpghtarg.eq.2 .and. nsource.gt.0 .and. ntarg.gt.0) then
          if(ifcharge.eq.1 .and. ifdipole.eq.1) then
            call h3ddirectcdg(nd*3, zk, source, charge, dipvec, nsource,
     $                        targ, ntarg, cpottarg, cgradtarg, thresh)
          endif
          if(ifcharge.eq.1 .and. ifdipole.eq.0) then
            call h3ddirectcg(nd*3, zk, source, charge, nsource,
     $                       targ, ntarg, cpottarg, cgradtarg, thresh)
          endif
          if(ifcharge.eq.0 .and. ifdipole.eq.1) then
            call h3ddirectdg(nd*3, zk, source, dipvec, nsource,
     $                       targ, ntarg, cpottarg, cgradtarg, thresh)
          endif
        endif
        !unpack from cpot and cgrad
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(srci,dimi,ndi,cnt)
        do srci=1,nsource
        cnt=1
        do dimi=1,3
        do ndi=1,nd
        if(ifeh.eq.3 .or. ifeh.eq.5 .or. ifeh.eq.4 .or. ifeh.eq.6) then
        h(ndi,dimi,srci) = cpot(cnt,srci)
        endif
        if(ifeh.eq.4 .or. ifeh.eq.6) then
        gradh(ndi,dimi,1,srci) = cgrad(cnt,1,srci)
        gradh(ndi,dimi,2,srci) = cgrad(cnt,2,srci)
        gradh(ndi,dimi,3,srci) = cgrad(cnt,3,srci)
        endif
        cnt=cnt+1
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
        !unpack targ output from cpottarg and cgradtarg
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(trgi,dimi,ndi,cnt)
        do trgi=1,ntarg
        cnt=1
        do dimi=1,3
        do ndi=1,nd
        if(ifehtarg.eq.3 .or. ifehtarg.eq.5 .or.
     $     ifehtarg.eq.4 .or. ifehtarg.eq.6) then
        htarg(ndi,dimi,trgi) = cpottarg(cnt,trgi)
        endif
        if(ifehtarg.eq.4 .or. ifehtarg.eq.6) then
        gradhtarg(ndi,dimi,1,trgi) = cgradtarg(cnt,1,trgi)
        gradhtarg(ndi,dimi,2,trgi) = cgradtarg(cnt,2,trgi)
        gradhtarg(ndi,dimi,3,trgi) = cgradtarg(cnt,3,trgi)
        endif
        cnt=cnt+1
        enddo
        enddo
        enddo
C$OMP END PARALLEL DO
      endif

      return
      end
