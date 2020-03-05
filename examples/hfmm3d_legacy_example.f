cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is the third release of the FMM3D library, together with
c     associated subroutines, which computes N-body interactions
c     governed by the Laplace or Helmholtz equations.
c  
c
      program testhelm
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: source(:,:),dipvec(:,:)
      complex *16, allocatable :: charge(:),dipstr(:)
      complex *16, allocatable :: pot(:),fld(:,:)

      complex *16, allocatable :: pot2(:),fld2(:,:)
        
      real *8, allocatable :: targ(:,:)
      complex *16, allocatable :: pottarg(:),fldtarg(:,:)
      complex *16, allocatable :: pottarg2(:),fldtarg2(:,:)
        
      complex *16 ptemp,ftemp(3)
c       
      complex *16 ima
      complex *16 zk
      data ima/(0.0d0,1.0d0)/
c
      done=1
      pi=4*atan(done)

      write(*,*) "============================="
      write(*,*) "Example fortran driver for legacy drivers"
      write(*,*) "On output tests potential at a few sources and targs"
c
c     Initialize simple printing routines. The parameters to prini
c     define output file numbers using standard Fortran conventions.
c
c     Calling prini(6,13) causes printing to the screen and to 
c     file fort.13.     
c
      call prini(6,13)
c
      nsource= 10000

      allocate(source(3,nsource),charge(nsource))
      allocate(dipstr(nsource),dipvec(3,nsource))
      allocate(pot(nsource),fld(3,nsource))
c
c     set Helmholtz parameter  (used only by H3D prefix routines)
c
      zk = 1.0d0 + ima*0.1d0
c
c     construct randomly located charge distribution on a unit sphere
c 
      d=hkrand(0)
      do i=1,nsource
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        source(1,i)=.5d0*cos(phi)*sin(theta)
        source(2,i)=.5d0*sin(phi)*sin(theta)
        source(3,i)=.5d0*cos(theta)
      enddo
c
c     construct target distribution on a target unit sphere 
c
      ntarg=nsource
      allocate(targ(3,ntarg),pottarg(ntarg),fldtarg(3,ntarg))
      do i=1,ntarg
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        targ(1,i)=.5d0*cos(phi)*sin(theta) 
        targ(2,i)=.5d0*sin(phi)*sin(theta)
        targ(3,i)=.5d0*cos(theta)
      enddo
c
c       
c     set precision flag
c
      iprec=1
c       
c     set source type flags and output flags
c
      ifpot=1
      iffld=0
c
      ifcharge=1
      ifdipole=1
c
      ifpottarg=1
      iffldtarg=0
c
c       set source strengths
c
      if (ifcharge .eq. 1 ) then
        do i=1,nsource
          charge(i)=hkrand(0) + ima*hkrand(0)
        enddo
      endif
c
      if (ifdipole .eq. 1) then
        do i=1,nsource
          dipstr(i)=hkrand(0) + ima*hkrand(0)
          dipvec(1,i)=hkrand(0)
          dipvec(2,i)=hkrand(0)
          dipvec(3,i)=hkrand(0)
        enddo
      endif

c
c     initialize timing call
c
        call cpu_time(t1)
C$        t1=omp_get_wtime()
c       
c     call FMM3D routine for sources and targets
c
        call hfmm3dparttarg(ier,iprec, zk,
     $     nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ntarg,targ,
     $     ifpottarg,pottarg,iffldtarg,fldtarg)

c       
c     get time for FMM call
c
        call cpu_time(t2)
C$        t2=omp_get_wtime()
c       
c       
      call prinf('nsource=*',nsource,1)
      call prinf('ntarg=*',ntarg,1)
      call prin2('after fmm, time (sec)=*',t2-t1,1)
      call prin2('after fmm, speed (points+targets/sec)=*',
     $     (nsource+ntarg)/(t2-t1),1)
c       
c     call direct calculation with subset of points to assess accuracy
c
      m=min(nsource,100)
c
c     ifprint=0 suppresses printing of source locations
c     ifprint=1 turns on printing of source locations
c
      ifprint=0
      if (ifprint .eq. 1) then
        call prin2('source=*',source,3*6)
      endif

      if (ifprint .eq. 1) then
        if( ifpot.eq.1 ) call prin2('after fmm, pot=*',pot,2*m)
        if( iffld.eq.1 ) call prin2('after fmm, fld=*',fld,3*2*m)
      endif

c
c       for direct calculation, initialize pot2,fld2 arrays to zero.
c
      allocate(pot2(m),fld2(3,m))
      do i=1,m
        if (ifpot .eq. 1) pot2(i)=0
        if (iffld .eq. 1) then
          fld2(1,i)=0
          fld2(2,i)=0
          fld2(3,i)=0
        endif
      enddo


      ifpottmp = 0
      iffldtmp = 0

        call cpu_time(t1)
C$        t1=omp_get_wtime()
     
      
      call prinf('m=*',m,1)
      call prinf('nsource=*',nsource,1)
      call prin2('source=*',source,12)
      call prinf('ifcharge=*',ifcharge,1)
      call prinf('ifdipole=*',ifdipole,1)
      call prinf('m=*',m,1)
      call prinf('ifpot=*',ifpot,1)
      call prinf('iffld=*',iffld,1)
      call h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,
     1       dipstr,dipvec,ifpottmp,pot2,iffldtmp,fld2,m,source,
     2       ifpot,pot2,iffld,fld2)

      
        call cpu_time(t2)
C$        t2=omp_get_wtime()
c
c       ifprint=1 turns on printing of first m values of potential and field
c
        if (ifprint .eq. 1) then
           if( ifpot.eq.1 ) call prin2('directly, pot=*',pot2,2*m)
           if( iffld.eq.1 ) call prin2('directly, fld=*',fld2,3*2*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(nsource)/dble(m),1)
        call prin2('directly, estimated speed (points/sec)=*',
     $     m/(t2-t1),1)
c       
        if (ifpot .eq. 1)  then
           call h3derror(pot,pot2,m,aerr,rerr)
           call prin2('relative L2 error in potential=*',rerr,1)
        endif
c
        if (iffld .eq. 1) then
           call h3derror(fld,fld2,3*m,aerr,rerr)
           call prin2('relative L2 error in field=*',rerr,1)
        endif

c
c
c          now evaluate potential directly at targets
c

      m=min(ntarg,100)
      allocate(pottarg2(m),fldtarg2(3,m))

      do i=1,m
        if (ifpottarg .eq. 1) pottarg2(i)=0
        if (iffldtarg .eq. 1) then
          fldtarg2(1,i)=0
          fldtarg2(2,i)=0
          fldtarg2(3,i)=0
        endif

      enddo




        call cpu_time(t1)
C$        t1=omp_get_wtime()
c
      call h3dpartdirect(zk,nsource,source,ifcharge,charge,ifdipole,
     1       dipstr,dipvec,ifpottmp,pot2,iffldtmp,fld2,m,targ,
     2       ifpottarg,pottarg2,iffldtarg,fldtarg2)
        call cpu_time(t2)
C$        t2=omp_get_wtime()
c
c
c
        if (ifprint .eq. 1) then
           if( ifpottarg.eq.1 ) 
     $        call prin2('after fmm, pottarg=*',pottarg,2*m)
           if( iffldtarg.eq.1 ) 
     $        call prin2('after fmm, fldtarg=*',fldtarg,3*2*m)
        endif

        if (ifprint .eq. 1) then
           if (ifpottarg .eq. 1) 
     $        call prin2('directly, pottarg=*',pottarg2,2*m)
           if( iffldtarg.eq.1 ) 
     $        call prin2('directly, fldtarg=*',fldtarg2,3*2*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(ntarg)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',
     $     m/(t2-t1),1)
c       
        if (ifpottarg .eq. 1) then
           call h3derror(pottarg,pottarg2,m,aerr,rerr)
           call prin2('relative L2 error in target potential=*',rerr,1)
        endif
c
        if (iffldtarg .eq. 1) then
           call h3derror(fldtarg,fldtarg2,3*m,aerr,rerr)
           call prin2('relative L2 error in target field=*',rerr,1)
        endif
c       
        stop
        end
c
c
c
c
        subroutine h3derror(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        complex *16 pot1(n),pot2(n)
c
        d=0
        a=0
c       
        do i=1,n
           d=d+abs(pot1(i)-pot2(i))**2
           a=a+abs(pot1(i))**2
        enddo
c       
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c       
        ae=d
        re=d/a
c       
        return
        end
c
c
c

