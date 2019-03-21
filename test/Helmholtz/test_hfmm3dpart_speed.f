      implicit none
      integer ns,nt,icase
      double precision, allocatable :: source(:,:),targ(:,:)
      double complex, allocatable :: charge(:),dipstr(:)
      double precision, allocatable :: dipvec(:,:)
      double complex, allocatable :: pot(:),pottarg(:)
      double complex, allocatable :: grad(:,:),gradtarg(:,:)

      double precision eps
      double complex eye,zk,tt,ph,pi,done
      integer i,j,k,ntest
      integer ifcharge,ifdipole,ifpgh,ifpghtarg
      double precision err,hkrand
      

      data eye/(0.0d0,1.0d0)/

      done = 1.0d0
      pi = atan(done)*4

c
cc      initialize printing routine
c
      call prini(6,13)

      zk = 4.0d0 + eye*0.002d0

      ns = 300000
      nt = ns

      allocate(source(3,ns),targ(3,nt))
      allocate(charge(ns),dipstr(ns),dipvec(3,ns))
      allocate(pot(ns))
      allocate(grad(3,ns))

      allocate(pottarg(nt))
      allocate(gradtarg(3,nt))

     
      
      ntest = 10

      do icase = 1,10

        ns = 20000*icase
        nt = ns


c
cc      generate sources uniformly in the unit cube 
c
c
        do i=1,ns
          tt = pi*hkrand(0)
          ph = 2*pi*hkrand(0)
          source(1,i) = sin(tt)*cos(ph)
          source(2,i) = sin(tt)*sin(ph)
          source(3,i) = cos(tt)

          charge(i) = hkrand(0) + eye*hkrand(0)
          dipstr(i) = hkrand(0) + eye*hkrand(0)

          dipvec(1,i) = hkrand(0)
          dipvec(2,i) = hkrand(0)
          dipvec(3,i) = hkrand(0)

          pot(i) = 0
          grad(1,i) = 0
          grad(2,i) = 0
          grad(3,i) = 0
        enddo

c
cc      generate targets uniformly in the unit cube
c
        do i=1,nt
          targ(1,i) = hkrand(0)
          targ(2,i) = hkrand(0)
          targ(3,i) = hkrand(0)

          pottarg(i) = 0
          gradtarg(1,i) = 0
          gradtarg(2,i) = 0
          gradtarg(3,i) = 0 
        enddo

        eps = 0.5d-6
c
cc     now test source to source, charge, 
c      with potentials
c
         write(6,*) 'testing source to source'
         write(6,*) 'interaction: charges'
         write(6,*) 'output: potentials'
         write(6,*) 
         write(6,*) 

         call hfmm3dpartstoscp(eps,zk,ns,source,charge,
     1        pot)

         ifcharge = 1
         ifdipole = 0
         ifpgh = 1
         ifpghtarg = 0

       
         call comperr(zk,ns,source,ifcharge,charge,ifdipole,dipstr,
     1     dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2     ntest,err)

        call prin2('l2 rel err=*',err,1)
        write(6,*)
        write(6,*)
        write(6,*) '================'
      enddo


      stop
      end
c----------------------------------------------------------
c
cc
c
c
c
c
      subroutine comperr(zk,ns,source,ifcharge,charge,ifdipole,dipstr,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

      implicit none
      double complex zk
      integer ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg
      
      double precision source(3,*),targ(3,*)
      double precision dipvec(3,*)
      double complex dipstr(*),charge(*)

      double complex pot(*),pottarg(*),grad(3,*),gradtarg(3,*)

      integer i,j,ntest,nd

      double precision err,ra
      
      double complex potex(ntest),gradex(3,ntest),pottargex(ntest),
     1                  gradtargex(3,ntest)

      double precision thresh

      nd = 1
      err = 0 
      do i=1,ntest
        potex(i) = 0
        pottargex(i) = 0

        gradex(1,i) = 0
        gradex(2,i) = 0
        gradex(3,i) = 0

        gradtargex(1,i) = 0
        gradtargex(2,i) = 0
        gradtargex(3,i) = 0
      enddo

      thresh = 1.0d-16

      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        if(ifpgh.eq.1) then
          call h3ddirectcp(nd,zk,source,charge,ns,source,ntest,
     1       potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call h3ddirectcg(nd,zk,source,charge,ns,source,ntest,
     1       potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call h3ddirectcp(nd,zk,source,charge,ns,targ,ntest,
     1       pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call h3ddirectcg(nd,zk,source,charge,ns,targ,ntest,
     1       pottargex,gradtargex,thresh)
        endif
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        if(ifpgh.eq.1) then
          call h3ddirectdp(nd,zk,source,dipstr,dipvec,
     1      ns,source,ntest,potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call h3ddirectdg(nd,zk,source,dipstr,dipvec,
     1       ns,source,ntest,potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call h3ddirectdp(nd,zk,source,dipstr,dipvec,
     1       ns,targ,ntest,pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call h3ddirectdg(nd,zk,source,dipstr,dipvec,
     1       ns,targ,ntest,pottargex,gradtargex,thresh)
        endif
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        if(ifpgh.eq.1) then
          call h3ddirectcdp(nd,zk,source,charge,dipstr,dipvec,
     1      ns,source,ntest,potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call h3ddirectcdg(nd,zk,source,charge,dipstr,dipvec,
     1       ns,source,ntest,potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call h3ddirectcdp(nd,zk,source,charge,dipstr,dipvec,
     1       ns,targ,ntest,pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call h3ddirectcdg(nd,zk,source,charge,dipstr,dipvec,
     1       ns,targ,ntest,pottargex,gradtargex,thresh)
        endif
      endif

      err = 0
      ra = 0

      if(ifpgh.eq.1) then
        do i=1,ntest
          ra = ra + abs(potex(i))**2
          err = err + abs(pot(i)-potex(i))**2
        enddo
      endif

      if(ifpgh.eq.2) then
        do i=1,ntest
          ra = ra + abs(potex(i))**2
          ra = ra + abs(gradex(1,i))**2
          ra = ra + abs(gradex(2,i))**2
          ra = ra + abs(gradex(3,i))**2

          err = err + abs(pot(i)-potex(i))**2
          err = err + abs(grad(1,i)-gradex(1,i))**2
          err = err + abs(grad(2,i)-gradex(2,i))**2
          err = err + abs(grad(3,i)-gradex(3,i))**2
        enddo
      endif


      if(ifpghtarg.eq.1) then
        do i=1,ntest
          ra = ra + abs(pottargex(i))**2
          err = err + abs(pottarg(i)-pottargex(i))**2
        enddo
      endif

      if(ifpghtarg.eq.2) then
        do i=1,ntest
          ra = ra + abs(pottargex(i))**2
          ra = ra + abs(gradtargex(1,i))**2
          ra = ra + abs(gradtargex(2,i))**2
          ra = ra + abs(gradtargex(3,i))**2

          err = err + abs(pottarg(i)-pottargex(i))**2
          err = err + abs(gradtarg(1,i)-gradtargex(1,i))**2
          err = err + abs(gradtarg(2,i)-gradtargex(2,i))**2
          err = err + abs(gradtarg(3,i)-gradtargex(3,i))**2
        enddo
      endif

      err = sqrt(err/ra)
      return
      end
      
