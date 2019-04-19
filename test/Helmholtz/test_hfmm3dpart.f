      implicit none
      integer ns,nt
      double precision, allocatable :: source(:,:),targ(:,:)
      double complex, allocatable :: charge(:)
      double complex, allocatable :: dipvec(:,:)
      double complex, allocatable :: pot(:),pottarg(:)
      double complex, allocatable :: grad(:,:),gradtarg(:,:)

      double precision eps
      double complex eye,zk
      integer i,j,k,ntest
      integer ifcharge,ifdipole,ifpgh,ifpghtarg
      integer ipass(18),len1,ntests,isum
      double precision err,hkrand
      character(len=72) str1
      

      data eye/(0.0d0,1.0d0)/

c
cc      initialize printing routine
c
      call prini(6,13)

      zk = 2.2d0

      ns = 2000
      nt = 1999
      
      ntest = 10

      allocate(source(3,ns),targ(3,nt))
      allocate(charge(ns),dipvec(3,ns))
      allocate(pot(ns))
      allocate(grad(3,ns))

      allocate(pottarg(nt))
      allocate(gradtarg(3,nt))

      eps = 0.5d-9

      write(*,*) "=========================================="
      write(*,*) "Testing suite for hfmm3dpart"
      write(*,'(a,e11.5)') "Requested precision = ",eps

      open(unit=33,file='print_testres.txt',access='append')

      ntests = 18
      do i=1,ntests
        ipass(i) = 0
      enddo

c
cc      generate sources uniformly in the unit cube 
c
c
      do i=1,ns
        source(1,i) = hkrand(0)**2
        source(2,i) = hkrand(0)**2
        source(3,i) = hkrand(0)**2

        charge(i) = hkrand(0) + eye*hkrand(0)

        dipvec(1,i) = hkrand(0) + eye*hkrand(0) 
        dipvec(2,i) = hkrand(0) + eye*hkrand(0)  
        dipvec(3,i) = hkrand(0) + eye*hkrand(0)

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
     1      pot)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 1
       ifpghtarg = 0

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'

      if(err.lt.eps) ipass(1) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c

c
cc     now test source to source, charge, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstoscg(eps,zk,ns,source,charge,
     1      pot,grad)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(2) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 
      

c
cc     now test source to source, dipole, 
c      with potentials
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstosdp(eps,zk,ns,source,dipvec,
     1      pot)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 0

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(3) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 



c
cc     now test source to source, dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstosdg(eps,zk,ns,source,dipvec,
     1      pot,grad)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(4) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c
cc     now test source to source, charge + dipole, 
c      with potentials
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstoscdp(eps,zk,ns,source,charge,dipvec,
     1      pot)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 0 

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(5) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c
cc     now test source to source, charge + dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstoscdg(eps,zk,ns,source,charge,dipvec,
     1      pot,grad)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(6) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 



c
cc     now test source to target, charge, 
c      with potentials
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstotcp(eps,zk,ns,source,charge,
     1      nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(7) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c
cc     now test source to target, charge, 
c      with potentials and gradients
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstotcg(eps,zk,ns,source,charge,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(8) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 
      


c
cc     now test source to target, dipole, 
c      with potentials
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstotdp(eps,zk,ns,source,dipvec,
     1      nt,targ,pottarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(9) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c
cc     now test source to target, dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstotdg(eps,zk,ns,source,dipvec,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(10) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 

c
cc     now test source to target, charge + dipole, 
c      with potentials
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstotcdp(eps,zk,ns,source,charge,dipvec,
     1      nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(11) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c
cc     now test source to target, charge + dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstotcdg(eps,zk,ns,source,charge,dipvec,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(12) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 

c
cc     now test source to source + target, charge, 
c      with potentials
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstostcp(eps,zk,ns,source,charge,
     1      pot,nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(13) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c
cc     now test source to source + target, charge, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstostcg(eps,zk,ns,source,charge,
     1      pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(14) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 
      


c
cc     now test source to source + target, dipole, 
c      with potentials
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstostdp(eps,zk,ns,source,dipvec,
     1      pot,nt,targ,pottarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(15) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c
cc     now test source to source + target, dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstostdg(eps,zk,ns,source,dipvec,
     1      pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(16) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 

c
cc     now test source to source + target, charge + dipole, 
c      with potentials
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstostcdp(eps,zk,ns,source,charge,dipvec,
     1      pot,nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(17) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 


c
cc     now test source to source + target, charge + dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call hfmm3dpartstostcdg(eps,zk,ns,source,charge,dipvec,
     1      pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      if(err.lt.eps) ipass(18) = 1
      call gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(err.ge.eps) write(33,*) str1(1:len1) 

      isum = 0
      do i=1,ntests
        isum = isum+ipass(i)
      enddo

      write(*,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in hfmm3dpart testing suite'
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in hfmm3dpart testing suite'
      close(33)
      

      stop
      end
c----------------------------------------------------------
c
cc
c
c
c
c
      subroutine comperr(zk,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

      implicit none
      double complex zk
      integer ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg
      
      double precision source(3,*),targ(3,*)
      double complex dipvec(3,*)
      double complex charge(*)

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
          call h3ddirectdp(nd,zk,source,dipvec,
     1      ns,source,ntest,potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call h3ddirectdg(nd,zk,source,dipvec,
     1       ns,source,ntest,potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call h3ddirectdp(nd,zk,source,dipvec,
     1       ns,targ,ntest,pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call h3ddirectdg(nd,zk,source,dipvec,
     1       ns,targ,ntest,pottargex,gradtargex,thresh)
        endif
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        if(ifpgh.eq.1) then
          call h3ddirectcdp(nd,zk,source,charge,dipvec,
     1      ns,source,ntest,potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call h3ddirectcdg(nd,zk,source,charge,dipvec,
     1       ns,source,ntest,potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call h3ddirectcdp(nd,zk,source,charge,dipvec,
     1       ns,targ,ntest,pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call h3ddirectcdg(nd,zk,source,charge,dipvec,
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
c
c
c
c
c-------------------------------------------------------
      subroutine gererrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      implicit real *8 (a-h,o-z)
      character(len=*) str1
      character(len=13) str2
      character(len=14) str3
      character(len=19) str4
      character(len=18) str5

      str2 = "Failed src to"
      len1 = 13
      if(ifpgh.gt.0.and.ifpghtarg.eq.0) then
        str3 = " src,"
        len1 = len1+5  
      endif
      if(ifpgh.eq.0.and.ifpghtarg.gt.0) then
        str3 = " targ,"
        len1 = len1+6
      endif
      if(ifpgh.gt.0.and.ifpghtarg.gt.0) then
        str3 = " src and targ,"
        len1 = len1+14
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        str4=" charge,"
        len1 = len1+8
      endif
      
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        str4=" dipole,"
        len1 = len1+8
      endif
      
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        str4=" charge and dipole,"
        len1 = len1+19
      endif

      if(ifpgh.eq.1.or.ifpghtarg.eq.1) then
        str5=" pot test"
        len1 = len1 + 9
      endif
      
      if(ifpgh.eq.2.or.ifpghtarg.eq.2) then
        str5=" pot and grad test"
        len1 = len1 + 18
      endif

      str1 = str2//trim(str3)//trim(str4)//trim(str5)

      return
      end
