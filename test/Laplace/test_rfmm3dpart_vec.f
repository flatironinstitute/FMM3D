      implicit none
      integer ns,nt,nd
      double precision, allocatable :: source(:,:),targ(:,:)
      double precision, allocatable :: charge(:,:)
      double precision, allocatable :: dipvec(:,:,:)
      double precision, allocatable :: pot(:,:),pottarg(:,:)
      double precision, allocatable :: grad(:,:,:),gradtarg(:,:,:)

      double precision eps
      double complex eye
      integer i,j,k,ntest,idim
      integer ifcharge,ifdipole,ifpgh,ifpghtarg
      double precision err,hkrand
      integer ipass(18),len1,ntests,isum
      character(len=72) str1
      

      data eye/(0.0d0,1.0d0)/

c
cc      initialize printing routine
c
      call prini(6,13)

      ns = 2000
      nt = 1999

      ntest = 10

      nd = 3

      allocate(source(3,ns),targ(3,nt))
      allocate(charge(nd,ns),dipvec(nd,3,ns))
      allocate(pot(nd,ns))
      allocate(grad(nd,3,ns))

      allocate(pottarg(nd,nt))
      allocate(gradtarg(nd,3,nt))


      eps = 0.5d-9

      write(*,*) "=========================================="
      write(*,*) "Testing suite for hfmm3dpart_vec"
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

        do idim=1,nd
          charge(idim,i) = hkrand(0) 

          dipvec(idim,1,i) = hkrand(0) 
          dipvec(idim,2,i) = hkrand(0)
          dipvec(idim,3,i) = hkrand(0)

          pot(idim,i) = 0
          grad(idim,1,i) = 0
          grad(idim,2,i) = 0
          grad(idim,3,i) = 0
        enddo
      enddo


c
cc      generate targets uniformly in the unit cube
c
      do i=1,nt
        targ(1,i) = hkrand(0)
        targ(2,i) = hkrand(0)
        targ(3,i) = hkrand(0)

        do idim=1,nd
          pottarg(idim,i) = 0
          gradtarg(idim,1,i) = 0
          gradtarg(idim,2,i) = 0
          gradtarg(idim,3,i) = 0 
        enddo
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

       call rfmm3dpartstoscp_vec(nd,eps,ns,source,charge,pot)


       ifcharge = 1
       ifdipole = 0
       ifpgh = 1
       ifpghtarg = 0

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstoscg_vec(nd,eps,ns,source,charge,
     1      pot,grad)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstosdp_vec(nd,eps,ns,source,dipvec,
     1      pot)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 0

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstosdg_vec(nd,eps,ns,source,dipvec,
     1      pot,grad)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstoscdp_vec(nd,eps,ns,source,charge,dipvec,
     1      pot)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 0 

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstoscdg_vec(nd,eps,ns,source,charge,dipvec,
     1      pot,grad)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstotcp_vec(nd,eps,ns,source,charge,
     1      nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstotcg_vec(nd,eps,ns,source,charge,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstotdp_vec(nd,eps,ns,source,dipvec,
     1      nt,targ,pottarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstotdg_vec(nd,eps,ns,source,dipvec,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstotcdp_vec(nd,eps,ns,source,charge,dipvec,
     1      nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstotcdg_vec(nd,eps,ns,source,charge,dipvec,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstostcp_vec(nd,eps,ns,source,charge,
     1      pot,nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstostcg_vec(nd,eps,ns,source,charge,
     1      pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstostdp_vec(nd,eps,ns,source,dipvec,
     1      pot,nt,targ,pottarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstostdg_vec(nd,eps,ns,source,dipvec,
     1      pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstostcdp_vec(nd,eps,ns,source,charge,
     1       dipvec,pot,nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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

       call rfmm3dpartstostcdg_vec(nd,eps,ns,source,charge,
     1    dipvec,pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
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
     1   ' out of ',ntests,' tests in rfmm3dpart vec testing suite'
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in rfmm3dpart vec testing suite'
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
      subroutine comperr_vec(nd,ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,
     2   gradtarg,ntest,err)

      implicit none
      double complex zk
      integer ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg
      
      double precision source(3,*),targ(3,*)
      double precision dipvec(nd,3,*)
      double precision charge(nd,*)

      double precision pot(nd,*),pottarg(nd,*),grad(nd,3,*),
     1   gradtarg(nd,3,*)

      integer i,j,ntest,nd,idim

      double precision err,ra
      
      double precision potex(nd,ntest),gradex(nd,3,ntest),
     1     pottargex(nd,ntest),gradtargex(nd,3,ntest)

      double precision thresh

      err = 0 
      do i=1,ntest
        do idim=1,nd
          potex(idim,i) = 0
          pottargex(idim,i) = 0

          gradex(idim,1,i) = 0
          gradex(idim,2,i) = 0
          gradex(idim,3,i) = 0

          gradtargex(idim,1,i) = 0
          gradtargex(idim,2,i) = 0
          gradtargex(idim,3,i) = 0
        enddo
      enddo

      thresh = 1.0d-16

      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        if(ifpgh.eq.1) then
          call l3ddirectcp(nd,source,charge,ns,source,ntest,
     1       potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call l3ddirectcg(nd,source,charge,ns,source,ntest,
     1       potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call l3ddirectcp(nd,source,charge,ns,targ,ntest,
     1       pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call l3ddirectcg(nd,source,charge,ns,targ,ntest,
     1       pottargex,gradtargex,thresh)
        endif
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        if(ifpgh.eq.1) then
          call l3ddirectdp(nd,source,dipvec,
     1      ns,source,ntest,potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call l3ddirectdg(nd,source,dipvec,
     1       ns,source,ntest,potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call l3ddirectdp(nd,source,dipvec,
     1       ns,targ,ntest,pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call l3ddirectdg(nd,source,dipvec,
     1       ns,targ,ntest,pottargex,gradtargex,thresh)
        endif
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        if(ifpgh.eq.1) then
          call l3ddirectcdp(nd,source,charge,dipvec,
     1      ns,source,ntest,potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call l3ddirectcdg(nd,source,charge,dipvec,
     1       ns,source,ntest,potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call l3ddirectcdp(nd,source,charge,dipvec,
     1       ns,targ,ntest,pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call l3ddirectcdg(nd,source,charge,dipvec,
     1       ns,targ,ntest,pottargex,gradtargex,thresh)
        endif
      endif

      err = 0
      ra = 0

      if(ifpgh.eq.1) then
        do i=1,ntest
          do idim=1,nd
            ra = ra + abs(potex(idim,i))**2
            err = err + abs(pot(idim,i)-potex(idim,i))**2
          enddo
        enddo
      endif

      if(ifpgh.eq.2) then
        do i=1,ntest
          do idim=1,nd
            ra = ra + abs(potex(idim,i))**2
            ra = ra + abs(gradex(idim,1,i))**2
            ra = ra + abs(gradex(idim,2,i))**2
            ra = ra + abs(gradex(idim,3,i))**2

            err = err + abs(pot(idim,i)-potex(idim,i))**2
            err = err + abs(grad(idim,1,i)-gradex(idim,1,i))**2
            err = err + abs(grad(idim,2,i)-gradex(idim,2,i))**2
            err = err + abs(grad(idim,3,i)-gradex(idim,3,i))**2
          enddo
        enddo
      endif


      if(ifpghtarg.eq.1) then
        do i=1,ntest
          do idim=1,nd
            ra = ra + abs(pottargex(idim,i))**2
            err = err + abs(pottarg(idim,i)-pottargex(idim,i))**2
          enddo
        enddo
      endif

      if(ifpghtarg.eq.2) then
        do i=1,ntest
          do idim=1,nd
            ra = ra + abs(pottargex(idim,i))**2
            ra = ra + abs(gradtargex(idim,1,i))**2
            ra = ra + abs(gradtargex(idim,2,i))**2
            ra = ra + abs(gradtargex(idim,3,i))**2

            err = err + abs(pottarg(idim,i)-pottargex(idim,i))**2
            err = err + abs(gradtarg(idim,1,i)-gradtargex(idim,1,i))**2
            err = err + abs(gradtarg(idim,2,i)-gradtargex(idim,2,i))**2
            err = err + abs(gradtarg(idim,3,i)-gradtargex(idim,3,i))**2
          enddo
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
