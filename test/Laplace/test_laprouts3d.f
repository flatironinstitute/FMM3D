c
c      TESTING SUITE FOR Laplace SUBROUTINE LIBRARY:  Mar 20, 2019
c
c
c
c         Sources at S, Targets at T.
c
c         l3dformmp      forms expansion about c0
c         l3dmpeval      evaluates expansion at T
c         l3dmpmp        shifts expansion to c1
c         l3dmploc       converts to local about y1
c         l3dtaeval      evaluates expansion at T
c         l3dlocloc      shifts to local about c0
c         l3dtaeval      evaluates expansion at T
c         l3dformta      forms local expansion about c0 
c                        directly from S
c         l3dtaeval      evaluates expansion at T
c
c         -------------------------------------------------
c         |       |   S   |       |       |T      |       |
c         |       |   c0  |       |       |  c3   |       |
c         |       |       |       |       |       |       |
c         -------c1-------------------------------c2------
c         |       |       |       |       |       |       |
c         |       |       |       |       |       |       |
c         |       |       |       |       |       |       |
c         -------------------------------------------------
c
      implicit real *8 (a-h,o-z)
      real *8 ztrg(3),sources(3,10)
	  real *8 c0(3),c1(3),c2(3),c3(3)
      real *8 xnodes(2000),wts(2000)
      real *8, allocatable :: dc(:,:)
c
      real *8 pot,fld(3),opot,ofld(3)
      complex *16, allocatable :: mpole1(:,:),mpole2(:,:)
      complex *16, allocatable :: locexp1(:,:),locexp2(:,:)
      real *8 charge(100)
      real *8 dipvec(3,100),wlege(100000),errs(2,5),err_exp(5)
      integer ipass(5)
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1
      pi=4.0*atan(done)
      call prini(6,13)
      write(*,*) "=========================================="
      write(*,*) "Testing suite for laprouts3d"

      open(unit=33,file='print_testres.txt',access='append')

c
c
c       set scale parameter neq 1 to make sure it is used correctly.
c
      bsize = 1.0d0
      rscale1 = bsize/4
      rscale2 = bsize/2

      rscale1 = hkrand(0)
      rscale2 = rscale1*2

cc      rscale1 = 1
cc      rscale2 = 1
c
c       set center for original multipole expansion
c
      c0(1)=0.248d0
      c0(2)=0.251d0
      c0(3)=0.249d0
c
c       create two charge sources.
c
      sources(1,1)=c0(1)+0.24d0
      sources(2,1)=c0(2)+ 0.25d0
      sources(3,1)=c0(3)+ 0.25d0
      charge(1)= 1.0d0 
      dipvec(1,1) = hkrand(0)
      dipvec(2,1) = hkrand(0)
      dipvec(3,1) = hkrand(0)


      ns = 1
      sources(1,2)=c0(1)+0.25d0
      sources(2,2)=c0(2)-0.25d0
      sources(3,2)=c0(3)-0.25d0
      charge(2)= -1.0d0  
      dipvec(1,2) = hkrand(0)
      dipvec(2,2) = hkrand(0)
      dipvec(3,2) = hkrand(0)
      ns = 2
c
c       create center for shifted expansions
c       (location jiggled for good measure)
c
      c1(1)= 0.015d0
      c1(2)= 0.012d0
      c1(3)= 0.013d0
      c2(1)= 2.015d0
      c2(2)= 0.012d0
      c2(3)= 0.013d0
      c3(1)= c2(1)-0.25d0
      c3(2)= c2(2)-0.249d0
      c3(3)= c2(3)-0.251d0
c
c       create target
c
      ztrg(1)=c3(1)-0.245d0
      ztrg(2)=c3(2)-0.25d0
      ztrg(3)=c3(3)-0.25d0
      nt = 1

      ntest = 5
      do i=1,ntest
        ipass(i) = 0 
      enddo
c
c       direct calculation:
c
      opot = 0
      ofld(1) = 0
      ofld(2) = 0
      ofld(3) = 0
      thresh = 1.0d-15
      nd = 1
      call l3ddirectcdg(nd,sources,charge,dipvec,ns,ztrg,
     1    nt,opot,ofld,thresh)

      eps = 0.5d-12
      call l3dterms(eps, nterms)
      call l3dterms(eps, nterms2)
      call l3dterms(eps, nterms3)

      rconv1 = 1.0d0/sqrt(3.0d0)
      rconv2 = sqrt(3.0d0)/2.0d0
      rconv3 = 0.75d0


      nlege = 100
      lw7 = 100000
      call ylgndrfwini(nlege,wlege,lw7,lused7)
c
c
c       create multipole expansion:
c
      allocate(mpole1(0:nterms,-nterms:nterms))

      nd = 1
      call mpzero(nd,mpole1,nterms)
      call l3dformmpcd(nd,rscale1,sources,charge,dipvec,
     1       ns,c0,nterms,mpole1,wlege,nlege)


      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      nd = 1

      call l3dmpevalg(nd,rscale1,c0,mpole1,nterms,ztrg,nt,pot,
     1      fld,wlege,nlege,thresh)

      err_exp(1) = 10*max(rconv1**(nterms),eps)
      write(*,'(a,e11.4)') 'Testing formmp and mpeval, expected error='
     1   ,err_exp(1)
      call errprint(pot,opot,fld,ofld,errs(1,1))
      if(max(errs(1,1),errs(2,1)).lt.err_exp(1)) ipass(1) = 1
      write(*,*)

      if(ipass(1).ne.1) write(33,*) 'Failed formmp and mpeval test' 




c
c    mpmp shift
c


      nn = nterms
      nn = max(nn,nterms2)
      nn = max(nn,nterms3)
      nn = 2*nn + 10
      allocate(dc(0:nn,0:nn))
      call getsqrtbinomialcoeffs(nn,dc)

      allocate(mpole2(0:nterms2,-nterms2:nterms2))
      call mpzero(nd,mpole2,nterms2)
      call l3dmpmp(nd,rscale1,c0,mpole1,nterms,
     1       rscale2,c1,mpole2,nterms2,dc,nn)

      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      nd = 1
      call l3dmpevalg(nd,rscale2,c1,mpole2,nterms,ztrg,nt,pot,
     1      fld,wlege,nlege,thresh)
      err_exp(2) = 10*max(rconv2**(nterms),eps)
      write(*,'(a,e11.4)') 'Testing mpmp, expected error='
     1   ,err_exp(2)
      call errprint(pot,opot,fld,ofld,errs(1,2))

      if(max(errs(1,2),errs(2,2)).lt.err_exp(2)) ipass(2) = 1
      write(*,*)
      if(ipass(2).ne.1) write(33,*) 'Failed mpmp test' 
c
cc    convert to local
c

      allocate(locexp2(0:nterms2,-nterms2:nterms2))
      call mpzero(nd,locexp2,nterms2)
      call l3dmploc(nd,rscale2,c1,mpole2,nterms2,
     1      rscale2,c2,locexp2,nterms2,dc,nn)


      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      nd = 1
      call l3dtaevalg(nd,rscale2,c2,locexp2,nterms2,ztrg,nt,
     1      pot,fld,wlege,nlege)
      err_exp(3) = 10*max(rconv3**(nterms),eps)
      write(*,'(a,e11.4)') 'Testing mploc, expected error='
     1   ,err_exp(3)
      call errprint(pot,opot,fld,ofld,errs(1,3))
      if(max(errs(1,3),errs(2,3)).lt.err_exp(3)) ipass(3) = 1
      write(*,*)
      if(ipass(3).ne.1) write(33,*) 'Failed mploc test' 



c
cc    shift local 
      allocate(locexp1(0:nterms3,-nterms3:nterms3))
      nd = 1
      call mpzero(nd,locexp1,nterms3)

      call l3dlocloc(nd,rscale2,c2,locexp2,nterms2,
     1      rscale1,c3,locexp1,nterms3,dc,nn)


      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      call l3dtaevalg(nd,rscale1,c3,locexp1,nterms3,ztrg,nt,
     1      pot,fld,wlege,nlege)

      err_exp(4) = 10*max(rconv3**(nterms),eps)
      write(*,'(a,e11.4)') 'Testing locloc, expected error='
     1   ,err_exp(4)
      call errprint(pot,opot,fld,ofld,errs(1,4))
      if(max(errs(1,4),errs(2,4)).lt.err_exp(4)) ipass(4) = 1
      write(*,*)
      if(ipass(4).ne.1) write(33,*) 'Failed locloc test' 

c
cc
cc
c    create local exp from sources

      nd = 1
      call mpzero(nd,locexp1,nterms3)
      call l3dformtacd(nd,rscale1,sources,charge,dipvec,
     1       ns,c3,nterms3,locexp1,wlege,nlege)

      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      call l3dtaevalg(nd,rscale1,c3,locexp1,nterms3,ztrg,nt,
     1      pot,fld,wlege,nlege)

      err_exp(5) = 10*max(rconv3**(nterms),eps)
      write(*,'(a,e11.4)') 'Testing formta and taeval, expected error='
     1   ,err_exp(5)
      call errprint(pot,opot,fld,ofld,errs(1,5))
      if(max(errs(1,5),errs(2,5)).lt.err_exp(5)) ipass(5) = 1
      write(*,*)
      if(ipass(5).ne.1) write(33,*) 'Failed formta and taeval test'

      isum = 0
      do i=1,ntest
        isum = isum+ipass(i)
      enddo

      write(*,'(a,i1,a,i1,a)') 'Successfully completed ',isum,
     1   ' out of ',ntest,' tests in laprouts3d testing suite'
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',isum,
     1   ' out of ',ntest,' tests in laprouts3d testing suite'
      close(33)
      


      stop
      end
c
c
c
c
c
c
      subroutine errprint(pot,opot,fld,ofld,errs)
      implicit real *8 (a-h,o-z)
      real *8 pot,opot,fld(3),ofld(3)
      real *8 errs(2)
 1000  format(4D15.5) 
      err = 0
      ddd = 0
      err = err + abs(fld(1)-ofld(1))**2
      err = err + abs(fld(2)-ofld(2))**2
      err = err + abs(fld(3)-ofld(3))**2
      ddd = ddd + abs(ofld(1))**2
      ddd = ddd + abs(ofld(2))**2
      ddd = ddd + abs(ofld(3))**2
      err = sqrt(err)
      ddd = sqrt(ddd)

      err1 = abs(pot-opot)/abs(opot)
      write(*,'(a,e11.4,a,e11.4)') 
     1     'pot error=',err1,'   grad error=',err/ddd

      errs(1) = err1 
      errs(2) = err/ddd

      return
      end
c
