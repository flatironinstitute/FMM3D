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
      real *8 ztrgs(3,2)
	  real *8 c0(3),c1(3),c2(3),c3(3)
      real *8 xnodes(2000),wts(2000)
      real *8, allocatable :: dc(:,:)
c
      real *8, allocatable :: pots(:,:),flds(:,:,:),hesss(:,:,:)
      real *8, allocatable :: opots(:,:),oflds(:,:,:),ohesss(:,:,:)
      complex *16, allocatable :: mpole1(:,:,:),mpole2(:,:,:)
      complex *16, allocatable :: locexp1(:,:,:),locexp2(:,:,:)
      real *8, allocatable :: charge(:,:),dipvec(:,:,:)
      real *8  errs(3,5),err_exp(5)
      real *8, allocatable :: scarray_loc(:),scarray_mp(:)
      real *8, allocatable :: scarray_loc2(:),scarray_mp2(:)
      real *8, allocatable :: scarray_loc3(:),scarray_mp3(:)
      real *8, allocatable :: wlege(:)
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

      lw7 = 200000
      ll = lw7
      allocate(wlege(ll))
      allocate(scarray_mp(ll),scarray_loc(ll))
      allocate(scarray_mp2(ll),scarray_loc2(ll))
      allocate(scarray_mp3(ll),scarray_loc3(ll))

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

      nd = 3

      call prini(6,13)

      allocate(charge(nd,10),dipvec(nd,3,10))

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
      do idim=1,nd
        charge(idim,1) = 1.0d0
        dipvec(idim,1,1) = hkrand(0)
        dipvec(idim,2,1) = hkrand(0)
        dipvec(idim,3,1) = hkrand(0)
      enddo


      ns = 1
      sources(1,2)=c0(1)+0.25d0
      sources(2,2)=c0(2)-0.25d0
      sources(3,2)=c0(3)-0.25d0
      do idim=1,nd
        charge(idim,2)= -1.0d0  
        dipvec(idim,1,2) = hkrand(0)
        dipvec(idim,2,2) = hkrand(0)
        dipvec(idim,3,2) = hkrand(0)
      enddo
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
      ztrgs(1,1)= ztrg(1)
      ztrgs(2,1)= ztrg(2)
      ztrgs(3,1)= ztrg(3)
      ztrgs(1,2)= ztrg(1)+0.05d0
      ztrgs(2,2)= ztrg(2)+0.05d0
      ztrgs(3,2)= ztrg(3)+0.05d0
      nt = 2

      allocate(opots(nd,nt),pots(nd,nt))
      allocate(oflds(nd,3,nt),flds(nd,3,nt))
      allocate(ohesss(nd,6,nt),hesss(nd,6,nt))

      ntest = 5
      do i=1,ntest
        ipass(i) = 0 
      enddo
c
c       direct calculation:
c
      do i=1,nt
        do idim=1,nd
          opots(idim,i) = 0
          oflds(idim,1,i) = 0
          oflds(idim,2,i) = 0
          oflds(idim,3,i) = 0
          ohesss(idim,1,i) = 0
          ohesss(idim,2,i) = 0
          ohesss(idim,3,i) = 0
          ohesss(idim,4,i) = 0
          ohesss(idim,5,i) = 0
          ohesss(idim,6,i) = 0
        enddo
      enddo

      thresh = 1.0d-15
      call l3ddirectcdh(nd,sources,charge,dipvec,ns,ztrgs,
     1    nt,opots,oflds,ohesss,thresh)

      eps = 0.5d-12
      call l3dterms(eps, nterms)
      call l3dterms(eps, nterms2)
      call l3dterms(eps, nterms3)

      nterms2 = nterms - 5
      nterms3 = nterms + 10

      nmax = max(nterms,nterms2)
      nmax = max(nmax,nterms3) 

      print *, "nmax=",nmax


      rconv1 = 1.0d0/sqrt(3.0d0)
      rconv2 = sqrt(3.0d0)/2.0d0
      rconv3 = 0.75d0


      nlege = 100
      call ylgndrfwini(nlege,wlege,lw7,lused7)

      call l3dmpevalhessdini(nterms,scarray_mp)
      call l3dtaevalhessdini(nterms,scarray_loc)

      call l3dmpevalhessdini(nterms2,scarray_mp2)
      call l3dtaevalhessdini(nterms2,scarray_loc2)

      call l3dmpevalhessdini(nterms3,scarray_mp3)
      call l3dtaevalhessdini(nterms3,scarray_loc3)
c
c
c       create multipole expansion:
c
      allocate(mpole1(nd,0:nterms,-nterms:nterms))

      call mpzero(nd,mpole1,nterms)
      call l3dformmpcd(nd,rscale1,sources,charge,dipvec,
     1       ns,c0,nterms,mpole1,wlege,nlege)

      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
          hesss(idim,1,i) = 0
          hesss(idim,2,i) = 0
          hesss(idim,3,i) = 0
          hesss(idim,4,i) = 0
          hesss(idim,5,i) = 0
          hesss(idim,6,i) = 0
        enddo
      enddo

      call l3dmpevalh(nd,rscale1,c0,mpole1,nterms,ztrgs,nt,pots,
     1      flds,hesss,thresh,scarray_mp)

      err_exp(1) = 10*max(rconv1**(nterms),eps)
      write(*,'(a,e11.4)') 'Testing formmp and mpeval, expected error='
     1   ,err_exp(1)
      call errprinth(nd,nt,pots,opots,flds,oflds,hesss,ohesss,
     1   errs(1,1))
      if(max(errs(1,1),errs(2,1),errs(3,1)).lt.err_exp(1)) ipass(1) = 1
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

      allocate(mpole2(nd,0:nterms2,-nterms2:nterms2))
      call mpzero(nd,mpole2,nterms2)
      call l3dmpmp(nd,rscale1,c0,mpole1,nterms,
     1       rscale2,c1,mpole2,nterms2,dc,nn)


      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
          hesss(idim,1,i) = 0
          hesss(idim,2,i) = 0
          hesss(idim,3,i) = 0
          hesss(idim,4,i) = 0
          hesss(idim,5,i) = 0
          hesss(idim,6,i) = 0
        enddo
      enddo

      call l3dmpevalh(nd,rscale2,c1,mpole2,nterms2,ztrgs,nt,pots,
     1      flds,hesss,thresh,scarray_mp2)
      nn2 = min(nterms,nterms2)
      err_exp(2) = 10*max(rconv2**(nn2),eps)
      write(*,'(a,e11.4)') 'Testing mpmp, expected error='
     1   ,err_exp(2)
      call errprinth(nd,nt,pots,opots,flds,oflds,hesss,ohesss,
     1   errs(1,2))

      if(max(errs(1,2),errs(2,2)).lt.err_exp(2)) ipass(2) = 1
      write(*,*)
      if(ipass(2).ne.1) write(33,*) 'Failed mpmp test' 
c
cc    convert to local
c

      allocate(locexp2(nd,0:nterms2,-nterms2:nterms2))
      call mpzero(nd,locexp2,nterms2)
      call l3dmploc(nd,rscale2,c1,mpole2,nterms2,
     1      rscale2,c2,locexp2,nterms2,dc,nn)


      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
          hesss(idim,1,i) = 0
          hesss(idim,2,i) = 0
          hesss(idim,3,i) = 0
          hesss(idim,4,i) = 0
          hesss(idim,5,i) = 0
          hesss(idim,6,i) = 0
        enddo
      enddo

      call l3dtaevalh(nd,rscale2,c2,locexp2,nterms2,ztrgs,nt,
     1      pots,flds,hesss,scarray_loc2)
      nn2 = min(nterms,nterms2)

      err_exp(3) = 10*max(rconv3**(nn2),eps)
      write(*,'(a,e11.4)') 'Testing mploc, expected error='
     1   ,err_exp(3)
      call errprinth(nd,nt,pots,opots,flds,oflds,hesss,ohesss,
     1   errs(1,3))
      if(max(errs(1,3),errs(2,3)).lt.err_exp(3)) ipass(3) = 1
      write(*,*)
      if(ipass(3).ne.1) write(33,*) 'Failed mploc test' 



c
cc    shift local 
      allocate(locexp1(nd,0:nterms3,-nterms3:nterms3))
      call mpzero(nd,locexp1,nterms3)

      call l3dlocloc(nd,rscale2,c2,locexp2,nterms2,
     1      rscale1,c3,locexp1,nterms3,dc,nn)


      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
          hesss(idim,1,i) = 0
          hesss(idim,2,i) = 0
          hesss(idim,3,i) = 0
          hesss(idim,4,i) = 0
          hesss(idim,5,i) = 0
          hesss(idim,6,i) = 0
        enddo
      enddo

      call l3dtaevalh(nd,rscale1,c3,locexp1,nterms3,ztrgs,nt,
     1      pots,flds,hesss,scarray_loc3)

      nn2 = min(nterms,nterms3)
      err_exp(4) = 10*max(rconv3**(nn2),eps)
      write(*,'(a,e11.4)') 'Testing locloc, expected error='
     1   ,err_exp(4)

      call errprinth(nd,nt,pots,opots,flds,oflds,hesss,ohesss,
     1   errs(1,4))

      if(max(errs(1,4),errs(2,4)).lt.err_exp(4)) ipass(4) = 1
      write(*,*)
      if(ipass(4).ne.1) write(33,*) 'Failed locloc test' 

c
cc
cc
c    create local exp from sources

      call mpzero(nd,locexp1,nterms3)
      call l3dformtacd(nd,rscale1,sources,charge,dipvec,
     1       ns,c3,nterms3,locexp1,wlege,nlege)


      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
          hesss(idim,1,i) = 0
          hesss(idim,2,i) = 0
          hesss(idim,3,i) = 0
          hesss(idim,4,i) = 0
          hesss(idim,5,i) = 0
          hesss(idim,6,i) = 0
        enddo
      enddo

      call l3dtaevalh(nd,rscale1,c3,locexp1,nterms3,ztrgs,nt,
     1      pots,flds,hesss,scarray_loc3)
      nn = min(nterms,nterms3)
      err_exp(5) = 10*max(rconv3**(nn),eps)
      write(*,'(a,e11.4)') 'Testing formta and taeval, expected error='
     1   ,err_exp(5)
      call errprinth(nd,nt,pots,opots,flds,oflds,hesss,ohesss,
     1   errs(1,5))
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
c
c
      subroutine errprinth(nd,nt,pot,opot,fld,ofld,hess,ohess,errs)
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,nt),opot(nd,nt),fld(nd,3,nt)
      real *8 ofld(nd,3,nt),hess(nd,6,nt),ohess(nd,6,nt)
      real *8 errs(3)
 1000  format(4D15.5) 
      
      errp = 0
      dddp = 0
      errf = 0
      dddf = 0
      errh = 0
      dddh = 0

      do i=1,nt
        do idim=1,nd
          errp = errp + abs(pot(idim,i)-opot(idim,i))**2
          errf = errf + abs(fld(idim,1,i)-ofld(idim,1,i))**2
          errf = errf + abs(fld(idim,2,i)-ofld(idim,2,i))**2
          errf = errf + abs(fld(idim,3,i)-ofld(idim,3,i))**2

          do j=1,6
            errh = errh + abs(hess(idim,j,i)-ohess(idim,j,i))**2
            dddh = dddh + abs(ohess(idim,j,i))**2
          enddo

          dddp = dddp + abs(opot(idim,i))**2
          dddf = dddf + abs(ofld(idim,1,i))**2
          dddf = dddf + abs(ofld(idim,2,i))**2
          dddf = dddf + abs(ofld(idim,3,i))**2
        enddo
      enddo

c

      errs(1) = sqrt(errp/dddp)
      errs(2) = sqrt(errf/dddf)
      errs(3) = sqrt(errh/dddh)
      write(*,'(a,e11.4,a,e11.4,a,e11.4)') 
     1     'pot error=',errs(1),'   grad error=',errs(2),
     2     '   hess error=',errs(3)   
      return
      end
c
