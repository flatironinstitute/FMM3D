c
c      TESTING SUITE FOR HELMHOLTZ SUBROUTINE LIBRARY:  Mar 20, 2019
c
c
c
c         Sources at S, Targets at T.
c
c         h3dformmp      forms expansion about c0
c         h3dmpeval      evaluates expansion at T
c         h3dmpmp        shifts expansion to c1
c         h3dmploc       converts to local about y1
c         h3dtaeval      evaluates expansion at T
c         h3dlocloc      shifts to local about c0
c         h3dtaeval      evaluates expansion at T
c         h3dformta      forms local expansion about c0 
c                        directly from S
c         h3dtaeval      evaluates expansion at T
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
c
      complex *16 pot,fld(3),opot,ofld(3)
      complex *16, allocatable :: mpole1(:,:),mpole2(:,:)
      complex *16, allocatable :: locexp1(:,:),locexp2(:,:)
      complex *16 charge(100),dipvec(3,100)
      real *8 wlege(100000)
      complex *16 zk,eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1
      pi=4.0*atan(done)
      zk = 0.2d0
c
      call prini(6,13)
c
c       set scale parameter neq 1 to make sure it is used correctly.
c
      bsize = 1.0d0
      rscale1 = abs(zk)*bsize/4
      rscale2 = abs(zk)*bsize/2

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
      charge(1)= 1.0d0 + eye*0.2d0
      dipvec(1,1) = hkrand(0) + eye*hkrand(0)
      dipvec(2,1) = hkrand(0) + eye*hkrand(0)
      dipvec(3,1) = hkrand(0) + eye*hkrand(0)


      ns = 1
      sources(1,2)=c0(1)+0.25d0
      sources(2,2)=c0(2)-0.25d0
      sources(3,2)=c0(3)-0.25d0
      charge(2)= -1.0d0 + eye*0.1d0
      dipvec(1,2) = hkrand(0) + eye*hkrand(0) 
      dipvec(2,2) = hkrand(0) + eye*hkrand(0)
      dipvec(3,2) = hkrand(0) + eye*hkrand(0)
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
c
c       direct calculation:
c
      opot = 0
      ofld(1) = 0
      ofld(2) = 0
      ofld(3) = 0
      thresh = 1.0d-15
      nd = 1
      call h3ddirectcdg(nd,zk,sources,charge,
     1       dipvec,ns,ztrg,nt,opot,
     1       ofld,thresh)



      eps = 0.5d-12
      call h3dterms(bsize/2, zk, eps, nterms)
      call h3dterms(bsize, zk, eps, nterms2)
      call h3dterms(bsize/2, zk, eps, nterms3)

      nterms = 90
      nterms2 = 90
      nterms3 = 90

      nlege = 100
      lw7 = 100000
      call ylgndrfwini(nlege,wlege,lw7,lused7)
c
c
c       create h-expansion:
c
      call prinf('calling formmp and mpeval =*', lused,0)
      ifinit = 1
      nquad = 2.0*max(nterms2,nterms3)
      call legewhts(nquad,xnodes,wts,ifinit)
      allocate(mpole1(0:nterms,-nterms:nterms))

      nd = 1
      call mpzero(nd,mpole1,nterms)
      call h3dformmpcd(nd,zk,rscale1,sources,charge,
     1       dipvec,ns,c0,nterms,mpole1,wlege,nlege)

      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      nd = 1

      call h3dmpevalg(nd,zk,rscale1,c0,mpole1,nterms,ztrg,nt,pot,
     1      fld,wlege,nlege,thresh)


      call errprint(pot,opot,fld,ofld)

c
c    mpmp shift
c
      call prinf('calling mpmp and mpeval *',nterms,0)
      radius = sqrt(3.0d0)/2.0d0
      allocate(mpole2(0:nterms2,-nterms2:nterms2))
      call mpzero(nd,mpole2,nterms2)
      call h3dmpmp(nd,zk,rscale1,c0,mpole1,nterms,
     1       rscale2,c1,mpole2,nterms2,radius,xnodes,wts,nquad)

      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      nd = 1
      call h3dmpevalg(nd,zk,rscale2,c1,mpole2,nterms,ztrg,nt,pot,
     1      fld,wlege,nlege,thresh)
       call errprint(pot,opot,fld,ofld)
c
c    convert to local
c
      radius = sqrt(3.0d0)/2.0d0
      call prin2('calling mploc and taeval *',wavek,0)
      nquad = 2.2*nterms2
      call legewhts(nquad,xnodes,wts,ifinit)

      allocate(locexp2(0:nterms2,-nterms2:nterms2))
      call mpzero(nd,locexp2,nterms2)
      call h3dmploc(nd,zk,rscale2,c1,mpole2,nterms2,
     1      rscale2,c2,locexp2,nterms2,radius,xnodes,wts,nquad)

cc      call prinm(locexp2,nterms2)

      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      nd = 1
      call h3dtaevalg(nd,zk,rscale2,c2,locexp2,nterms2,ztrg,nt,
     1      pot,fld,wlege,nlege)

      call errprint(pot,opot,fld,ofld)

c
c    shift local 
c
      radius = sqrt(3.0d0)/4.0d0
      call prinf('calling locloc and taeval *',nterms3,0)

      allocate(locexp1(0:nterms3,-nterms3:nterms3))
      call mpzero(nd,locexp1,nterms3)

      call h3dlocloc(nd,zk,rscale2,c2,locexp2,nterms2,
     1      rscale1,c3,locexp1,nterms3,radius,xnodes,wts,nquad)


      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      nd = 1
      call h3dtaevalg(nd,zk,rscale1,c3,locexp1,nterms3,ztrg,nt,
     1      pot,fld,wlege,nlege)

      call errprint(pot,opot,fld,ofld)

c
c
c    create local exp from sources

           call prin2('calling l3dformta and taeval *',nterms,0)

      call mpzero(nd,locexp1,nterms3)
      nd = 1
      call h3dformtacd(nd,zk,rscale1,sources,charge,
     1       dipvec,ns,c3,nterms3,locexp1,wlege,nlege)
      call errprint(pot,opot,fld,ofld)
      pot = 0
      fld(1) = 0
      fld(2) = 0
      fld(3) = 0
      nd = 1
      call h3dtaevalg(nd,zk,rscale1,c3,locexp1,nterms3,ztrg,nt,
     1      pot,fld,wlege,nlege)

      call errprint(pot,opot,fld,ofld)


      stop
      end
c
c
c
c
c
C
C
      subroutine errprint(pot,opot,fld,ofld)
      implicit real *8 (a-h,o-z)
      complex *16 pot,opot,fld(3),ofld(3)
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
      write(6,*) 
     1     ' POT error,       rel error,    FLD error,       rel error'
      write(6,1000) abs(pot-opot),
     1 		abs(pot-opot)/abs(opot),err,err/ddd
      write(13,*) 
     1     ' POT error,       rel error,    FLD error,       rel error'
      write(13,1000) abs(pot-opot),
     1 		abs(pot-opot)/abs(opot),err,err/ddd
      return
      end
c
