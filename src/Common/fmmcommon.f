c
c
c        a collection of common subroutines used in the fmm
c           
c          mpalloc - allocate workspace for multipole and local
c                     expansions
c
c          dreorderf - permute a real array with given permutation
c
c          dreroderi - pertmute a real array with inverse of
c                      given permutation
c 
c          drescale - rescale a vector with a scalar
c 
c          mpzero - zero out a multipole/local expansion
c  
c          mpadd - add a multipole expansion to an existing one
c
c          mpscale - scale multipole expansion coeffs 
c
c          cart2polar - cartesian to polar coordinates
c
c          getpwrotmat - get a collection of rotation matrices
c                  for pw stuff
c
c          getsqrtbinomialcoeffs - get square root of binomial
c                         coeffs array
c          bnlcft - binomial coeffs and square roots
c          fstrtn - generate rotation matrix
c          geterrstr - error string generator for something
c
c
c
      subroutine mpalloc(nd,laddr,iaddr,nlevels,lmptot,nterms)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i and iaddr(2,i) points to the local
c     expansion of box i
c  
c     Input arguments
c     nd          in: Integer
c                 number of expansions per box            
c 
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array provinding access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     nterms      in: Integer(0:nlevels)
c                 Number of terms requried in expansions at each
c                 level
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr       out: Integer *8(2,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer *8
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels),nd
      integer *8 iaddr(2,*), lmptot 
      integer laddr(2,0:nlevels)
      integer ibox,i,iptr
      integer *8 istart,nn,itmp

      istart = 1
      do i = 0,nlevels

        nn = (2*nterms(i)+1)*2*(nterms(i)+1)*nd 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)

         do ibox = laddr(1,i),laddr(2,i)
c            Allocate memory for the multipole expansion         

             itmp = ibox-laddr(1,i)
             iaddr(1,ibox) = istart + itmp*2*nn

c            Allocate memory for the local expansion
             iaddr(2,ibox) = istart + itmp*2*nn + nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*2*nn
      enddo
      lmptot = istart

      return
      end
c----------------------------------------------------------------     

      subroutine dreorderf(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the sorting order defined by
c        iarr
c
c        arrsort(j,i) = arr(j,iarr(i)), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer ndim,idim,i,n
      double precision arr(ndim,*),arrsort(ndim,*)
      integer iarr(*)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,i) = arr(idim,iarr(i))
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c--------------------------------------------------
      subroutine dreorderi(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the inverse of the
c        sorting order defined by
c        iarr.
c
c        Note that this subroutine is the inverse of 
c        dreorderf
c
c        arrsort(j,iarr(i)) = arr(j,i), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer i,idim,ndim,n
      double precision arr(ndim,*),arrsort(ndim,*)
      integer iarr(*)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,iarr(i)) = arr(idim,i)
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c----------------------------------------------------------      





      subroutine ireorderf(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the sorting order defined by
c        iarr
c
c        arrsort(j,i) = arr(j,iarr(i)), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer ndim,idim,i,n
      integer arr(ndim,*),arrsort(ndim,*)
      integer iarr(*)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,i) = arr(idim,iarr(i))
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c--------------------------------------------------
      subroutine ireorderi(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the inverse of the
c        sorting order defined by
c        iarr.
c
c        Note that this subroutine is the inverse of 
c        dreorderf
c
c        arrsort(j,iarr(i)) = arr(j,i), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer i,idim,ndim,n
      integer arr(ndim,*),arrsort(ndim,*)
      integer iarr(*)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,iarr(i)) = arr(idim,i)
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c
c
      subroutine drescale(n,a,r)
      implicit none
      real *8 a(n),r
      integer i,n

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,n
        a(i) = a(i)*r
      enddo
C$OMP END PARALLEL DO

      


      return
      end





      subroutine mpzero(nd,mpole,nterms)
c
cc      this subroutine zeros out a collection of 
c       multipole expansions
c
c       input:
c       nd     in:integer
c              number of expansions
c
c       nterms in:integer
c              number of terms in the expansions
c
c       inout:
c       mpole  inout: double complex(nd,0:nterms,-nterms:nterms)
c              multipole expansions to be zeroed out
c
      implicit none
      integer nd, i,j,nterms,idim
      double complex mpole(nd,0:nterms,-nterms:nterms)

      do i=-nterms,nterms
        do j=0,nterms
          do idim=1,nd
            mpole(idim,j,i) = 0
          enddo
        enddo
      enddo

      return
      end
c
c
c
c
c
      subroutine mpadd(nd,mpolein,mpoleout,nterms)
c
cc      this subroutine add mpolein to mpoleout 
c
c       input:
c       nd     in:integer
c              number of expansions
c
c       nterms in:integer
c              number of terms in the expansions
c
c       mpolein in: complex *16(nd,0;nterms,-nterms:nterms)
c               multipole expansion to be added
c
c       inout:
c       mpoleout  inout: double complex(nd,0:nterms,-nterms:nterms)
c              multipole expansions to which mpolein is to be added
c
      implicit none
      integer nd, i,j,nterms,idim
      double complex mpolein(nd,0:nterms,-nterms:nterms)
      double complex mpoleout(nd,0:nterms,-nterms:nterms)

      do i=-nterms,nterms
        do j=0,nterms
          do idim=1,nd
            mpoleout(idim,j,i) = mpoleout(idim,j,i) + mpolein(idim,j,i)
          enddo
        enddo
      enddo

      return
      end
c----------------------------------------------------------------

      subroutine mpscale(nd,nterms,mpolein,rsc,mpoleout)
c
cc      this subroutine rescales a multipole 
c       expansion where mpoleout(i,j) = mpolein(i,j)*rscalepow(i)
c
c       input
c       nd      in: integer
c               number of multipole expansions
c
c       nterms  in: integer
c               number of terms in the expansion
c    
c       rsc     in: double precision(0:nterms)
c               scaling factor for the multipole expansions            
c
c       mpolein   in: double complex(nd,0:nterms,-nterms:nterms)
c                input multipole expansions
c
c       output
c       mpoleout  out: double complex(nd,0:nterms,-nterms:nterms)
c                output multipole expansions
c
      implicit none
      integer nd,nterms,i,j,idim
      double precision rsc(0:nterms)
      double complex mpolein(nd,0:nterms,-nterms:nterms)
      double complex mpoleout(nd,0:nterms,-nterms:nterms)

      do j=-nterms,nterms
        do i=0,nterms
          do idim=1,nd
            mpoleout(idim,i,j) = mpolein(idim,i,j)*rsc(i)
          enddo
        enddo
      enddo

      return
      end
c
c
c
c
c**********************************************************************
      subroutine cart2polar(zat,r,theta,phi)
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c	zat   :  Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c	r     :  |zat|
c	theta : angle subtended with respect to z-axis
c	phi   : angle of (zat(1),zat(2)) subtended with 
c               respect to x-axis
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 zat(3)
      complex *16 ephi,eye
      data eye/(0.0d0,1.0d0)/
c
c 
      r= sqrt(zat(1)**2+zat(2)**2+zat(3)**2)
      proj = sqrt(zat(1)**2+zat(2)**2)
c
      theta = datan2(proj,zat(3))
      if( abs(zat(1)) .eq. 0 .and. abs(zat(2)) .eq. 0 ) then
      phi = 0
      else
      phi = datan2(zat(2),zat(1))
      endif
      return
      end
c
c
c
c
c
      subroutine getsqrtbinomialcoeffs(n,dc)
      implicit none
      integer i,j,n
      real *8 dc(0:n,0:n)
      real *8, allocatable :: d(:,:)
      allocate(d(0:n,0:n))

      do i=0,n
        do j=0,n
          d(i,j) = 0
          dc(i,j) = 0
        enddo
      enddo

      do i=0,n
        d(i,0) = 1
        dc(i,0) = 1
      enddo
      do i=1,n
        d(i,i) = 1
        dc(i,i) = 1
        do j=i+1,n
          d(j,i) =d(j-1,i)+d(j-1,i-1)
          dc(j,i) = sqrt(d(j,i))
        enddo
      enddo

      return
      end



c***********************************************************************
      subroutine getpwrotmat(nterms,carray,rdpi2,rdmpi2,rdsq3,rdmsq3,
     1   dc)
c***********************************************************************
      implicit double precision (a-h,o-z)
      double precision pi
      double precision carray(0:4*nterms,0:4*nterms)
      double precision dc(0:4*nterms,0:4*nterms)
      double precision rdpi2(0:nterms,0:nterms,-nterms:nterms) 
      double precision rdmpi2(0:nterms,0:nterms,-nterms:nterms) 
      double precision rdsq3(0:nterms,0:nterms,-nterms:nterms) 
      double precision rdmsq3(0:nterms,0:nterms,-nterms:nterms) 
      integer nterms
c
c----- call initialization routines
c
      call prini(6,13)
      pi = 4*datan(1.0d0)
c
      call bnlcft(carray,dc,4*nterms)
      theta = pi/2.0d0
      call fstrtn(nterms,rdpi2,dc,theta)
      theta = -pi/2.0d0
      call fstrtn(nterms,rdmpi2,dc,theta)
      theta = dacos(dsqrt(3.0d0)/3.0d0)
      call fstrtn(nterms,rdsq3,dc,theta)
      theta = dacos(-dsqrt(3.0d0)/3.0d0)
      call fstrtn(nterms,rdmsq3,dc,theta)
      return
      end
c
c***********************************************************************
      subroutine zflip(nterms,mpole,mrotate)
c***********************************************************************
      implicit double precision (a-h,o-z)
      double complex mpole(0:nterms,0:nterms)
      double complex mrotate(0:nterms,0:nterms)
c
      do 200 l=0,nterms,2
         do 100 m=0,l
            mrotate(l,m)=dconjg(mpole(l,m))
100      continue
200   continue
      do 400 l=1,nterms,2
         do 300 m=0,l
            mrotate(l,m)=-dconjg(mpole(l,m))
300      continue
400   continue
ccc      call prinf(' mrotate is *',nterms,0)
ccc      call prinm(mrotate,nterms)
      return
      end
c***********************************************************************
      subroutine prinout(mpole,ll,nterms)
c***********************************************************************
c
c     printing routine displays multipole moments one degree (l)
c     at a time , from m = 0 to m = l  .
c
c***********************************************************************
      double complex mpole(0:nterms,0:nterms)
      integer nterms
c
      do 100 l = 0,ll
	 write(6,1000)(mpole(l,m),m=0,ll)
	 write(13,1000)(mpole(l,m),m=0,ll)
ccc	 write(6,1001)
ccc	 write(13,1001)
100   continue
1000  format(6d12.5)
1001  format(/)
      return
      end
c
c***********************************************************************
      subroutine bnlcft(c, sqc, nterms)
c***********************************************************************
c
c     usage:
c           computes the binomial coefficients c_nterms^n, where
c           n=0,1,2,...,nterms.
c     input:
c           nterms: an integer indicates the number we are going
c                 to choose from.
c     output:
c           c:    an array consists of the binomial coefficients.
c           sqc:   an array consists of the square root of the
c                 binormial coefficients.
c
c***********************************************************************
c
      implicit double precision (a-h,o-z)
      integer nterms,n,m
      double precision c(0:nterms,0:nterms)
      double precision sqc(0:nterms,0:nterms)
c
      do 100 n=0,nterms
        c(n,0)=1.0d0
        sqc(n,0)=1.0d0
100   continue
      do 300 m=1,nterms
         c(m,m)=1.0d0
         sqc(m,m)=1.0d0
         do 200 n=m+1,nterms
            c(n,m)=c(n-1,m)+c(n-1,m-1)
            sqc(n,m)=dsqrt(c(n,m))
200      continue
300   continue
c
      return
      end
c
c***********************************************************************
      subroutine fstrtn(nterms,d,sqc,theta)
c***********************************************************************
c
c     usage:
c           implement the fast version of rotation matrices from
c           the recurrences formulas.
c     input:
c           nterms: an integer indicates the dimension of d.
c           sqc:    an array contains the square root of the
c                   binormial coefficients.
c           theta:  the rotate angle about the y-axis.
c     output:
c           d:      an array which contains the rotation matrix.
c                   note: only half of d are evaluated, the other
c                   half can be obtained by using the symmetricity.
c
c***********************************************************************
c
c---------------declares of fstrtn--------------------------------------
c
      implicit double precision (a-h,o-z)
      integer nterms
c
ccc     parameter (nterms=8)
ccc     double precision d(0:nterms,-nterms:nterms,-nterms:nterms)
ccc     double precision c(0:4*nterms, 0:4*nterms)
c
      double precision d(0:nterms,0:nterms,-nterms:nterms)
      double precision sqc(0:4*nterms, 0:4*nterms)
      double precision ctheta, stheta, hsthta
      double precision cthtap, cthtan, theta, precis
c
      data precis/1.0d-19/
      data ww/0.7071067811865476d+00/
c
c---------------executions of fstrtn------------------------------------
c
      ctheta=dcos(theta)
      if (dabs(ctheta).le.precis) ctheta=0.0d0
      stheta=dsin(-theta)
      if (dabs(stheta).le.precis) stheta=0.0d0
c     stheta=dsqrt(1.0d0-ctheta*ctheta)
      hsthta=ww*stheta
      cthtap=ww*(1.0d0+ctheta)
      cthtan=-ww*(1.0d0-ctheta)
c
c     initial setup for some coefficient matrix.
c
c     call bnlcft(c, sqc, 4*nterms)
c
      d(0,0,0)=1.0d0
c
      do 1000 ij=1,nterms
c        if (ij.ge.3) stop
c
c     compute the result for m=0 case, use formula (1).
c
         do 200 im=-ij,-1
            d(ij,0,im)=-sqc(ij-im,2)*d(ij-1,0,im+1)
            if (im.gt.(1-ij)) then
               d(ij,0,im)=d(ij,0,im)+sqc(ij+im,2)*d(ij-1,0,im-1)
            endif
            d(ij,0,im)=d(ij,0,im)*hsthta
            if (im.gt.-ij) then
               d(ij,0,im)=d(ij,0,im)+
     1           d(ij-1,0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1)
            endif
            d(ij,0,im)=d(ij,0,im)/ij
200      continue
c
         d(ij,0,0)=d(ij-1,0,0)*ctheta
         if (ij.gt.1) then
            d(ij,0,0)=d(ij,0,0)+hsthta*sqc(ij,2)*(d(ij-1,0,-1)+
     1                 d(ij-1,0,1))/ij
         endif
c
         do 400 im=1,ij
            d(ij,0,im)=-sqc(ij+im,2)*d(ij-1,0,im-1)
            if (im.lt.(ij-1)) then
               d(ij,0,im)=d(ij,0,im)+sqc(ij-im,2)*d(ij-1,0,im+1)
            endif 
            d(ij,0,im)=d(ij,0,im)*hsthta
            if (im.lt.ij) then
               d(ij,0,im)=d(ij,0,im)+
     1           d(ij-1,0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1)
            endif
            d(ij,0,im)=d(ij,0,im)/ij
400      continue
c
c     compute the result for 0<m<=j case, use formula (2).
c
         do 800 imp=1,ij
            do 500 im=-ij,-1
               d(ij,imp,im)=d(ij-1,imp-1,im+1)*cthtan*sqc(ij-im,2)
               if (im.gt.(1-ij)) then
                  d(ij,imp,im)=d(ij,imp,im)-
     1             d(ij-1,imp-1,im-1)*cthtap*sqc(ij+im,2)
               endif
               if (im.gt.-ij) then
                  d(ij,imp,im)=d(ij,imp,im)+
     1             d(ij-1,imp-1,im)*stheta*sqc(ij+im,1)*sqc(ij-im,1)
               endif 
               d(ij,imp,im)=d(ij,imp,im)*ww/sqc(ij+imp,2)
500         continue
c
            d(ij,imp,0)=ij*stheta*d(ij-1,imp-1,0)
            if (ij.gt.1) then
               d(ij,imp,0)=d(ij,imp,0)-sqc(ij,2)*(
     1            d(ij-1,imp-1,-1)*cthtap+d(ij-1,imp-1,1)*cthtan)
            endif
            d(ij,imp,0)=d(ij,imp,0)*ww/sqc(ij+imp,2)
c
            do 600 im=1,ij
               d(ij,imp,im)=d(ij-1,imp-1,im-1)*cthtap*sqc(ij+im,2)
               if (im.lt.(ij-1)) then
                  d(ij,imp,im)=d(ij,imp,im)-
     1             d(ij-1,imp-1,im+1)*cthtan*sqc(ij-im,2)
               endif
               if (im.lt.ij) then
                  d(ij,imp,im)=d(ij,imp,im)+
     1             d(ij-1,imp-1,im)*stheta*sqc(ij+im,1)*sqc(ij-im,1)
               endif
               d(ij,imp,im)=d(ij,imp,im)*ww/sqc(ij+imp,2)
600         continue
c
c     use symmetricity, i.e. formula (3.80) in biedenharn & loucks
c     book, to compute the lower part of the matrix
c
c           do im=-ij,ij
c              d(ij,-imp,im)=d(ij,imp,-im)
c           enddo
c
800      continue
1000  continue
c
c
      return
      end
c
c-------------------------------------------------------
      subroutine geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      implicit real *8 (a-h,o-z)
      character(len=*) str1
      character(len=13) str2
      character(len=14) str3
      character(len=19) str4
      character(len=30) str5

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

      if(ipgh.eq.3.or.ifpghtarg.eq.3) then
        str5=" pot, grad, and hess test"
        len1 = len1+25
      endif

      str1 = str2//trim(str3)//trim(str4)//trim(str5)

      return
      end
