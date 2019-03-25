c***********************************************************************
      subroutine rotgen(nterms,carray,rdpi2,rdmpi2,rdsq3,rdmsq3,dc)
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
c----------end of fstrtn------------------------------------------------
c
      return
      end
c
