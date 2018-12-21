c
cc      plane wave routines for Helmholtz 3D FMM
c 
c***********************************************************************
      subroutine rlscini(rlsc,nlambs,rlams,zk,nterms)
c***********************************************************************
c
c       this subroutine computes p_{n,m}(i\sqrt(\lambda^2 - k^2)/k) (-i)^n
c       for all combinations of \lambda,n,m required for plane wave
c       routines in a Helmholtz fmm
c
c       input:
c       nlambs -     number of discretization nodes in lambda quadrature
c       rlams -      location of discretization nodes in lambda
c                    quadrature. Note here rlams stores sqrt(\lambda^2 -
c                    k^2)
c       zk -         Helmholtz parameter(k)
c       nterms -     max n for which p_{n,m}'s need to be computed
c
c       output:
c       rlsc(0:nterms,0:nterms,nlambs) -  rlsc(n,m,nl) = p_{n,m}(
c                                             i rlams(nl)/zk) (-i)^n
c
c
      implicit real *8 (a-h,o-z)
      complex *16 rlsc(0:nterms,0:nterms,nlambs),rlams(nlambs)
      complex *16 zmult,ipow(0:nterms),ima,zk

      done = 1
      ima = cmplx(0,done)
      zmult = cmplx(0,-done)
      ipow(0) = done

      do 900 i=1,nterms

      ipow(i) = ipow(i-1)*zmult

900   continue      
c
      do 2000 nl = 1,nlambs

      zmult = ima*rlams(nl)/zk
      call zylgndr(nterms,zmult,rlsc(0,0,nl))

      do 1900 iin=1,nterms
      do 1800 im=0,nterms

      rlsc(iin,im,nl) = rlsc(iin,im,nl)*ipow(iin)

 1800 continue
 1900 continue
      

2000  continue
      return
      end
c***********************************************************************
      subroutine mkexps(rlams,nlambs,numphys,nexptotp,zk,xs,ys,zs)
      implicit real *8 (a-h,o-z)
      complex *16 ima,zk,rk
      complex *16 xs(-5:5,nexptotp)
      complex *16 ys(-5:5,nexptotp)
      complex *16 zs(5,nexptotp),rlams(nlambs)
      real *8   u
      integer  nlambs,numphys(nlambs),nexptotp
      data ima/(0.0d0,1.0d0)/
c
c     this subroutine computes the tables of exponentials needed
c     for translating exponential representations of harmonic
c     functions, discretized via normans quadratures.
c
c     u   = \int_0^\infty e^{-\lambda z}
c           \int_0^{2\pi} e^{i\lambda(x cos(u)+y sin(u))}
c           mexpphys(lambda,u) du dlambda
c
c     mexpphys(*):  discrete values of the moment function 
c                   m(\lambda,u), ordered as follows.
c
c         mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),..., 
c              m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
c         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
c              m(\lambda_2,0),...,
c                  m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
c         etc.
c
c     on input:
c
c     rlams(nlambs)   discretization points in lambda integral.
c                     Note rlams stores \sqrt(\lambda^2 - k^2)
c     nlambs          number of discret. pts. in lambda integral
c     numphys(j)     number of nodes in u integral needed 
c                    for corresponding lambda =  lambda_j. 
c     nexptotp        sum_j numphys(j)
c     zk              Helmholtz parameter (k)
c
c     on output:
c
c     xs(-1,nexptotp)   e^{-i \lambda_j (cos(u_k)}  in above ordering
c     xs(-2,nexptotp)   e^{-i \lambda_j (2 cos(u_k)}  in above ordering.
c     xs(-3,nexptotp)   e^{-i \lambda_j (3 cos(u_k)}  in above ordering.
c     xs(-4,nexptotp)   e^{-i \lambda_j (4 cos(u_k)}  in above ordering.
c     xs(-5,nexptotp)   e^{-i \lambda_j (5 cos(u_k)}  in above ordering.
c     xs(1,nexptotp)   e^{i \lambda_j (cos(u_k)}  in above ordering
c     xs(2,nexptotp)   e^{i \lambda_j (2 cos(u_k)}  in above ordering.
c     xs(3,nexptotp)   e^{i \lambda_j (3 cos(u_k)}  in above ordering.
c     xs(4,nexptotp)   e^{i \lambda_j (4 cos(u_k)}  in above ordering.
c     xs(5,nexptotp)   e^{i \lambda_j (5 cos(u_k)}  in above ordering.
c     ys(-1,nexptotp)   e^{i \lambda_j (sin(u_k)}  in above ordering.
c     ys(-2,nexptotp)   e^{i \lambda_j (2 sin(u_k)}  in above ordering.
c     ys(-3,nexptotp)   e^{i \lambda_j (3 sin(u_k)}  in above ordering.
c     ys(-4,nexptotp)   e^{i \lambda_j (4 sin(u_k)}  in above ordering.
c     ys(-5,nexptotp)   e^{i \lambda_j (5 sin(u_k)}  in above ordering.
c     ys(1,nexptotp)   e^{i \lambda_j (sin(u_k)}  in above ordering.
c     ys(2,nexptotp)   e^{i \lambda_j (2 sin(u_k)}  in above ordering.
c     ys(3,nexptotp)   e^{i \lambda_j (3 sin(u_k)}  in above ordering.
c     ys(4,nexptotp)   e^{i \lambda_j (4 sin(u_k)}  in above ordering.
c     ys(5,nexptotp)   e^{i \lambda_j (5 sin(u_k)}  in above ordering.
c     zs(1,nexptotp)   e^{-\sqrt(\lambda_j^2-k^2)}     in above ordering.
c     zs(2,nexptotp)    e^{-2 \sqrt(\lambda_j^2-k^2)}   in above ordering. 
c     zs(3,nexptotp)    e^{-3 \sqrt(\lambda_j^2-k^2)}   in above ordering. 
c     zs(4,nexptotp)    e^{-4 \sqrt(\lambda_j^2-k^2)}   in above ordering. 
c     zs(5,nexptotp)    e^{-5 \sqrt(\lambda_j^2-k^2)}   in above ordering. 
c------------------------------------------------------------
c      
c     loop over each lambda value 
c
      pi = 4*datan(1.0d0)
      ntot = 0
      do 400 nl = 1,nlambs
         hu=2*pi/numphys(nl)

         rk = sqrt(rlams(nl)**2 + zk**2)

         do 200 mth = 1,numphys(nl)
            u = (mth-1)*hu
            ncurrent = ntot+mth
            zs(1,ncurrent) = cdexp(-rlams(nl) )
            zs(2,ncurrent) = cdexp(-2.0d0*rlams(nl) )
            zs(3,ncurrent) = cdexp(-3.0d0*rlams(nl) )
            zs(4,ncurrent) = cdexp(-4.0d0*rlams(nl) )
            zs(5,ncurrent) = cdexp(-5.0d0*rlams(nl) )
            xs(-1,ncurrent) = cdexp(-ima*rk*cos(u))
            xs(-2,ncurrent) = cdexp(-ima*rk*2.0d0*cos(u))
            xs(-3,ncurrent) = cdexp(-ima*rk*3.0d0*cos(u))
            xs(-4,ncurrent) = cdexp(-ima*rk*4.0d0*cos(u))
            xs(-5,ncurrent) = cdexp(-ima*rk*5.0d0*cos(u))
            xs(0,ncurrent) = 1
            xs(1,ncurrent) = cdexp(ima*rk*cos(u))
            xs(2,ncurrent) = cdexp(ima*rk*2.0d0*cos(u))
            xs(3,ncurrent) = cdexp(ima*rk*3.0d0*cos(u))
            xs(4,ncurrent) = cdexp(ima*rk*4.0d0*cos(u))
            xs(5,ncurrent) = cdexp(ima*rk*5.0d0*cos(u))
            ys(-1,ncurrent) = cdexp(-ima*rk*dsin(u))
            ys(-2,ncurrent) = cdexp(-ima*rk*2.0d0*dsin(u))
            ys(-3,ncurrent) = cdexp(-ima*rk*3.0d0*dsin(u))
            ys(-4,ncurrent) = cdexp(-ima*rk*4.0d0*dsin(u))
            ys(-5,ncurrent) = cdexp(-ima*rk*5.0d0*dsin(u))
            ys(0,ncurrent) = 1
            ys(1,ncurrent) = cdexp(ima*rk*dsin(u))
            ys(2,ncurrent) = cdexp(ima*rk*2.0d0*dsin(u))
            ys(3,ncurrent) = cdexp(ima*rk*3.0d0*dsin(u))
            ys(4,ncurrent) = cdexp(ima*rk*4.0d0*dsin(u))
            ys(5,ncurrent) = cdexp(ima*rk*5.0d0*dsin(u))
200      continue
         ntot = ntot+numphys(nl)
400   continue
      return
      end
c***********************************************************************
      subroutine mkfexp(nlambs,numfour,numphys,fexp,fexp2)
      implicit real *8 (a-h,o-z)
      complex *16 ima
      complex *16 fexp(1),fexp2(1)
      integer  nlambs,numphys(nlambs),numfour(nlambs)
      data ima/(0.0d0,1.0d0)/
c
c     this subroutine computes the tables of exponentials needed
c     for mapping from fourier to physical domain. 
c     in order to minimize storage, they are organized in a 
c     one-dimenional array corresponding to the order in which they
c     are accessed by subroutine ftophys.
c    
c***********************************************************************
      pi = 4*datan(1.0d0)

      next = 1
      do 1600 i=1,nlambs

	  nalpha = numphys(i)
      halpha=2*pi/nalpha

      do 1400 j=1,nalpha
	  do 1200 mm = 1,numfour(i)

      alpha=(j-1)*halpha
      fexp(next) = cdexp(ima*mm*alpha)
	  next = next + 1

1200  continue
1400  continue
1600  continue


      next = 1
      do 2600 i=1,nlambs

	  nalpha = numphys(i)
      halpha=2*pi/nalpha

	  do 2400 mm = 1,numfour(i)
      do 2200 j=1,nalpha

      alpha=(j-1)*halpha
      fexp2(next) = cdexp(-ima*mm*alpha)
	  next = next + 1

2200  continue
2400  continue
2600  continue

      return
      end
c***********************************************************************

      subroutine mpoletoexp(mpole,nterms,nlambs,numtets,nexptot,
     1                mexpupf,mexpdownf,rlsc)

c     This subroutine converts a multipole expansion into the
c     corresponding exponential moment function mexp for
c     both the +z direction and the -z direction
c
c     U(x,y,z) = \sum_{n=0}^{nterms} \sum_{m=-n,n} mpole(n,m)
c                P_n^m (\cos(\theta)) e^{i m \phi} h_{n} (kr)
c  
c              = (1/2\pi) \int_{0}^{\infty} e^{-\sqrt(\lambda^2-k^2) z}
c                \int_{0}^{2\pi} e^{i \lambda (x \cos(\alpha) +
c                y \sin(\alpha))} mexpup(\lambda,\alpha)
c                d\alpha d\lambda
c 
c     for the +z direction and
c
c              = (1/2\pi) \int_{0}^{\infty} e^{\sqrt(\lambda^2-k^2) z}
c                \int_{0}^{2\pi} e^{-i \lambda (x \cos(\alpha) +
c                y \sin(\alpha))} mexpdown(\lambda,\alpha)
c                d\alpha d\lambda
c 
c     for the -z direction
c
c     NOTE: The expression for -z corresponds to the mapping
c     (x,y,z) -> (-x,-y,-z), ie reflection through
c     the origin.
c
c     NOTE 2: The multipole expansion is assumed to have been
c     rescaled so that the box containing the sources has unit
c     dimension
c
c     NOTE 3: Since we store the exponential moment function in
c     Fourier domain (w.r.t the \alpha variable) we compute
c 
c     M_\lambda(m) = (i)**m (-i)**n \sum_{n=m}^{N} P_{n,m}(i
c                       \sqrt(\lambda^2 -k^2)/k) mpole(n,m)
c      
c
c     INPUT arguments
c     mpole       in: complex *16 (0:nterms, -nterms:nterms)
c                 The multipole expansion 
c  
c     nterms:     in: integer
c                 Order of the multipole expansion
c
c     nlambs      in: integer
c                 number of discretization points in the \lambda
c                 integral
c
c     numtets     in: integer(nlambs)
c                 number of fourier modes needed in expansion
c                 of \alpha variable for each \lambda variable
c
c     nexptot     in: integer
c                 nexptot = \sum_{j} numtets(j)
c
c     rlsc        in: real *8(nlambs, 0:nterms, 0:nterms)
c                 scaled discretization points in the \lambda
c                 integral
c
c     OUTPUT 
c     mexpupf     out: complex *16 (nexptot)
c                 Fourier coefficients of the function
c                 mexpup(\lambda,\alpha) for successive
c                 discrete lambda values. They are ordered as
c                 follows
c
c                 mexpupf(1,...., 2*numtets(1)+1) = fourier modes
c                             for \lambda_1
c                 a_{0},a_{1},a_{-1},a_{2},a_{-2},...,a_{numtets(1)},
c                 a_{-numtets(1)}
c
c                 mexpupf(2*numtets(1)+2,...., numtets(2) = fourier
c                 modes for \lambda_2
c
c                 ETC
c
c     mexpdownf   out: complex *16 (nexptot)
c                 Fourier coefficients of the function 
c                 mexpdown(\lambda,\alpha) for successive
c                 discrete \lambda values
c---------------------------------------------------------------

      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mexpupf(nexptot)
      complex *16 mexpdownf(nexptot)
      complex *16 zeyep,ztmp1,ztmp2,ztmp3,ztmp4
      complex *16 rlsc(0:nterms,0:nterms,nlambs)

c     Temp variables
      real *8 sgn
      integer ntot,ncurrent,nl,mth,nm

      ntot = 0
      do nl=1,nlambs
         sgn = 1.0d0
 
         ztmp1 = 0
         ztmp2 = 0

         mth = 0 
         ncurrent = ntot + 1

         do nm = mth,nterms,2
            ztmp1 = ztmp1 + rlsc(nm,mth,nl)*mpole(nm,mth)
         enddo

         do nm=mth+1,nterms,2
            ztmp2 = ztmp2 + rlsc(nm,mth,nl)*mpole(nm,mth)
         enddo

         mexpupf(ncurrent) = ztmp1+ztmp2
         mexpdownf(ncurrent) = sgn*(ztmp1-ztmp2)

         do mth = 1,numtets(nl)
            ncurrent = ntot + 2*mth
            ztmp1 = 0
            ztmp2 = 0
            ztmp3 = 0
            ztmp4 = 0

            sgn = -sgn
            do nm = mth,nterms,2
               ztmp1 = ztmp1 + rlsc(nm,mth,nl)*mpole(nm,mth)
            enddo

            do nm=mth+1,nterms,2
               ztmp2 = ztmp2 + rlsc(nm,mth,nl)*mpole(nm,mth)
            enddo

            do nm = mth,nterms,2
               ztmp3 = ztmp3 + rlsc(nm,mth,nl)*mpole(nm,-mth)
            enddo

            do nm=mth+1,nterms,2
               ztmp4 = ztmp4 + rlsc(nm,mth,nl)*mpole(nm,-mth)
            enddo

            mexpupf(ncurrent) = ztmp1+ztmp2
            mexpdownf(ncurrent) = sgn*(ztmp1-ztmp2)

            mexpupf(ncurrent+1) = ztmp3+ztmp4
            mexpdownf(ncurrent+1) = sgn*(ztmp3-ztmp4)
         enddo
         ntot = ntot + 2*numtets(nl)+1
      enddo

      return
      end

c-----------------------------------------------------------------      
      subroutine exptolocal(local,nterms,zk,rlambs,whts,nlambs,
     1       numtets,nthmax,nexptot,lexp1f,lexp2f,scale,rlsc)
c-----------------------------------------------------------------
c     This sburoutine converts the Fourier representation of
c     an exponential moment function into a local
c     multipole expansion (with respect to the same box center)
c
c     (+z direction and -z combined)
c
c     u(x,y,z) = \int_{0}^{\infty} e^{-\sqrt(\lambda-k^2) z}
c                \int_{0}^{2\pi} e^{i\lambda(x\cos(\alpha) +
c                y\sin(\alpha))} lexp1 (\lambda,\alpha) 
c                P_{n,m}(i \sqrt(\lambda^2 -k^2)/k) \lambda/
c                \sqrt(\lambda^2 - k^2)
c                d\lambda
c                d\alpha
c
c              = \sum_{0}^{nterms} \sum_{m=-n,n} local(n,m) Y_n^m
c                 (\cos(\theta)) e^{i m \phi} r^{n}
c
c     INPUT arguments
c     nterms           in: integer
c                      Order of local expansion
c
c     zk               in: complex *16
c                      Helmholtz parameter
c
c     rlambs           in: complex *16(nlambs)
c                      discretization points in the \lambda integral
c
c     whts             in: complex *16(nlambs)
c                      quadrature weights in \lambda integral
c
c     nlambs           in: integer
c                      number of discretization points in \lambda
c                      integral
c
c     numtets          in: integer(nlambs)
c                      number of fourier modes in expansion of
c                      \alpha variable for \lambda_j
c
c     nthmax           in: integer
c                      max_j numtets(j)
c
c     nexptot          in: integer
c                      sum_j numtets(j)
c                      
c
c     lexp1f(nexptot)  complex *16(nexptot)
c                      Fourier coefficients of the function 
c                      lexp1 for discrete \lambda values in the +z
c                      direction.
c
c                      They are ordered as follows:
c
c                      lexp1f(1,...,2*numtets(1)+1) = Fourier modes
c                      for \lambda_1
c                      a_{0},a_{1},a_{-1},a_{2},a_{-2},...a_{numtets(1)},
c                      a_{-numtets(1)}
c                      lexp1f(2*numtets(1)+2,...,numtets(2) = Fourier
c                      modes for \lambda_2 etc
c
c     lexp2f(nexptot)  complex *16(nexptot)
c                      Fourier coefficients of the function 
c                      lexp2 for discrete \lambda values in the -z
c                      direction.
c
c                      They are ordered as follows:
c
c                      lexp2f(1,...,2*numtets(1)+1) = Fourier modes
c                      for \lambda_1
c                      a_{0},a_{1},a_{-1},a_{2},a_{-2},...a_{numtets(1)},
c                      a_{-numtets(1)}
c                      lexp2f(2*numtets(1)+2,...,numtets(2) = Fourier
c                      modes for \lambda_2 etc
c
c     scale            in: real *8
c                      scaling parameter for local expansion
c
c     OUTPUT
c     local(0:nterms,-nterms:nterms): output local expansion of order
c                                     nterms
        
      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot,nthmax
      integer ncurrent,ntot,nl,ncurrent2
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 ima,zmult
      complex *16 lexp1f(nexptot),lexp2f(nexptot)
      complex *16 rlambs(nlambs), whts(nlambs)
      complex *16 rlsc(0:nterms,0:nterms,nlambs)
      complex *16 zk,zk2
      real *8 scale, rscale(0:nterms)
    
c     Temporary variables
      integer i, nm, mth, j, mmax
      real *8 done

      done = 1
      ima = cmplx(0,done)
      zk2 = 1/(ima*zk)

      rscale(0) = 1

      do nm=0,nterms

         if(nm.gt.0) then

         rscale(nm) = rscale(nm-1)*scale

         endif

         do mth = -nterms,nterms
            local(nm,mth) = 0.0d0
         enddo
      enddo

      ntot = 0
      do nl=1,nlambs
c        Add contributions to local expansion
         do nm=0,nterms,2
            mmax = numtets(nl)
            if(mmax.gt.nm) mmax = nm

            ncurrent = ntot+1
            local(nm,0) = local(nm,0)+(lexp2f(ncurrent)+
     1          lexp1f(ncurrent))*rlsc(nm,0,nl)*whts(nl)
            do mth=1,mmax
               ncurrent = ntot+2*mth
               ncurrent2 = ntot+2*mth+1
               local(nm,mth) = local(nm,mth)+(lexp2f(ncurrent)+
     1           lexp1f(ncurrent))*rlsc(nm,mth,nl)*whts(nl)
               local(nm,-mth) = local(nm,-mth)+(lexp2f(ncurrent2)+
     1           lexp1f(ncurrent2))*rlsc(nm,mth,nl)*whts(nl)
            enddo
         enddo
         do nm=1,nterms,2
            mmax = numtets(nl) 
            if(mmax.gt.nm) mmax = nm
            ncurrent = ntot+1
            local(nm,0) = local(nm,0)+(lexp2f(ncurrent)-
     1         lexp1f(ncurrent))*rlsc(nm,0,nl)*whts(nl)
            do mth =1,mmax
               ncurrent = ntot+2*mth
               ncurrent2 = ntot+2*mth+1
               local(nm,mth) = local(nm,mth)+(lexp2f(ncurrent)-
     1            lexp1f(ncurrent))*rlsc(nm,mth,nl)*whts(nl)
               local(nm,-mth) = local(nm,-mth)+(lexp2f(ncurrent2)-
     1             lexp1f(ncurrent2))*rlsc(nm,mth,nl)*whts(nl)
            enddo
         enddo
         ntot = ntot + 2*numtets(nl)+1
      enddo

      do nm=0,nterms
         do mth = 0,nm
            local(nm,mth) = local(nm,mth)*rscale(nm)*zk2
            if(mth.gt.0) local(nm,-mth) = local(nm,-mth)
     1          *rscale(nm)*zk2

         enddo
      enddo

      return
      end
c***********************************************************************
      subroutine phystof(mexpf,nlambs,numfour,numphys,mexpphys,fexp2)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      complex *16 mexpf(*)
      complex *16 mexpphys(*),ima
      complex *16 fexp2(*)
      real *8     alphas(0:2000)
      integer  nlambs,numfour(nlambs),numphys(nlambs)
      data ima/(0.0d0,1.0d0)/
c***********************************************************************
c
c     this subroutine converts the discretized exponential moment function
c     into its fourier expansion.
c
c     on input:
c
c     mexpphys(*):  discrete values of the moment function 
c                   m(\lambda,\alpha), ordered as follows.
c
c         mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),..., 
c              m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
c         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
c              m(\lambda_2,0),...,
c                  m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
c         etc.
c
c     nlambs:        number of discretization pts. in lambda integral
c     numfour(j):   number of fourier modes in the expansion
c                      of the function m(\lambda_j,\alpha)
c     on output:
c
c     mexpf(*):     fourier coefficients of the function 
c                   mexp(lambda,alpha) for discrete lambda values. 
c                   they are ordered as follows:
c
c               mexpf(1,...,2*numfour(1)+1) = fourier modes for lambda_1
c               a_{0},a_{1},a_{-1},a_{2},a_{-2},...a_{numfour(1)},
c                   a_{-numfour(1)}
c               mexpf(2*numfour(1)+2,...,numfour(2)) = fourier modes
c                                              for lambda_2
c               etc.
c
c------------------------------------------------------------
      done=1.0d0
c
c
      pi=datan(done)*4
      nftot = 0
      nptot  = 0
      next  = 1

      do 2000 i=1,nlambs
        nalpha = numphys(i)
        halpha=2*pi/nalpha

        do 200 j=1,nalpha
          alphas(j)=(j-1)*halpha
200     continue

        mexpf(nftot+1) = 0.0d0
        do 400 ival=1,nalpha
           mexpf(nftot+1) = mexpf(nftot+1) + mexpphys(nptot+ival) 
400     continue
        mexpf(nftot+1) = mexpf(nftot+1)/nalpha
        do 800 mm = 1,numfour(i)
           mexpf(nftot+2*mm) = 0.0d0
           mexpf(nftot+2*mm+1) = 0.0d0

           do 600 ival=1,nalpha
              mexpf(nftot+2*mm) = mexpf(nftot+2*mm) +
     1          fexp2(next)*mexpphys(nptot+ival)
              mexpf(nftot+2*mm+1) = mexpf(nftot+2*mm+1) +
     1          conjg(fexp2(next))*mexpphys(nptot+ival)
                next = next+1
600        continue
           mexpf(nftot+2*mm) = mexpf(nftot+2*mm)/nalpha
           mexpf(nftot+2*mm+1) = mexpf(nftot+2*mm+1)/nalpha
800     continue
        nftot = nftot+2*numfour(i)+1
        nptot = nptot+numphys(i)
 2000 continue

      return
      end
c
csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c
      subroutine ftophys(mexpf,nlambs,numfour,numphys,mexpphys,fexp)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      complex *16 mexpf(*)
      complex *16 mexpphys(*),ima,ctmp
      complex *16 fexp(*)
      real *8  alphas(0:1000)
      integer  nlambs,numfour(nlambs),numphys(nlambs)
      data ima/(0.0d0,1.0d0)/
c***********************************************************************
c
c     this subroutine evaluates the fourier expansion of the
c     exponential moment function m(\lambda,\alpha) at equispaced
c     nodes.
c
c     on input:
c
c     mexpf(*):     fourier coefficients of the function 
c                   mexp(lambda,alpha) for discrete lambda values. 
c                   they are ordered as follows:
c
c               mexpf(1,...,2*numfour(1)+1) = fourier modes for lambda_1
c               a_{0},a_{1},a_{-1},a_{2},a_{-2},...,a_{numfour(1)},
c               a_{-numfour(1)}
c               mexpf(2*numfour(1)+2,...,numfour(2)) = fourier modes
c                                              for lambda_2
c               etc.
c
c     nlambs:        number of discretization pts. in lambda integral
c     numfour(j):   number of fourier modes in the expansion
c                      of the function m(\lambda_j,\alpha)
c     fexp =      precomputed array of exponentials needed for
c                  fourier series evaluation
c
c     on output:
c
c     mexpphys(*):  discrete values of the moment function 
c                   m(\lambda,\alpha), ordered as follows.
c
c         mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),..., 
c              m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
c         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
c              m(\lambda_2,0),...,
c                  m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
c         etc.
c
c------------------------------------------------------------
      done=1.0d0
c
c
      pi=datan(done)*4
      nftot = 0
      nptot  = 0
      next = 1
      do 2000 i=1,nlambs
        do 1200 ival=1,numphys(i)
           mexpphys(nptot+ival) = mexpf(nftot+1)
           do 200 mm = 1,numfour(i)

              ctmp = mexpf(nftot+2*mm)*fexp(next)
              ctmp = ctmp + mexpf(nftot+2*mm+1)*conjg(fexp(next))
              next = next + 1
              mexpphys(nptot+ival) = mexpphys(nptot+ival) + ctmp

200    continue
1200   continue

       nftot = nftot+2*numfour(i)+1
       nptot = nptot+numphys(i)

2000  continue

      return
      end
c----------------------------------------------------------      

      subroutine processudexp(zk2,ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nuall,uall,
     3           nu1234,u1234,ndall,dall,nd5678,d5678,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpuall,mexpu5678,mexpdall,
     5           mexpd1234,xs,ys,zs,fexpback,rlsc)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer nuall,ndall,nu1234,nd5678
       integer uall(*),dall(*),u1234(*),d5678(*)
       real *8 rscale
       complex *16 zk2
       complex *16 rlams(*),whts(*)
       complex *16 tloc(0:nterms,-nterms:nterms)
       complex *16 mexp(nexptotp,nboxes,6)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nexptot),mexpdown(nexptot)
       complex *16 mexpupphys(nexptotp),mexpdownphys(nexptotp)
       complex *16 mexpuall(nexptotp),mexpdall(nexptotp)
       complex *16 mexpd1234(nexptotp),mexpu5678(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp),zs(5,nexptotp)
       complex *16 rlsc(0:nterms,0:nterms,nlams)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j
       complex *16 ztmp,zmul
     
       real *8 ctmp(3)


       do i=1,nexptotp

       mexpuall(i) = 0
       mexpdall(i) = 0
       mexpu5678(i) = 0
       mexpd1234(i) = 0

       enddo
      
   
       ctmp(1) = centers(1,ibox) - rscale/2.0d0
       ctmp(2) = centers(2,ibox) - rscale/2.0d0
       ctmp(3) = centers(3,ibox) - rscale/2.0d0
       
       do i=1,nuall
          jbox = uall(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
             mexpdall(j) = mexpdall(j) + mexp(j,jbox,2)*zmul
          enddo
       enddo

       do i=1,nu1234
          jbox = u1234(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

          do j=1,nexptotp
             zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
             mexpd1234(j) = mexpd1234(j) + mexp(j,jbox,2)*zmul
          enddo
       enddo

       do i=1,ndall
          jbox = dall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         
          do j=1,nexptotp
             zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
             mexpuall(j) = mexpuall(j) + mexp(j,jbox,1)*zmul
          enddo
       enddo

       do i=1,nd5678
          jbox = d5678(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
             mexpu5678(j) = mexpu5678(j) + mexp(j,jbox,1)*zmul
          enddo
       enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
       jbox = ichild(1,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpuall(i)
          mexpdownphys(i) = mexpdall(i) + mexpd1234(i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call h3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpuall(i)*xs(1,i)
          mexpdownphys(i) = (mexpdall(i) + mexpd1234(i))*xs(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call h3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpuall(i)*ys(1,i)
          mexpdownphys(i) = (mexpdall(i) + mexpd1234(i))*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call h3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpuall(i)*ys(1,i)*xs(1,i)
          mexpdownphys(i) = (mexpdall(i) + mexpd1234(i))*ys(-1,i)*
     1                        xs(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call h3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpuall(i)+mexpu5678(i))*zs(1,i)
          mexpdownphys(i) = mexpdall(i)/zs(1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call h3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpuall(i)+mexpu5678(i))*zs(1,i)*xs(1,i)
          mexpdownphys(i) = mexpdall(i)/zs(1,i)*xs(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call h3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpuall(i)+mexpu5678(i))*zs(1,i)*ys(1,i)
          mexpdownphys(i) = mexpdall(i)/zs(1,i)*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call h3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpuall(i)+mexpu5678(i))*zs(1,i)*ys(1,i)*
     1                       xs(1,i) 
          mexpdownphys(i) = mexpdall(i)/zs(1,i)*ys(-1,i)*xs(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call h3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processnsexp(zk2,ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nnall,nall,
     3           nn1256,n1256,nn12,n12,nn56,n56,
     4           nsall,sall,ns3478,s3478,ns34,s34,ns78,s78,mexpup,
     5           mexpdown,mexpupphys,mexpdownphys,
     6           mexpnall,mexpn3478,mexpn34,mexpn78,mexpsall,
     7           mexps1256,mexps12,mexps56,rdplus,
     8           xs,ys,zs,fexpback,rlsc)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer nnall,nsall,nn1256,ns3478,nn12,nn56,ns34,ns78
       integer nall(*),sall(*),n1256(*),s3478(*)
       integer n12(*),n56(*),s34(*),s78(*)
       real *8 rscale
       complex *16 zk2
       complex *16 rlams(*),whts(*)
       complex *16 tloc(0:nterms,-nterms:nterms)
       complex *16 tloc2(0:nterms,-nterms:nterms)
       complex *16 mexp(nexptotp,nboxes,6)
       real *8 rdplus(0:nterms,0:nterms,-nterms:nterms)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nexptot),mexpdown(nexptot)
       complex *16 mexpupphys(nexptotp),mexpdownphys(nexptotp)
       complex *16 mexpnall(nexptotp),mexpsall(nexptotp)
       complex *16 mexps1256(nexptotp),mexpn3478(nexptotp)
       complex *16 mexps12(nexptotp),mexps56(nexptotp)
       complex *16 mexpn34(nexptotp),mexpn78(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp),zs(5,nexptotp)
       complex *16 rlsc(0:nterms,0:nterms,nlams)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j
       complex *16 ztmp,zmul
     
       real *8 ctmp(3)


       do i=1,nexptotp

       mexpnall(i) = 0
       mexpsall(i) = 0
       mexpn3478(i) = 0
       mexpn34(i) = 0
       mexpn78(i) = 0
       mexps1256(i) = 0
       mexps12(i) = 0
       mexps56(i) = 0

       enddo
      
   
       ctmp(1) = centers(1,ibox) - rscale/2.0d0
       ctmp(2) = centers(2,ibox) - rscale/2.0d0
       ctmp(3) = centers(3,ibox) - rscale/2.0d0
       
       do i=1,nnall
          jbox = nall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             mexpsall(j) = mexpsall(j) + mexp(j,jbox,4)*zmul
          enddo

       enddo

       do i=1,nn1256
          jbox = n1256(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             mexps1256(j) = mexps1256(j) + mexp(j,jbox,4)*zmul
          enddo
       enddo

       do i=1,nn12
          jbox = n12(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             mexps12(j) = mexps12(j) + mexp(j,jbox,4)*zmul
          enddo
       enddo


       do i=1,nn56
          jbox = n56(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             mexps56(j) = mexps56(j) + mexp(j,jbox,4)*zmul
          enddo
       enddo


       do i=1,nsall
          jbox = sall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             mexpnall(j) = mexpnall(j) + mexp(j,jbox,3)*zmul
          enddo
       enddo

       do i=1,ns3478
          jbox = s3478(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             mexpn3478(j) = mexpn3478(j) + mexp(j,jbox,3)*zmul
          enddo
       enddo

       do i=1,ns34
          jbox = s34(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             mexpn34(j) = mexpn34(j) + mexp(j,jbox,3)*zmul
          enddo
       enddo

       do i=1,ns78
          jbox = s78(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             mexpn78(j) = mexpn78(j) + mexp(j,jbox,3)*zmul
          enddo
       enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
       jbox = ichild(1,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpnall(i)
          mexpdownphys(i) = mexpsall(i) + mexps1256(i) + mexps12(i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpnall(i)*ys(1,i)
          mexpdownphys(i) = (mexpsall(i) + mexps1256(i) + mexps12(i))*
     1                        ys(-1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpnall(i)+mexpn34(i)+mexpn3478(i))*
     1         zs(1,i)
          mexpdownphys(i) = mexpsall(i)/zs(1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpnall(i)+mexpn34(i)+mexpn3478(i))*
     1         ys(1,i)*zs(1,i)
          mexpdownphys(i) = mexpsall(i)*ys(-1,i)/zs(1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)
        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpnall(i)*xs(1,i)
          mexpdownphys(i) = (mexpsall(i) + mexps1256(i) + mexps56(i))*
     1                        xs(-1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)
        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpnall(i)*ys(1,i)*xs(1,i)
          mexpdownphys(i) = (mexpsall(i) + mexps1256(i) + mexps56(i))*
     1                        ys(-1,i)*xs(-1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)
        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpnall(i)+mexpn78(i)+mexpn3478(i))*
     1         xs(1,i)*zs(1,i)
          mexpdownphys(i) = mexpsall(i)*xs(-1,i)/zs(1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)
        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpnall(i)+mexpn78(i)+mexpn3478(i))*
     1         ys(1,i)*zs(1,i)*xs(1,i)
          mexpdownphys(i) = mexpsall(i)*ys(-1,i)/zs(1,i)*xs(-1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processewexp(zk2,ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,neall,eall,
     3           ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,ne3,e3,ne5,e5,
     4           ne7,e7,nwall,wall,nw2468,w2468,nw24,w24,nw68,w68,
     5           nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6           mexpup,mexpdown,mexpupphys,mexpdownphys,
     7           mexpeall,mexpe2468,mexpe24,mexpe68,mexpe2,mexpe4,
     8           mexpe6,mexpe8,mexpwall,mexpw1357,mexpw13,mexpw57,
     9           mexpw1,mexpw3,mexpw5,mexpw7,rdminus,
     9           xs,ys,zs,fexpback,rlsc)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer neall,nwall,ne1357,nw2468,ne13,ne57,nw24,nw68
       integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
       integer eall(*),wall(*),e1357(*),w2468(*)
       integer e13(*),e57(*),w24(*),w68(*)
       integer e1(*),e3(*),e5(*),e7(*),w2(*),w4(*),w6(*),w8(*)
       real *8 rscale
       complex *16 zk2
       complex *16 rlams(*),whts(*)
       complex *16 tloc(0:nterms,-nterms:nterms)
       complex *16 tloc2(0:nterms,-nterms:nterms)
       complex *16 mexp(nexptotp,nboxes,6)
       real *8 rdminus(0:nterms,0:nterms,-nterms:nterms)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nexptot),mexpdown(nexptot)
       complex *16 mexpupphys(nexptotp),mexpdownphys(nexptotp)
       complex *16 mexpeall(nexptotp),mexpwall(nexptotp)
       complex *16 mexpw1357(nexptotp),mexpe2468(nexptotp)
       complex *16 mexpw13(nexptotp),mexpw57(nexptotp)
       complex *16 mexpe24(nexptotp),mexpe68(nexptotp)
       complex *16 mexpw1(nexptotp),mexpw3(nexptotp)
       complex *16 mexpw5(nexptotp),mexpw7(nexptotp)
       complex *16 mexpe2(nexptotp),mexpe4(nexptotp)
       complex *16 mexpe6(nexptotp),mexpe8(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp),zs(5,nexptotp)
       complex *16 rlsc(0:nterms,0:nterms,nlams)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j,l
       complex *16 ztmp,zmul
     
       real *8 ctmp(3)


       do i=1,nexptotp

       mexpeall(i) = 0
       mexpwall(i) = 0
       mexpe2468(i) = 0
       mexpe24(i) = 0
       mexpe68(i) = 0
       mexpe2(i) = 0
       mexpe4(i) = 0
       mexpe6(i) = 0
       mexpe8(i) = 0
       mexpw1357(i) = 0
       mexpw13(i) = 0
       mexpw57(i) = 0
       mexpw1(i) = 0
       mexpw3(i) = 0
       mexpw5(i) = 0
       mexpw7(i) = 0

       enddo
      
   
       ctmp(1) = centers(1,ibox) - rscale/2.0d0
       ctmp(2) = centers(2,ibox) - rscale/2.0d0
       ctmp(3) = centers(3,ibox) - rscale/2.0d0
       
       do i=1,neall
          jbox = eall(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpwall(j) = mexpwall(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,ne1357
          jbox = e1357(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw1357(j) = mexpw1357(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,ne13
          jbox = e13(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw13(j) = mexpw13(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo


       do i=1,ne57
          jbox = e57(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw57(j) = mexpw57(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,ne1
          jbox = e1(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw1(j) = mexpw1(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo


       do i=1,ne3
          jbox = e3(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw3(j) = mexpw3(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,ne5
          jbox = e5(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw5(j) = mexpw5(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo


       do i=1,ne7
          jbox = e7(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw7(j) = mexpw7(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,nwall
          jbox = wall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpeall(j) = mexpeall(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw2468
          jbox = w2468(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe2468(j) = mexpe2468(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw24
          jbox = w24(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe24(j) = mexpe24(j) + mexp(j,jbox,5)*zmul
          enddo

       enddo

       

       do i=1,nw68
          jbox = w68(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe68(j) = mexpe68(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw2
          jbox = w2(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe2(j) = mexpe2(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo


       do i=1,nw4
          jbox = w4(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe4(j) = mexpe4(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw6
          jbox = w6(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe6(j) = mexpe6(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw8
          jbox = w8(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe8(j) = mexpe8(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
       jbox = ichild(1,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpeall(i)
          mexpdownphys(i) = mexpwall(i)+mexpw1357(i)+mexpw13(i)+
     1                        mexpw1(i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nterms,tloc,tloc2,rdminus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpeall(i)+mexpe2468(i)+mexpe24(i)+
     1                       mexpe2(i))*zs(1,i)      
          mexpdownphys(i) = mexpwall(i)/zs(1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nterms,tloc,tloc2,rdminus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpeall(i)*ys(1,i)
          mexpdownphys(i) = (mexpwall(i)+mexpw1357(i)+mexpw13(i)+
     1                        mexpw3(i))*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpeall(i)+mexpe2468(i)+mexpe24(i)+
     1                       mexpe4(i))*zs(1,i)*ys(1,i)      
          mexpdownphys(i) = mexpwall(i)/zs(1,i)*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)
        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpeall(i)*xs(-1,i)
          mexpdownphys(i) = (mexpwall(i)+mexpw1357(i)+mexpw57(i)+
     1                        mexpw5(i))*xs(1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)
        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpeall(i)+mexpe2468(i)+mexpe68(i)+
     1                       mexpe6(i))*zs(1,i)*xs(-1,i)      
          mexpdownphys(i) = mexpwall(i)/zs(1,i)*xs(1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)
        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpeall(i)*xs(-1,i)*ys(1,i)
          mexpdownphys(i) = (mexpwall(i)+mexpw1357(i)+mexpw57(i)+
     1                        mexpw7(i))*xs(1,i)*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpeall(i)+mexpe2468(i)+mexpe68(i)+
     1                       mexpe8(i))*zs(1,i)*xs(-1,i)*ys(1,i)      
          mexpdownphys(i) = mexpwall(i)/zs(1,i)*xs(1,i)*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        call h3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      
