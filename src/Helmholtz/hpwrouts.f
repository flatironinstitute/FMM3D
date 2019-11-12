c
cc      plane wave routines for Helmholtz 3D FMM
c 
c***********************************************************************
      subroutine hrlscini(rlsc,nlambs,rlams,rsc,zk,nterms)
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
c       rsc -        scaling parameter for ynm
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
      real *8 rsc,rsctmp

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
      call zylgndrsc(nterms,zmult,rsc,rlsc(0,0,nl))
cc      call zylgndr(nterms,zmult,rlsc(0,0,nl))

      do 1900 iin=1,nterms
      do 1800 im=0,nterms

      rlsc(iin,im,nl) = rlsc(iin,im,nl)*ipow(iin)

 1800 continue
 1900 continue
      

2000  continue
      return
      end
c***********************************************************************
      subroutine hmkexps(rlams,nlambs,numphys,nexptotp,zk,xs,ys,zs)
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
      subroutine hmkfexp(nlambs,numfour,numphys,fexp,fexp2)
      implicit real *8 (a-h,o-z)
      double complex ima
      double complex fexp(*),fexp2(*)
      integer  nlambs,numphys(*),numfour(*)
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

      do i=1,nlambs
        nalpha = numphys(i)
        halpha=2*pi/nalpha

        do j=1,nalpha
          do mm = 1,numfour(i)
            alpha=(j-1)*halpha
            fexp(next) = exp(ima*mm*alpha)
            next = next + 1
          enddo
        enddo
      enddo


      next = 1
      do i=1,nlambs
 
        nalpha = numphys(i)
        halpha=2*pi/nalpha
        do mm = 1,numfour(i)
          do j=1,nalpha
            alpha=(j-1)*halpha
            fexp2(next) = exp(-ima*mm*alpha)
	        next = next + 1
          enddo
        enddo
      enddo


      return
      end
c***********************************************************************

      subroutine hmpoletoexp(nd,mpole,nterms,nlambs,numtets,nexptot,
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
c     nd          in: integer
c                 number of multipole expansions       
c     
c     mpole       in: complex *16 (nd,0:nterms, -nterms:nterms)
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
c     mexpupf     out: complex *16 (nd,nexptot)
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
c     mexpdownf   out: complex *16 (nd,nexptot)
c                 Fourier coefficients of the function 
c                 mexpdown(\lambda,\alpha) for successive
c                 discrete \lambda values
c---------------------------------------------------------------

      implicit none
      integer nd
      integer nterms,nlambs,numtets(nlambs),nexptot
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      complex *16 mexpupf(nd,nexptot)
      complex *16 mexpdownf(nd,nexptot)
      complex *16, allocatable :: ztmp1(:),ztmp2(:),ztmp3(:),ztmp4(:)
      complex *16 zeyep
      complex *16 rlsc(0:nterms,0:nterms,nlambs)

c     Temp variables
      real *8 sgn
      integer ntot,ncurrent,nl,mth,nm,idim

      allocate(ztmp1(nd),ztmp2(nd),ztmp3(nd),ztmp4(nd))

      ntot = 0
      do nl=1,nlambs
         sgn = 1.0d0

         do idim=1,nd
           ztmp1(idim) = 0
           ztmp2(idim) = 0
         enddo

         mth = 0 
         ncurrent = ntot + 1

         do nm = mth,nterms,2
            do idim=1,nd
              ztmp1(idim) = ztmp1(idim) + 
     1            rlsc(nm,mth,nl)*mpole(idim,nm,mth)
            enddo
         enddo

         do nm=mth+1,nterms,2
            do idim=1,nd
              ztmp2(idim) = ztmp2(idim) + 
     1            rlsc(nm,mth,nl)*mpole(idim,nm,mth)
            enddo
         enddo

         do idim=1,nd
           mexpupf(idim,ncurrent) = ztmp1(idim)+ztmp2(idim)
           mexpdownf(idim,ncurrent) = sgn*(ztmp1(idim)-ztmp2(idim))
         enddo

         do mth = 1,numtets(nl)
            ncurrent = ntot + 2*mth

            do idim=1,nd
              ztmp1(idim) = 0
              ztmp2(idim) = 0
              ztmp3(idim) = 0
              ztmp4(idim) = 0
            enddo

            sgn = -sgn
            do nm = mth,nterms,2
               do idim=1,nd
                  ztmp1(idim) = ztmp1(idim) + 
     1                rlsc(nm,mth,nl)*mpole(idim,nm,mth)
                  ztmp3(idim) = ztmp3(idim) + 
     1                rlsc(nm,mth,nl)*mpole(idim,nm,-mth)
               enddo
            enddo

            do nm=mth+1,nterms,2
               do idim=1,nd
                 ztmp2(idim) = ztmp2(idim) + 
     1              rlsc(nm,mth,nl)*mpole(idim,nm,mth)
                 ztmp4(idim) = ztmp4(idim) + 
     1              rlsc(nm,mth,nl)*mpole(idim,nm,-mth)
               enddo
            enddo

            do idim=1,nd
              mexpupf(idim,ncurrent) = ztmp1(idim)+ztmp2(idim)
              mexpdownf(idim,ncurrent) = sgn*(ztmp1(idim)-ztmp2(idim))

              mexpupf(idim,ncurrent+1) = ztmp3(idim)+ztmp4(idim)
              mexpdownf(idim,ncurrent+1) = sgn*(ztmp3(idim)-ztmp4(idim))
            enddo
         enddo
         ntot = ntot + 2*numtets(nl)+1
      enddo

      return
      end

c-----------------------------------------------------------------      
      subroutine hexptolocal(nd,local,nterms,zk,rlambs,whts,nlambs,
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
c     nd               in: integer
c                      number of local expansions
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
c     lexp1f           complex *16(nd,nexptot)
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
c     lexp2f           complex *16(nd,nexptot)
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
c     local(nd,0:nterms,-nterms:nterms): output local expansion of order
c                                     nterms
        
      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot,nthmax,nd
      integer ncurrent,ntot,nl,ncurrent2
      complex *16 local(nd,0:nterms,-nterms:nterms)
      complex *16 ima,zmult
      complex *16 lexp1f(nd,nexptot),lexp2f(nd,nexptot)
      complex *16 rlambs(nlambs), whts(nlambs)
      complex *16 rlsc(0:nterms,0:nterms,nlambs),ztmp
      complex *16 zk,zk2
      real *8 scale, rscale(0:nterms)
    
c     Temporary variables
      integer i, nm, mth, j, mmax,idim
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
            do idim=1,nd
              local(idim,nm,mth) = 0.0d0
            enddo
         enddo
      enddo

      ntot = 0
      do nl=1,nlambs
c        Add contributions to local expansion
         do nm=0,nterms,2
            mmax = numtets(nl)
            if(mmax.gt.nm) mmax = nm

            ncurrent = ntot+1
            ztmp = rlsc(nm,0,nl)*whts(nl)
            do idim=1,nd
               local(idim,nm,0) = local(idim,nm,0)+
     1            (lexp2f(idim,ncurrent)+
     1            lexp1f(idim,ncurrent))*ztmp
            enddo
            do mth=1,mmax
               ncurrent = ntot+2*mth
               ncurrent2 = ntot+2*mth+1

               ztmp = rlsc(nm,mth,nl)*whts(nl)
               do idim=1,nd
                  local(idim,nm,mth) = local(idim,nm,mth)+
     1             (lexp2f(idim,ncurrent)+
     1             lexp1f(idim,ncurrent))*ztmp
                  local(idim,nm,-mth) = local(idim,nm,-mth)+
     1              (lexp2f(idim,ncurrent2)+
     1             lexp1f(idim,ncurrent2))*ztmp
               enddo
            enddo
         enddo
         do nm=1,nterms,2
            mmax = numtets(nl) 
            if(mmax.gt.nm) mmax = nm
            ncurrent = ntot+1

            ztmp = rlsc(nm,0,nl)*whts(nl)
            do idim=1,nd
              local(idim,nm,0) = local(idim,nm,0)+
     1           (lexp2f(idim,ncurrent)-
     1           lexp1f(idim,ncurrent))*ztmp
            enddo
            do mth =1,mmax
               ncurrent = ntot+2*mth
               ncurrent2 = ntot+2*mth+1
               ztmp = rlsc(nm,mth,nl)*whts(nl)
               do idim=1,nd
                 local(idim,nm,mth) = local(idim,nm,mth)+
     1            (lexp2f(idim,ncurrent)-
     1            lexp1f(idim,ncurrent))*ztmp
                 local(idim,nm,-mth) = local(idim,nm,-mth)+
     1             (lexp2f(idim,ncurrent2)-
     1             lexp1f(idim,ncurrent2))*ztmp
               enddo
            enddo
         enddo
         ntot = ntot + 2*numtets(nl)+1
      enddo

      ztmp = zk2
      do nm=0,nterms
         do mth = -nm,nm
              
            do idim=1,nd
              local(idim,nm,mth) = local(idim,nm,mth)*ztmp
            enddo
         enddo
      enddo

      return
      end
c***********************************************************************
      subroutine hphystof(nd,mexpf,nlambs,numfour,numphys,mexpphys,
     1   fexp2)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nd,idim
      complex *16 mexpf(nd,*)
      complex *16 mexpphys(nd,*),ima
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
c     nd:   number of expansions 
c     mexpphys(nd,*):  discrete values of the moment function 
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
c     mexpf(nd,*):     fourier coefficients of the function 
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

      do i=1,nlambs
        nalpha = numphys(i)
        halpha=2*pi/nalpha

        do j=1,nalpha
          alphas(j)=(j-1)*halpha
        enddo

        do idim=1,nd
          mexpf(idim,nftot+1) = 0.0d0
        enddo
        do ival=1,nalpha
           do idim=1,nd
              mexpf(idim,nftot+1) = mexpf(idim,nftot+1) + 
     1            mexpphys(idim,nptot+ival)
           enddo
        enddo

        do idim=1,nd
           mexpf(idim,nftot+1) = mexpf(idim,nftot+1)/nalpha
        enddo
        do mm = 1,numfour(i)

           do idim=1,nd
              mexpf(idim,nftot+2*mm) = 0.0d0
              mexpf(idim,nftot+2*mm+1) = 0.0d0
           enddo

           do ival=1,nalpha
              do idim=1,nd
                 mexpf(idim,nftot+2*mm) = mexpf(idim,nftot+2*mm) +
     1            fexp2(next)*mexpphys(idim,nptot+ival)
                mexpf(idim,nftot+2*mm+1) = mexpf(idim,nftot+2*mm+1) +
     1            conjg(fexp2(next))*mexpphys(idim,nptot+ival)
               enddo
               next = next+1
           enddo
           do idim=1,nd
              mexpf(idim,nftot+2*mm) = mexpf(idim,nftot+2*mm)/nalpha
              mexpf(idim,nftot+2*mm+1) = mexpf(idim,nftot+2*mm+1)/nalpha
           enddo
        enddo
        nftot = nftot+2*numfour(i)+1
        nptot = nptot+numphys(i)
      enddo

      return
      end
c
c********************************************************************
      subroutine hftophys(nd,mexpf,nlambs,numfour,numphys,mexpphys,fexp)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      complex *16 mexpf(nd,*)
      complex *16 mexpphys(nd,*),ima,ctmp
      complex *16 fexp(*)
      real *8  alphas(0:1000)
      integer  nlambs,numfour(nlambs),numphys(nlambs)
      integer nd,idim
      data ima/(0.0d0,1.0d0)/
c***********************************************************************
c
c     this subroutine evaluates the fourier expansion of the
c     exponential moment function m(\lambda,\alpha) at equispaced
c     nodes.
c
c     on input:
c
c     nd      :     number of expansions
c     mexpf(nd,*):     fourier coefficients of the function 
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
c     mexpphys(nd,*):  discrete values of the moment function 
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
      do i=1,nlambs
        do ival=1,numphys(i)
           do idim=1,nd
             mexpphys(idim,nptot+ival) = mexpf(idim,nftot+1)
           enddo
           do mm = 1,numfour(i)
              do idim=1,nd
                ctmp = mexpf(idim,nftot+2*mm)*fexp(next)
                ctmp = ctmp + mexpf(idim,nftot+2*mm+1)*conjg(fexp(next))
                mexpphys(idim,nptot+ival) = mexpphys(idim,nptot+ival) + 
     1              ctmp
              enddo
              next = next + 1
           enddo
         enddo

         nftot = nftot+2*numfour(i)+1
         nptot = nptot+numphys(i)
      enddo


      return
      end
c----------------------------------------------------------      

      subroutine hprocessudexp(nd,zk2,ibox,ilev,nboxes,centers,ichild,
     1           rscale,bs,nterms,iaddr,rmlexp,rlams,whts,nlams,
     2           nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nuall,uall,
     3           nu1234,u1234,ndall,dall,nd5678,d5678,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpuall,mexpu5678,mexpdall,
     5           mexpd1234,xs,ys,zs,fexpback,rlsc)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer nd,idim
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer *8 iaddr(2,nboxes)
       integer ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer nuall,ndall,nu1234,nd5678
       integer uall(*),dall(*),u1234(*),d5678(*)
       real *8 rscale,bs
       complex *16 zk2
       complex *16 rlams(*),whts(*)
       complex *16 tloc(nd,0:nterms,-nterms:nterms)
       complex *16 mexp(nd,nexptotp,nboxes,6)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nd,nexptot),mexpdown(nd,nexptot)
       complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
       complex *16 mexpuall(nd,nexptotp),mexpdall(nd,nexptotp)
       complex *16 mexpd1234(nd,nexptotp),mexpu5678(nd,nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp),zs(5,nexptotp)
       complex *16 rlsc(0:nterms,0:nterms,nlams)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j
       complex *16 ztmp,ztmp2,zmul

     
       real *8 ctmp(3)



       do i=1,nexptotp
         do idim=1,nd
           mexpuall(idim,i) = 0
           mexpdall(idim,i) = 0
           mexpu5678(idim,i) = 0
           mexpd1234(idim,i) = 0
         enddo
       enddo
      
   
       ctmp(1) = centers(1,ibox) - bs/2.0d0
       ctmp(2) = centers(2,ibox) - bs/2.0d0
       ctmp(3) = centers(3,ibox) - bs/2.0d0
       
       do i=1,nuall
          jbox = uall(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
             do idim=1,nd
               mexpdall(idim,j) = mexpdall(idim,j) + 
     1              mexp(idim,j,jbox,2)*zmul
             enddo
          enddo
       enddo

       do i=1,nu1234
          jbox = u1234(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

          do j=1,nexptotp
             zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
             do idim=1,nd
               mexpd1234(idim,j) = mexpd1234(idim,j) + 
     1            mexp(idim,j,jbox,2)*zmul
             enddo
          enddo
       enddo

       do i=1,ndall
          jbox = dall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

         
          do j=1,nexptotp
             zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
             do idim=1,nd
                mexpuall(idim,j) = mexpuall(idim,j) + 
     1              mexp(idim,j,jbox,1)*zmul
             enddo
          enddo
       enddo

       do i=1,nd5678
          jbox = d5678(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
             do idim=1,nd
               mexpu5678(idim,j) = mexpu5678(idim,j) + 
     1            mexp(idim,j,jbox,1)*zmul
             enddo
          enddo
       enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
       jbox = ichild(1,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)
            mexpdownphys(idim,i) = mexpdall(idim,i) + mexpd1234(idim,i)
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          do idim=1,nd
             mexpupphys(idim,i)  = mexpuall(idim,i)*xs(1,i)
             mexpdownphys(idim,i) = (mexpdall(idim,i) + 
     1            mexpd1234(idim,i))*xs(-1,i)
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpdall(idim,i) + 
     1           mexpd1234(idim,i))*ys(-1,i)
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = ys(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)
          do idim=1,nd
             mexpupphys(idim,i)  = mexpuall(idim,i)*ztmp
             mexpdownphys(idim,i) = (mexpdall(idim,i) + 
     1            mexpd1234(idim,i))*ztmp2
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = 1/zs(1,i)

          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1         mexpu5678(idim,i))*zs(1,i)
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = zs(1,i)*xs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1           mexpu5678(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)

          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1          mexpu5678(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)*xs(1,i)
          ztmp2 = xs(-1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1         mexpu5678(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      

      subroutine hprocessnsexp(nd,zk2,ibox,ilev,nboxes,centers,ichild,
     1           rscale,bs,nterms,iaddr,rmlexp,rlams,whts,nlams,
     2           nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nnall,nall,
     3           nn1256,n1256,nn12,n12,nn56,n56,
     4           nsall,sall,ns3478,s3478,ns34,s34,ns78,s78,mexpup,
     5           mexpdown,mexpupphys,mexpdownphys,
     6           mexpnall,mexpn3478,mexpn34,mexpn78,mexpsall,
     7           mexps1256,mexps12,mexps56,rdplus,
     8           xs,ys,zs,fexpback,rlsc)
c--------------------------------------------------------------------
c      process north south expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer nd,idim
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer *8 iaddr(2,nboxes)
       integer ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer nnall,nsall,nn1256,ns3478,nn12,nn56,ns34,ns78
       integer nall(*),sall(*),n1256(*),s3478(*)
       integer n12(*),n56(*),s34(*),s78(*)
       real *8 rscale,bs
       complex *16 zk2
       complex *16 rlams(*),whts(*)
       complex *16 tloc(nd,0:nterms,-nterms:nterms)
       complex *16, allocatable :: tloc2(:,:,:)
       complex *16 mexp(nd,nexptotp,nboxes,6)
       real *8 rdplus(0:nterms,0:nterms,-nterms:nterms)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nd,nexptot),mexpdown(nd,nexptot)
       complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
       complex *16 mexpnall(nd,nexptotp),mexpsall(nd,nexptotp)
       complex *16 mexps1256(nd,nexptotp),mexpn3478(nd,nexptotp)
       complex *16 mexps12(nd,nexptotp),mexps56(nd,nexptotp)
       complex *16 mexpn34(nd,nexptotp),mexpn78(nd,nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp),zs(5,nexptotp)
       complex *16 rlsc(0:nterms,0:nterms,nlams)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j
       complex *16 ztmp,zmul,ztmp2

     
       real *8 ctmp(3)
       allocate(tloc2(nd,0:nterms,-nterms:nterms))
       

       do i=1,nexptotp
         do idim=1,nd
           mexpnall(idim,i) = 0
           mexpsall(idim,i) = 0
           mexpn3478(idim,i) = 0
           mexpn34(idim,i) = 0
           mexpn78(idim,i) = 0
           mexps1256(idim,i) = 0
           mexps12(idim,i) = 0
           mexps56(idim,i) = 0
         enddo
       enddo
      
   
       ctmp(1) = centers(1,ibox) - bs/2.0d0
       ctmp(2) = centers(2,ibox) - bs/2.0d0
       ctmp(3) = centers(3,ibox) - bs/2.0d0
       
       do i=1,nnall
          jbox = nall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             do idim=1,nd
               mexpsall(idim,j) = mexpsall(idim,j) + 
     1            mexp(idim,j,jbox,4)*zmul
             enddo
          enddo

       enddo

       do i=1,nn1256
          jbox = n1256(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             do idim=1,nd
               mexps1256(idim,j) = mexps1256(idim,j) + 
     1             mexp(idim,j,jbox,4)*zmul
             enddo
          enddo
       enddo

       do i=1,nn12
          jbox = n12(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             do idim=1,nd
               mexps12(idim,j) = mexps12(idim,j) + 
     1           mexp(idim,j,jbox,4)*zmul
             enddo
          enddo
       enddo


       do i=1,nn56
          jbox = n56(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             do idim=1,nd
               mexps56(idim,j) = mexps56(idim,j) + 
     1             mexp(idim,j,jbox,4)*zmul
             enddo
          enddo
       enddo


       do i=1,nsall
          jbox = sall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             do idim=1,nd
               mexpnall(idim,j) = mexpnall(idim,j) + 
     1            mexp(idim,j,jbox,3)*zmul
             enddo
          enddo
       enddo

       do i=1,ns3478
          jbox = s3478(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             do idim=1,nd
               mexpn3478(idim,j) = mexpn3478(idim,j) + 
     1             mexp(idim,j,jbox,3)*zmul
             enddo
          enddo
       enddo

       do i=1,ns34
          jbox = s34(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             do idim=1,nd
               mexpn34(idim,j) = mexpn34(idim,j) + 
     1            mexp(idim,j,jbox,3)*zmul
             enddo
          enddo
       enddo

       do i=1,ns78
          jbox = s78(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             do idim=1,nd
               mexpn78(idim,j) = mexpn78(idim,j) + 
     1            mexp(idim,j,jbox,3)*zmul
             enddo
          enddo
       enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
       jbox = ichild(1,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)
            mexpdownphys(idim,i) = mexpsall(idim,i) + 
     1         mexps1256(idim,i) + mexps12(idim,i)
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpsall(idim,i) + 
     1           mexps1256(idim,i) + mexps12(idim,i))*ys(-1,i)      
         enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = 1/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+
     1         mexpn34(idim,i)+mexpn3478(idim,i))*zs(1,i)
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp      
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+
     1         mexpn34(idim,i)+mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp

          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*xs(1,i)
            mexpdownphys(idim,i) = (mexpsall(idim,i) + 
     1          mexps1256(idim,i) + mexps56(idim,i))*xs(-1,i)      
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
         ztmp = ys(1,i)*xs(1,i)
         ztmp2 = ys(-1,i)*xs(-1,i)
         do idim=1,nd
           mexpupphys(idim,i)  = mexpnall(idim,i)*ztmp
           mexpdownphys(idim,i) = (mexpsall(idim,i) + 
     1        mexps1256(idim,i) + mexps56(idim,i))*ztmp2
         enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
         ztmp = xs(1,i)*zs(1,i)
         ztmp2 = xs(-1,i)/zs(1,i)
         do idim=1,nd
           mexpupphys(idim,i)  = (mexpnall(idim,i)+
     1        mexpn78(idim,i)+mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
         enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
         ztmp = ys(1,i)*zs(1,i)*xs(1,i)
         ztmp2 = ys(-1,i)/zs(1,i)*xs(-1,i)
         do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+
     1         mexpn78(idim,i)+mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2
         enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      

      subroutine hprocessewexp(nd,zk2,ibox,ilev,nboxes,centers,ichild,
     1           rscale,bs,nterms,iaddr,rmlexp,rlams,whts,nlams,
     2           nfourier,
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
c      process east west expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer nd
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer *8 iaddr(2,nboxes)
       integer ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer neall,nwall,ne1357,nw2468,ne13,ne57,nw24,nw68
       integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
       integer eall(*),wall(*),e1357(*),w2468(*)
       integer e13(*),e57(*),w24(*),w68(*)
       integer e1(*),e3(*),e5(*),e7(*),w2(*),w4(*),w6(*),w8(*)
       real *8 rscale,bs
       complex *16 zk2
       complex *16 rlams(*),whts(*)
       complex *16 tloc(nd,0:nterms,-nterms:nterms)
       complex *16, allocatable :: tloc2(:,:,:)
       complex *16 mexp(nd,nexptotp,nboxes,6)
       real *8 rdminus(0:nterms,0:nterms,-nterms:nterms)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nd,nexptot),mexpdown(nd,nexptot)
       complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
       complex *16 mexpeall(nd,nexptotp),mexpwall(nd,nexptotp)
       complex *16 mexpw1357(nd,nexptotp),mexpe2468(nd,nexptotp)
       complex *16 mexpw13(nd,nexptotp),mexpw57(nd,nexptotp)
       complex *16 mexpe24(nd,nexptotp),mexpe68(nd,nexptotp)
       complex *16 mexpw1(nd,nexptotp),mexpw3(nd,nexptotp)
       complex *16 mexpw5(nd,nexptotp),mexpw7(nd,nexptotp)
       complex *16 mexpe2(nd,nexptotp),mexpe4(nd,nexptotp)
       complex *16 mexpe6(nd,nexptotp),mexpe8(nd,nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp),zs(5,nexptotp)
       complex *16 rlsc(0:nterms,0:nterms,nlams)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j,l,idim
       complex *16 ztmp,zmul,ztmp2

     
       real *8 ctmp(3)

       allocate(tloc2(nd,0:nterms,-nterms:nterms))


       do i=1,nexptotp
         do idim=1,nd
           mexpeall(idim,i) = 0
           mexpwall(idim,i) = 0
           mexpe2468(idim,i) = 0
           mexpe24(idim,i) = 0
           mexpe68(idim,i) = 0
           mexpe2(idim,i) = 0
           mexpe4(idim,i) = 0
           mexpe6(idim,i) = 0
           mexpe8(idim,i) = 0
           mexpw1357(idim,i) = 0
           mexpw13(idim,i) = 0
           mexpw57(idim,i) = 0
           mexpw1(idim,i) = 0
           mexpw3(idim,i) = 0
           mexpw5(idim,i) = 0
           mexpw7(idim,i) = 0
         enddo
       enddo
      
   
       ctmp(1) = centers(1,ibox) - bs/2.0d0
       ctmp(2) = centers(2,ibox) - bs/2.0d0
       ctmp(3) = centers(3,ibox) - bs/2.0d0
       
       do i=1,neall
          jbox = eall(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             do idim=1,nd
               mexpwall(idim,j) = mexpwall(idim,j) + 
     1             mexp(idim,j,jbox,6)*zmul
             enddo
          enddo
       enddo

       do i=1,ne1357
          jbox = e1357(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             do idim=1,nd
               mexpw1357(idim,j) = mexpw1357(idim,j) + 
     1             mexp(idim,j,jbox,6)*zmul
             enddo
          enddo
       enddo

       do i=1,ne13
          jbox = e13(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             do idim=1,nd
               mexpw13(idim,j) = mexpw13(idim,j) + 
     1            mexp(idim,j,jbox,6)*zmul
             enddo
          enddo
       enddo


       do i=1,ne57
          jbox = e57(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             do idim=1,nd
               mexpw57(idim,j) = mexpw57(idim,j) + 
     1            mexp(idim,j,jbox,6)*zmul
             enddo
          enddo
       enddo

       do i=1,ne1
          jbox = e1(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             do idim=1,nd
               mexpw1(idim,j) = mexpw1(idim,j) + 
     1            mexp(idim,j,jbox,6)*zmul
             enddo
          enddo
       enddo


       do i=1,ne3
          jbox = e3(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             do idim=1,nd
               mexpw3(idim,j) = mexpw3(idim,j) + 
     1            mexp(idim,j,jbox,6)*zmul
             enddo
          enddo
       enddo

       do i=1,ne5
          jbox = e5(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             do idim=1,nd
               mexpw5(idim,j) = mexpw5(idim,j) + 
     1            mexp(idim,j,jbox,6)*zmul
             enddo
          enddo
       enddo


       do i=1,ne7
          jbox = e7(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             do idim=1,nd
               mexpw7(idim,j) = mexpw7(idim,j) + 
     1            mexp(idim,j,jbox,6)*zmul
             enddo
          enddo
       enddo

       do i=1,nwall
          jbox = wall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             do idim=1,nd
               mexpeall(idim,j) = mexpeall(idim,j) + 
     1             mexp(idim,j,jbox,5)*zmul
             enddo
          enddo
       enddo

       do i=1,nw2468
          jbox = w2468(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             do idim=1,nd
               mexpe2468(idim,j) = mexpe2468(idim,j) + 
     1            mexp(idim,j,jbox,5)*zmul
             enddo
          enddo
       enddo

       do i=1,nw24
          jbox = w24(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             do idim=1,nd
               mexpe24(idim,j) = mexpe24(idim,j) + 
     1             mexp(idim,j,jbox,5)*zmul
             enddo
          enddo

       enddo

       

       do i=1,nw68
          jbox = w68(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             do idim=1,nd
               mexpe68(idim,j) = mexpe68(idim,j) + 
     1             mexp(idim,j,jbox,5)*zmul
             enddo
          enddo
       enddo

       do i=1,nw2
          jbox = w2(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             do idim=1,nd
               mexpe2(idim,j) = mexpe2(idim,j) + 
     1            mexp(idim,j,jbox,5)*zmul
             enddo
          enddo
       enddo


       do i=1,nw4
          jbox = w4(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             do idim=1,nd
               mexpe4(idim,j) = mexpe4(idim,j) + 
     1            mexp(idim,j,jbox,5)*zmul
             enddo
          enddo
       enddo

       do i=1,nw6
          jbox = w6(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             do idim=1,nd
               mexpe6(idim,j) = mexpe6(idim,j) + 
     1            mexp(idim,j,jbox,5)*zmul
             enddo
          enddo
       enddo

       do i=1,nw8
          jbox = w8(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             do idim=1,nd
               mexpe8(idim,j) = mexpe8(idim,j) + 
     1            mexp(idim,j,jbox,5)*zmul
             enddo
          enddo
       enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
       jbox = ichild(1,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          do idim=1,nd
             mexpupphys(idim,i)  = mexpeall(idim,i)
             mexpdownphys(idim,i) = mexpwall(idim,i)+mexpw1357(idim,i)+
     1                mexpw13(idim,i)+mexpw1(idim,i)
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = 1/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1          mexpe24(idim,i)+mexpe2(idim,i))*zs(1,i)      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpwall(idim,i)+mexpw1357(idim,i)+
     1                mexpw13(idim,i)+mexpw3(idim,i))*ys(-1,i)
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1               mexpe24(idim,i)+mexpe4(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*xs(-1,i)
            mexpdownphys(idim,i) = (mexpwall(idim,i)+mexpw1357(idim,i)+
     1            mexpw57(idim,i)+mexpw5(idim,i))*xs(1,i)
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = zs(1,i)*xs(-1,i)
          ztmp2 = xs(1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1              mexpe68(idim,i)+mexpe6(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = xs(-1,i)*ys(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)

          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*ztmp
            mexpdownphys(idim,i) = (mexpwall(idim,i)+mexpw1357(idim,i)+
     1            mexpw57(idim,i)+mexpw7(idim,i))*ztmp2
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          ztmp = zs(1,i)*xs(-1,i)*ys(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1            mexpe68(idim,i)+mexpe8(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
       enddo

       call hphystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call hphystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call hexptolocal(nd,tloc,nterms,zk2,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      
c
c
c--------------------------------------------------------------------
      subroutine hprocesslist3udexp(nd,ibox,nboxes,centers,
     1           bs,nterms,nexptotp,mexp,nuall,uall,
     3           ndall,dall,mexpuall,mexpdall,
     5           xs,ys,zs)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer idim,nd
      integer ibox,nboxes,nterms,nlams,nthmax
      integer nexptot,nexptotp
      integer nuall,ndall
      integer uall(*),dall(*)
      double precision bs
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision centers(3,nboxes)
      double complex mexpuall(nd,nexptotp),mexpdall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double complex zs(5,nexptotp)

c      temp variables
      integer jbox,i,ix,iy,iz,j
      double precision rtmp
      double complex ztmp,zmul,ztmp2
     
      double precision ctmp(3)


      do i=1,nexptotp
        do idim=1,nd
          mexpuall(idim,i) = 0
          mexpdall(idim,i) = 0
        enddo
      enddo
      
   
      ctmp(1) = centers(1,ibox) - bs/2.0d0
      ctmp(2) = centers(2,ibox) - bs/2.0d0
      ctmp(3) = centers(3,ibox) - bs/2.0d0
  
      
      do i=1,nuall
        jbox = uall(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
          do idim=1,nd
            mexpdall(idim,j) = mexpdall(idim,j) + 
     1                         mexp(idim,j,jbox,2)*zmul
          enddo
        enddo
      enddo

      do i=1,ndall
        jbox = dall(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

        do j=1,nexptotp
          zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
          do idim=1,nd
            mexpuall(idim,j) = mexpuall(idim,j) + 
     1                         mexp(idim,j,jbox,1)*zmul
          enddo
        enddo
      enddo


      return
      end
c--------------------------------------------------------------------      

      subroutine hprocesslist3nsexp(nd,ibox,nboxes,centers,
     1           bs,nterms,nexptotp,mexp,nnall,nall,
     3           nsall,sall,mexpnall,mexpsall,
     5           xs,ys,zs)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer nd
      integer ibox,nboxes,nterms,nlams,nthmax
      integer nexptotp
      integer nnall,nsall
      integer nall(*),sall(*)
      double precision bs
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision centers(3,*)
      double complex mexpnall(nd,nexptotp),mexpsall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double complex zs(5,nexptotp)

c      temp variables
      integer jbox,i,ix,iy,iz,j,idim
      double complex ztmp,zmul,ztmp2
      double precision rtmp
    
      double precision ctmp(3)


      do i=1,nexptotp
        do idim=1,nd
          mexpnall(idim,i) = 0
          mexpsall(idim,i) = 0
        enddo
      enddo
      
   
      ctmp(1) = centers(1,ibox) - bs/2.0d0
      ctmp(2) = centers(2,ibox) - bs/2.0d0
      ctmp(3) = centers(3,ibox) - bs/2.0d0
       
      do i=1,nnall
        jbox = nall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
           zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
           do idim=1,nd
             mexpsall(idim,j) = mexpsall(idim,j) + 
     1                          mexp(idim,j,jbox,4)*zmul
           enddo
        enddo

      enddo

      do i=1,nsall
        jbox = sall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
          do idim=1,nd
            mexpnall(idim,j) = mexpnall(idim,j) + 
     1                         mexp(idim,j,jbox,3)*zmul
          enddo
        enddo
      enddo

      return
      end
c--------------------------------------------------------------------      

      subroutine hprocesslist3ewexp(nd,ibox,nboxes,centers,
     1           bs,nterms,nexptotp,mexp,neall,eall,nwall,wall,
     4           mexpeall,mexpwall,xs,ys,zs)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer nd
      integer ibox,nboxes,nterms,nlams,nthmax
      integer nexptotp
      integer neall,nwall
      integer eall(*),wall(*)
      double precision bs
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision centers(3,*)
      double complex mexpeall(nd,nexptotp),mexpwall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double complex zs(5,nexptotp)

c      temp variables
      integer jbox,i,ix,iy,iz,j,l,idim
      double complex ztmp,zmul,ztmp2
      double complex rtmp
     
      double precision ctmp(3)


      do i=1,nexptotp
        do idim=1,nd
          mexpeall(idim,i) = 0
          mexpwall(idim,i) = 0
        enddo
      enddo
      
   
      ctmp(1) = centers(1,ibox) - bs/2.0d0
      ctmp(2) = centers(2,ibox) - bs/2.0d0
      ctmp(3) = centers(3,ibox) - bs/2.0d0
       
      do i=1,neall
        jbox = eall(i)
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs
         
        do j=1,nexptotp
          zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
          do idim=1,nd
            mexpwall(idim,j) = mexpwall(idim,j) + 
     1                         mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo

      do i=1,nwall
        jbox = wall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/bs
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/bs
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/bs

         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpeall(idim,j) = mexpeall(idim,j) + 
     1                         mexp(idim,j,jbox,5)*zmul
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
c--------------------------------------------------------------------     

      subroutine hpw_ud_eval_p(nd,zk2,center,boxsize,ntarg,targ,nlam,
     1   rlams,
     1   whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot)
      implicit none
      integer nd
      real *8 center(3),boxsize,targ(3,ntarg)
      complex *16 rlams(nlam),pot(nd,ntarg)
      complex *16 whts(nlam),zk2
      integer ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),cc2(:)
      integer itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr
      complex *16 rz,zsc,rk
      complex *16, allocatable :: zexp(:),zexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      zsc = -ima/zk2

      allocate(zexp(nlam),zexpinv(nlam),cc(nphmax),cc2(nphmax))



      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize


        do i=1,nlam
          zexp(i) = exp(-z*rlams(i))*whts(i)
          zexpinv(i) = exp(z*rlams(i))*whts(i)
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          rk = sqrt(rlams(il)**2 + zk2**2)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            cc(iphys) = exp(ima*rk*(x*cos(alpha) + y*sin(alpha)))
            cc2(iphys) = exp(-ima*rk*(x*cos(alpha) + y*sin(alpha)))
          enddo

          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz = 
     1            (mexpupphys(idim,ii)*zexp(il)*cc(iphys) + 
     2             mexpdownphys(idim,ii)*zexpinv(il)*cc2(iphys))*hh
              pot(idim,itarg) = pot(idim,itarg) + rz*zsc 
            enddo
          enddo
          istart = istart + nphys(il)
        enddo
      enddo

      return
      end
c
c
c
c
c
c--------------------------------------------------------------------     

      subroutine hpw_ns_eval_p(nd,zk2,center,boxsize,ntarg,targ,nlam,
     1   rlams,whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot)
      implicit none
      integer nd
      real *8 center(3),boxsize,targ(3,ntarg)
      complex *16 rlams(nlam),pot(nd,ntarg)
      complex *16 whts(nlam),zk2
      integer ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),cc2(:)
      integer itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr
      complex *16 rz,zsc,rk
      complex *16, allocatable :: zexp(:),zexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(zexp(nlam),zexpinv(nlam),cc(nphmax),cc2(nphmax))

      zsc = -ima/zk2


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

        do i=1,nlam
          zexp(i) = exp(-y*rlams(i))*whts(i)
          zexpinv(i) = exp(y*rlams(i))*whts(i)
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          rk = sqrt(rlams(il)**2 + zk2**2)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            cc(iphys) = exp(ima*rk*(z*cos(alpha) + x*sin(alpha)))
            cc2(iphys) = exp(-ima*rk*(z*cos(alpha) + x*sin(alpha)))
          enddo


          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz = 
     1            (mexpupphys(idim,ii)*zexp(il)*cc(iphys) + 
     2             mexpdownphys(idim,ii)*zexpinv(il)*cc2(iphys))*hh
              pot(idim,itarg) = pot(idim,itarg) + rz*zsc 
            enddo
          enddo
          istart = istart + nphys(il)

        enddo
      enddo

      return
      end
c
c
c
c
c
c--------------------------------------------------------------------     

      subroutine hpw_ew_eval_p(nd,zk2,center,boxsize,ntarg,targ,nlam,
     1   rlams,whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot)
      implicit none
      integer nd
      real *8 center(3),boxsize,targ(3,ntarg)
      complex *16 rlams(nlam),pot(nd,ntarg)
      complex *16 whts(nlam),zk2
      integer ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),cc2(:)
      integer itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr
      complex *16 rz,zsc,rk
      complex *16, allocatable :: zexp(:),zexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(zexp(nlam),zexpinv(nlam),cc(nphmax),cc2(nphmax))

      zsc = -ima/zk2

      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize


        do i=1,nlam
          zexp(i) = exp(-x*rlams(i))*whts(i)
          zexpinv(i) = exp(x*rlams(i))*whts(i)
        enddo



        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          rk = sqrt(rlams(il)**2 + zk2**2)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            cc(iphys) = exp(ima*rk*(-z*cos(alpha) + y*sin(alpha)))
            cc2(iphys) = exp(ima*rk*(z*cos(alpha) - y*sin(alpha)))
          enddo

          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz = 
     1            (mexpupphys(idim,ii)*zexp(il)*cc(iphys) + 
     2             mexpdownphys(idim,ii)*zexpinv(il)*cc2(iphys))*hh
              pot(idim,itarg) = pot(idim,itarg) + rz*zsc 
            enddo
          enddo
          istart = istart + nphys(il)

        enddo
      enddo

      return
      end
c
c
c
c
c
c--------------------------------------------------------------------     

      subroutine hpw_ud_eval_g(nd,zk2,center,boxsize,ntarg,targ,nlam,
     1   rlams,whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot,
     2   grad)
      implicit none
      integer nd
      real *8 center(3),boxsize,targ(3,ntarg)
      complex *16 rlams(nlam),pot(nd,ntarg)
      complex *16 grad(nd,3,ntarg),whts(nlam),zk2
      integer ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),crc(:),crs(:),cc2(:)
      integer itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr,binv
      complex *16 rz,rz1,rz2,zsc,rk
      complex *16, allocatable :: zexp(:),zexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(zexp(nlam),zexpinv(nlam),cc(nphmax),cc2(nphmax))
      allocate(crc(nphmax),crs(nphmax))

      binv = 1.0d0/boxsize
      zsc = -ima/zk2


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

        do i=1,nlam
          zexp(i) = exp(-z*rlams(i))*whts(i)
          zexpinv(i) = exp(z*rlams(i))*whts(i)
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          rk = sqrt(rlams(il)**2 + zk2**2)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            crc(iphys) = rk*cos(alpha)*ima
            crs(iphys) = rk*sin(alpha)*ima
            cc(iphys) = exp(crc(iphys)*x + crs(iphys)*y)
            cc2(iphys) = exp(-(crc(iphys)*x + crs(iphys)*y))
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz1 = mexpupphys(idim,ii)*zexp(il)*cc(iphys)*hh*zsc
              rz2 = mexpdownphys(idim,ii)*zexpinv(il)*cc2(iphys)*hh*zsc
              rz = rz1 + rz2
              pot(idim,itarg) = pot(idim,itarg) + rz
              grad(idim,1,itarg) = grad(idim,1,itarg) + 
     1            (rz1-rz2)*crc(iphys)*binv
              grad(idim,2,itarg) = grad(idim,2,itarg) + 
     1            (rz1-rz2)*crs(iphys)*binv
              grad(idim,3,itarg) = grad(idim,3,itarg) - 
     1            (rz1-rz2)*binv*rlams(il)
            enddo
          enddo
          istart = istart + nphys(il)
        enddo
      enddo

      return
      end
c
c
c
c
c
c
c
c
c--------------------------------------------------------------------     

      subroutine hpw_ns_eval_g(nd,zk2,center,boxsize,ntarg,targ,nlam,
     1  rlams,whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot,
     2  grad)
      implicit none
      integer nd
      real *8 center(3),boxsize,targ(3,ntarg)
      complex *16 rlams(nlam),pot(nd,ntarg)
      complex *16 grad(nd,3,ntarg),whts(nlam),zk2
      integer ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),crc(:),crs(:),cc2(:)
      integer itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr,binv
      complex *16 rz,rz1,rz2,zsc,rk
      complex *16, allocatable :: zexp(:),zexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(zexp(nlam),zexpinv(nlam),cc(nphmax),cc2(nphmax))
      allocate(crc(nphmax),crs(nphmax))

      binv = 1.0d0/boxsize

      zsc = -ima/zk2

      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

        do i=1,nlam
          zexp(i) = exp(-y*rlams(i))*whts(i)
          zexpinv(i) = exp(y*rlams(i))*whts(i)
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          rk = sqrt(rlams(il)**2 + zk2**2) 

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            crc(iphys) = rk*cos(alpha)*ima
            crs(iphys) = rk*sin(alpha)*ima
            cc(iphys) = exp(crc(iphys)*z + crs(iphys)*x)
            cc2(iphys) = exp(-(crc(iphys)*z + crs(iphys)*x))
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz1 = mexpupphys(idim,ii)*zexp(il)*cc(iphys)*hh*zsc
              rz2 = mexpdownphys(idim,ii)*zexpinv(il)*cc2(iphys)*hh*zsc
              rz = rz1 + rz2
              pot(idim,itarg) = pot(idim,itarg) + rz
              grad(idim,1,itarg) = grad(idim,1,itarg) + 
     1            (rz1-rz2)*crs(iphys)*binv
              grad(idim,2,itarg) = grad(idim,2,itarg) - 
     1            (rz1-rz2)*binv*rlams(il)
              grad(idim,3,itarg) = grad(idim,3,itarg) + 
     1            (rz1-rz2)*crc(iphys)*binv

            enddo
          enddo
          istart = istart + nphys(il)
        enddo
      enddo

      return
      end
c
c
c
c
c
c--------------------------------------------------------------------     

      subroutine hpw_ew_eval_g(nd,zk2,center,boxsize,ntarg,targ,nlam,
     1  rlams,whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot,
     2  grad)
      implicit none
      integer nd
      real *8 center(3),boxsize,targ(3,ntarg)
      complex *16 rlams(nlam),pot(nd,ntarg)
      complex *16 grad(nd,3,ntarg),whts(nlam),zk2
      integer ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),crc(:),crs(:),cc2(:)
      integer itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr,binv
      complex *16 rz,rz1,rz2,rk,zsc
      complex *16, allocatable :: zexp(:),zexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(zexp(nlam),zexpinv(nlam),cc(nphmax),cc2(nphmax))
      allocate(crc(nphmax),crs(nphmax))

      binv = 1.0d0/boxsize
      zsc = -ima/zk2


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

        do i=1,nlam
          zexp(i) = exp(-x*rlams(i))*whts(i)
          zexpinv(i) = exp(x*rlams(i))*whts(i)
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)
          
          rk = sqrt(rlams(il)**2 + zk2**2)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            crc(iphys) = rk*cos(alpha)*ima
            crs(iphys) = rk*sin(alpha)*ima
            cc(iphys) = exp(-crc(iphys)*z + crs(iphys)*y)
            cc2(iphys) = exp(crc(iphys)*z - crs(iphys)*y)
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz1 = mexpupphys(idim,ii)*zexp(il)*cc(iphys)*hh*zsc
              rz2 = mexpdownphys(idim,ii)*zexpinv(il)*cc2(iphys)*hh*zsc
              rz = rz1 + rz2
              pot(idim,itarg) = pot(idim,itarg) + rz
              grad(idim,1,itarg) = grad(idim,1,itarg) - 
     1            (rz1-rz2)*binv*rlams(il)
              grad(idim,2,itarg) = grad(idim,2,itarg) + 
     1            (rz1-rz2)*crs(iphys)*binv
              grad(idim,3,itarg) = grad(idim,3,itarg) - 
     1            (rz1-rz2)*crc(iphys)*binv

            enddo
          enddo
          istart = istart + nphys(il)
        enddo
      enddo

      return
      end
c
c
c
c
c
