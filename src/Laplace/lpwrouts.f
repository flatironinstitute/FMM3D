c
cc     Plane wave routines for Laplace 3D FMM
c-------------------------------------------------------------

      subroutine rlscini(rlsc,nlambs,rlams,nterms)
      implicit double precision (a-h,o-z)
      double precision rlsc(0:nterms,0:nterms,nlambs)
      double precision     rlams(nlambs),rlampow(0:100)
      double precision     facts(0:200)
c
      facts(0) = 1.0d0
      do i = 1,100
	    facts(i) = facts(i-1)*dsqrt(i+0.0d0)
      enddo
c
      do nl = 1,nlambs
c
c     compute powers of lambda_nl
c
        rlampow(0) = 1.0d0
        rmul = rlams(nl)
        do j = 1,nterms
          rlampow(j) = rlampow(j-1)*rmul
        enddo
        do j = 0,nterms
          do k = 0,j
            rlsc(j,k,nl) = rlampow(j)/(facts(j-k)*facts(j+k))
          enddo
        enddo
      enddo

      return
      end
c-------------------------------------------------------------
      subroutine mkexps(rlams,nlambs,numphys,nexptotp,xs,ys,zs)
      implicit double precision (a-h,o-z)
      integer(8) nlambs,nexptotp
      double complex ima
      double complex xs(-5:5,nexptotp)
      double complex ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision     rlams(nlambs),u
      integer(8) numphys(nlambs)
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
c     rlams(nlambs)   discretization points in lambda integral 
c     nlambs          number of discret. pts. in lambda integral
c     numphys(j)     number of nodes in u integral needed 
c                    for corresponding lambda =  lambda_j. 
c     nexptotp        sum_j numphys(j)
c
c     on output:
c
c     xs(m,nexptotp)   e^{i \lambda_j (m cos(u_k)}  in above ordering
c                         for m=-5,-4,-3,-2,-1,0,1,2,3,4,5
c     ys(m,nexptotp)   e^{i \lambda_j (m sin(u_k)}  in above ordering,
c                         for m=-5,-4,-3,-2,-1,0,1,2,3,4,5
c     zs(1,nexptotp)   e^{- m \lambda_j}     in above ordering,
c                         for m=1,2,3,4,5
c------------------------------------------------------------
c      
c     loop over each lambda value 
c
      pi = 4*datan(1.0d0)
      ntot = 0
      do nl = 1,nlambs
        hu=2*pi/numphys(nl)
        do mth = 1,numphys(nl)
          u = (mth-1)*hu
          ncurrent = ntot+mth
          zs(1,ncurrent) = exp( -rlams(nl) )
          zs(2,ncurrent) = exp( - 2.0d0*rlams(nl) )
          zs(3,ncurrent) = exp( - 3.0d0*rlams(nl) )
          zs(4,ncurrent) = exp( - 4.0d0*rlams(nl) )
          zs(5,ncurrent) = exp( - 5.0d0*rlams(nl) )
          xs(1,ncurrent) = exp(ima*rlams(nl)*cos(u))
          xs(2,ncurrent) = exp(ima*rlams(nl)*2.0d0*cos(u))
          xs(3,ncurrent) = exp(ima*rlams(nl)*3.0d0*cos(u))
          xs(4,ncurrent) = exp(ima*rlams(nl)*4.0d0*cos(u))
          xs(5,ncurrent) = exp(ima*rlams(nl)*5.0d0*cos(u))
          ys(1,ncurrent) = exp(ima*rlams(nl)*sin(u))
          ys(2,ncurrent) = exp(ima*rlams(nl)*2.0d0*sin(u))
          ys(3,ncurrent) = exp(ima*rlams(nl)*3.0d0*sin(u))
          ys(4,ncurrent) = exp(ima*rlams(nl)*4.0d0*sin(u))
          ys(5,ncurrent) = exp(ima*rlams(nl)*5.0d0*sin(u))

          xs(0,ncurrent)  = 1.0d0
          xs(-1,ncurrent) = exp(-ima*rlams(nl)*cos(u))
          xs(-2,ncurrent) = exp(-ima*rlams(nl)*2.0d0*cos(u))
          xs(-3,ncurrent) = exp(-ima*rlams(nl)*3.0d0*cos(u))
          xs(-4,ncurrent) = exp(-ima*rlams(nl)*4.0d0*cos(u))
          xs(-5,ncurrent) = exp(-ima*rlams(nl)*5.0d0*cos(u))

          ys(0,ncurrent) = 1.0d0
          ys(-1,ncurrent) = exp(-ima*rlams(nl)*sin(u))
          ys(-2,ncurrent) = exp(-ima*rlams(nl)*2.0d0*sin(u))
          ys(-3,ncurrent) = exp(-ima*rlams(nl)*3.0d0*sin(u))
          ys(-4,ncurrent) = exp(-ima*rlams(nl)*4.0d0*sin(u))
          ys(-5,ncurrent) = exp(-ima*rlams(nl)*5.0d0*sin(u))
        enddo
        ntot = ntot + numphys(nl)
      enddo

      return
      end
c***********************************************************************
      subroutine mkfexp(nlambs,numfour,numphys,fexpe,fexpo,fexpback)
      implicit double precision (a-h,o-z)
      double complex ima
      double complex fexpe(1)
      double complex fexpo(1)
      double complex fexpback(1)
      integer(8)  nlambs,numphys(nlambs),numfour(nlambs)
      data ima/(0.0d0,1.0d0)/
c
c     this subroutine computes the tables of exponentials needed
c     for mapping from fourier to physical domain. 
c     in order to minimize storage, they are organized in a 
c     one-dimenional array corresponding to the order in which they
c     are accessed by subroutine ftophys.
c    
c     size of fexpe, fexpo =          40000   for nlambs = 39
c     size of fexpe, fexpo =          15000   for nlambs = 30
c     size of fexpe, fexpo =           4000   for nlambs = 20
c     size of fexpe, fexpo =            400   for nlambs = 10
c
c***********************************************************************
      pi = 4*datan(1.0d0)
      nexte = 1
      nexto = 1
      do i=1,nlambs
	    nalpha = numphys(i)
        halpha=2*pi/nalpha
        do j=1,nalpha
          alpha=(j-1)*halpha
	      do mm = 2,numfour(i),2
            fexpe(nexte)  = cdexp(ima*(mm-1)*alpha)
	        nexte = nexte + 1
          enddo
	      do mm = 3,numfour(i),2
            fexpo(nexto)  = cdexp(ima*(mm-1)*alpha)
	        nexto = nexto + 1
          enddo
        enddo
      enddo

      next = 1
      do i=1,nlambs
	    nalpha = numphys(i)
        halpha=2*pi/nalpha
	    do mm = 2,numfour(i)
          do j=1,nalpha
            alpha=(j-1)*halpha
            fexpback(next)  = cdexp(-ima*(mm-1)*alpha)
	        next = next + 1
          enddo
        enddo
      enddo

      return
      end

c---------------------------------------------------
      subroutine mpoletoexp(nd,mpole,nterms,nlambs,numtets,nexptot,
     1                mexpupf,mexpdownf,rlsc)

c     This subroutine converts a multipole expansion into the
c     corresponding exponential moment function mexp for
c     both the +z direction and the -z direction
c
c     Note: this subroutine is the same as mpoletoexp
c     in mtxbothnew.f but just has a different data structure
c     for mpole
c
c     U(x,y,z) = \sum_{n=0}^{nterms} \sum_{m=-n,n} mpole(n,m)
c                Y_n^m (\cos(\theta)) e^{i m \phi}/r^{n+1}
c  
c              = (1/2\pi) \int_{0}^{\infty} e^{-\lambda z}
c                \int_{0}^{2\pi} e^{i \lambda (x \cos(\alpha) +
c                y \sin(\alpha))} mexpup(\lambda,\alpha)
c                d\alpha d\lambda
c 
c     for the +z direction and
c
c              = (1/2\pi) \int_{0}^{\infty} e^{\lambda z}
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
c     M_\lambda(m) = (i)**m \sum_{n=m}^{N} c(n,m) mpole(n,m)
c     lambda^n
c 
c     for m >=0 only, where c(n,m) = 1/sqrt((n+m)!(n-m)!)
c
c     For possible future reference, it should be noted that it
c     is NOT true that M_\lambda(-m) = dconjg(M_\lambda(m))
c
c     Inspection of the integral formula for Y_n^{-m} shows
c     that M_\lambda(-m) = dconjg(M_\lambda) * (-1)**m
c
c     INPUT arguments
c     nd          in: integer(8)
c                 number of multipole expansions
c
c     mpole       in: double complex (nd,0:nterms, -nterms:nterms)
c                 The multipole expansion 
c  
c     nterms:     in: integer(8)
c                 Order of the multipole expansion
c
c     nlambs      in: integer(8)
c                 number of discretization points in the \lambda
c                 integral
c
c     numtets     in: integer(8)(nlambs)
c                 number of fourier modes needed in expansion
c                 of \alpha variable for each \lambda variable
c
c     nexptot     in: integer(8)
c                 nexptot = \sum_{j} numtets(j)
c
c     rlsc        in: double precision(0:nterms, 0:nterms,nlambs)
c                 scaled discretization points in the \lambda
c                 integral
c
c     OUTPUT 
c     mexpupf     out: double complex (nd,nexptot)
c                 Fourier coefficients of the function
c                 mexpup(\lambda,\alpha) for successive
c                 discrete lambda values. They are ordered as
c                 follows
c
c                 mexpupf(1,...., numtets(1)) = fourier modes
c                             for \lambda_1
c
c                 mexpupf(numtets(1)+1,...., numters(2) = fourier
c                 modes for \lambda_2
c
c                 ETC
c
c     mexpdownf   out: double complex (nd,nexptot)
c                 Fourier coefficients of the function 
c                 mexpdown(\lambda,\alpha) for successive
c                 discrete \lambda values
c---------------------------------------------------------------

      implicit none
      integer(8) nd
      integer(8) nterms,nlambs,numtets(nlambs),nexptot
      double complex mpole(nd,0:nterms,-nterms:nterms)
      double complex mexpupf(nd,nexptot)
      double complex mexpdownf(nd,nexptot)
      double precision rlsc(0:nterms,0:nterms,nlambs)

c     Temp variables
      double complex, allocatable :: ztmp1(:),ztmp2(:)
      double complex zeyep
      double precision sgn
      integer(8) ntot,ncurrent,nl,mth,nm,idim

      allocate(ztmp1(nd),ztmp2(nd))

      ntot = 0
      do nl=1,nlambs
        sgn = -1.0d0
        zeyep = 1.0d0

        do mth = 0,numtets(nl)-1
          ncurrent = ntot + mth + 1

          do idim = 1,nd
            ztmp1(idim) = 0.0d0
            ztmp2(idim) = 0.0d0
          enddo

          sgn = -sgn
          do nm = mth,nterms,2
            do idim=1,nd
              ztmp1(idim) = ztmp1(idim) + 
     1            rlsc(nm,mth,nl)*mpole(idim,nm,mth)
            enddo
          enddo

          do nm=mth+1,nterms,2
            do idim=1,nd
              ztmp2(idim) = ztmp2(idim) + 
     1           rlsc(nm,mth,nl)*mpole(idim,nm,mth)
            enddo
          enddo

          do idim=1,nd
            mexpupf(idim,ncurrent) = (ztmp1(idim)+ztmp2(idim))*zeyep
            mexpdownf(idim,ncurrent) = 
     1         sgn*(ztmp1(idim)-ztmp2(idim))*zeyep
          enddo
          zeyep = zeyep*dcmplx(0.0d0,1.0d0)
        enddo
        ntot = ntot + numtets(nl)
      enddo

      return
      end

c -----------------------------------------------------------------
      subroutine exptolocal(nd,local,nterms,rlambs,whts,nlambs,numtets,
     1                     nthmax,nexptot,lexp1f,lexp2f,scale,rlsc)
c-----------------------------------------------------------------
c     INPUT arguments
c     nd               in: number of local expansions
c 
c     nterms           in: integer(8)
c                      Order of local expansion
c
c     rlambs           in: double precision(nlambs)
c                      discretization points in the \lambda integral
c
c     whts             in: double precision(nlambs)
c                      quadrature weights in \lambda integral
c
c     nlambs           in: integer(8)
c                      number of discretization points in \lambda
c                      integral
c
c     numtets          in: integer(8)(nlambs)
c                      number of fourier modes in expansion of
c                      \alpha variable for \lambda_j
c
c     nthmax           in: integer(8)
c                      max_j numtets(j)
c
c     nexptot          in: integer(8)
c                      sum_j numtets(j)
c                      
c
c     lexp1f(nd,nexptot)  double complex(nd,nexptot)
c                      Fourier coefficients of the function 
c                      lexp1 for discrete \lambda values
c                      in the +z direction
c                      They are ordered as follows:
c
c                      lexp1f(1,...,numtets(1)) = Fourier modes
c                      for \lambda_1
c                      lexp1f(numtets(1)+1,...,numtets(2) = Fourier
c                      modes for \lambda_2 etc
c
c
c     lexp2f(nd,nexptot)  double complex(nd,nexptot)
c                      Fourier coefficients of the function 
c                      lexp1 for discrete \lambda values
c                      in the -z direction
c                      They are ordered as follows:
c
c                      lexp1f(1,...,numtets(1)) = Fourier modes
c                      for \lambda_1
c                      lexp1f(numtets(1)+1,...,numtets(2) = Fourier
c                      modes for \lambda_2 etc
c
c     scale            in: double precision
c                      scaling parameter for local expansion
c
c     rlsc        in: double precision(nlambs, 0:nterms, 0:nterms)
c                 scaled discretization points in the \lambda
c                 integral
c
c     OUTPUT
c     local(nd,0:nterms,-nterms:nterms): output local expansion of order
c                                     nterms
        
      implicit none
      integer(8) nd
      integer(8) nterms,nlambs,numtets(nlambs),nexptot,nthmax
      integer(8) ncurrent,ntot,nl
      double complex local(nd,0:nterms,-nterms:nterms)
      double complex lexp1f(nd,nexptot),lexp2f(nd,nexptot)
      double complex zeye(0:nterms)
      double precision rlambs(nlambs), rlambpow(0:nterms) ,whts(nlambs)
      double precision rmul,rlsc(0:nterms,0:nterms,nlambs)
      double precision scale, rscale(0:nterms)
      double complex ima
    
c     Temporary variables
      integer(8) i, nm, mth, j, mmax,idim
      double precision dtmp

      data ima/(0.0d0,1.0d0)/


      zeye(0) = 1.0d0
      do i=1,nterms
        zeye(i) = zeye(i-1)*ima
      enddo

      rscale(0) = 1
      do nm=0,nterms
         if(nm.gt.0) rscale(nm) = rscale(nm-1)*scale
         do mth = -nterms,nterms
           do idim=1,nd
             local(idim,nm,mth) = 0.0d0
           enddo
         enddo
      enddo

      ntot = 1
      do nl=1,nlambs
c        Add contributions to local expansion
        do nm=0,nterms,2
          mmax = numtets(nl)-1
          if(mmax.gt.nm) mmax = nm
          do mth=0,mmax
            ncurrent = ntot+mth
            dtmp = rlsc(nm,mth,nl)*whts(nl)
            do idim=1,nd
              local(idim,nm,mth) = local(idim,nm,mth)+
     1          (lexp1f(idim,ncurrent)+lexp2f(idim,ncurrent))*dtmp
            enddo
          enddo
        enddo
        do nm=1,nterms,2
          mmax = numtets(nl) - 1
          if(mmax.gt.nm) mmax = nm
          do mth =0,mmax
            ncurrent = ntot+mth
            dtmp = -rlsc(nm,mth,nl)*whts(nl)
            do idim=1,nd
              local(idim,nm,mth) = local(idim,nm,mth)+
     1          (lexp1f(idim,ncurrent)-lexp2f(idim,ncurrent))*dtmp
            enddo
          enddo
        enddo
        ntot = ntot + numtets(nl)
      enddo

      do nm=0,nterms
        do idim=1,nd
          local(idim,nm,0) = local(idim,nm,0)*zeye(0)
        enddo
        do mth = 1,nm
          do idim=1,nd
            local(idim,nm,mth) = local(idim,nm,mth)*zeye(mth)
            local(idim,nm,-mth) = dconjg(local(idim,nm,mth))
          enddo
        enddo
      enddo

      return
      end
c------------------------------------------------
      subroutine phystof(nd,mexpf,nlambs,numfour,numphys,
     1                      mexpphys,fexpback)
      implicit double precision (a-h,o-z)
      integer(8) nd
      double complex mexpf(nd,*)
      double complex mexpphys(nd,*),ima
      double complex fexpback(*)
      double precision hh
      double precision, allocatable :: alphas(:)
      integer(8)  nlambs,numfour(nlambs),numphys(nlambs),nthmax
      data ima/(0.0d0,1.0d0)/
c
c     this subroutine converts the discretized exponential moment function
c     into its fourier expansion.
c
c     on input:
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
c     nlambs:        number of discretization pts. in lambda integral
c     numfour(j):   number of fourier modes in the expansion
c                      of the function m(\lambda_j,\alpha)
c     nthmax =      max_j numfour(j)
c
c     on output:
c
c     mexpf(nd,*):     fourier coefficients of the function 
c                   mexp(lambda,alpha) for discrete lambda values. 
c                   they are ordered as follows:
c
c               mexpf(1,...,numfour(1)) = fourier modes for lambda_1
c               mexpf(numfour(1)+1,...,numfour(2)) = fourier modes
c                                              for lambda_2
c               etc.
c
c------------------------------------------------------------
      done=1.0d0

      allocate(alphas(0:1000))
c
c
      pi=datan(done)*4
      nftot = 0
      nptot  = 0
      next  = 1
      do i=1,nlambs
        nalpha = numphys(i)
        hh = 1.0d0/nalpha
        halpha=2*pi*hh
        do j=1,nalpha
          alphas(j)=(j-1)*halpha
        enddo

        do idim=1,nd
          mexpf(idim,nftot+1) = 0.0d0
        enddo

        do ival=1,nalpha
          do idim=1,nd
            mexpf(idim,nftot+1) = mexpf(idim,nftot+1) + 
     1          mexpphys(idim,nptot+ival)*hh 
          enddo
        enddo

        do mm = 2,numfour(i)
          do idim=1,nd
            mexpf(idim,nftot+mm) = 0.0d0 
          enddo
          do ival=1,nalpha
            do idim=1,nd
              mexpf(idim,nftot+mm) = mexpf(idim,nftot+mm)+
     1          fexpback(next)*mexpphys(idim,nptot+ival)*hh
            enddo
            next = next+1
          enddo
        enddo
        nftot = nftot+numfour(i)
        nptot = nptot+numphys(i)
      enddo

      return
      end
c
c------------------------------------------------

c
      subroutine ftophys(nd,mexpf,nlambs,rlams,numfour,numphys,
     1                      nthmax,mexpphys,fexpe,fexpo)
      implicit double precision (a-h,o-z)
      integer(8) nd,nlambs
      double complex mexpf(nd,*)
      double complex mexpphys(nd,*),ima,ctmp
      double complex fexpe(*)
      double complex fexpo(*)
      double precision     rlams(nlambs)
      double precision, allocatable :: alphas(:) 
      integer(8) numfour(nlambs),numphys(nlambs),nthmax
      data ima/(0.0d0,1.0d0)/
c
c     this subroutine evaluates the fourier expansion of the
c     exponential moment function m(\lambda,\alpha) at equispaced
c     nodes.
c
c     on input:
c
c     mexpf(nd,*):     fourier coefficients of the function 
c                   mexp(lambda,alpha) for discrete lambda values. 
c                   they are ordered as follows:
c
c               mexpf(1,...,numfour(1)) = fourier modes for lambda_1
c               mexpf(numfour(1)+1,...,numfour(2)) = fourier modes
c                                              for lambda_2
c               etc.
c
c     nlambs:        number of discretization pts. in lambda integral
c     rlams(nlambs): discretization points in lambda integral.
c     numfour(j):   number of fourier modes in the expansion
c                      of the function m(\lambda_j,\alpha)
c     nthmax =      max_j numfour(j)
c     fexpe =      precomputed array of exponentials needed for
c                  fourier series evaluation
c     fexpo =      precomputed array of exponentials needed for
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

      allocate(alphas(0:1000))
c
c
      pi=datan(done)*4
      nftot = 0
      nptot  = 0
      nexte = 1
      nexto = 1
      do i=1,nlambs
        do ival=1,numphys(i)
          do idim=1,nd
            mexpphys(idim,nptot+ival) = mexpf(idim,nftot+1)
          enddo
          do mm = 2,numfour(i),2
            do idim=1,nd
              rtmp = 2*imag(fexpe(nexte)*mexpf(idim,nftot+mm))
              mexpphys(idim,nptot+ival) = mexpphys(idim,nptot+ival) +
     1                dcmplx(0.0d0,rtmp)
            enddo
            nexte = nexte + 1
          enddo
          do mm = 3,numfour(i),2
            do idim=1,nd
              rtmp = 2*real(fexpo(nexto)*mexpf(idim,nftot+mm))
              mexpphys(idim,nptot+ival) = mexpphys(idim,nptot+ival) +
     1                rtmp
            enddo
            nexto = nexto + 1
          enddo
        enddo
        nftot = nftot+numfour(i)
        nptot = nptot+numphys(i)
      enddo

      return
      end
c
c
c--------------------------------------------------------------------
      subroutine processudexp(nd,ibox,ilev,nboxes,centers,ichild,
     1           rscale,bs,nterms,iaddr,rmlexp,rlams,whts,nlams,
     2           nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nuall,uall,
     3           nu1234,u1234,ndall,dall,nd5678,d5678,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpuall,mexpu5678,mexpdall,
     5           mexpd1234,xs,ys,zs,fexpback,rlsc,rscpow,
     6           pgboxwexp,cntlist4,list4,nlist4s,ilist4,mnlist4)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) idim,nd
      integer(8) ibox,ilev,nboxes,nterms,nlams,nthmax
      integer(8) nphysical(nlams),nfourier(nlams)
      integer(8) iaddr(2,nboxes)
      integer(8) ichild(8,nboxes)
      integer(8) nexptot,nexptotp,nmax
      integer(8) nuall,ndall,nu1234,nd5678
      integer(8) uall(*),dall(*),u1234(*),d5678(*)
      double precision rscale,bs
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:)  
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision rmlexp(*),centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nd,nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpuall(nd,nexptotp),mexpdall(nd,nexptotp)
      double complex mexpd1234(nd,nexptotp),mexpu5678(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)
      integer(8) cntlist4,list4(*),nlist4s(*),ilist4(*),mnlist4
      integer(8) nlist4
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)

c      temp variables
      integer(8) jbox,ctr,ii,jj,i,ix,iy,iz,j,kbox
      double precision rtmp,rtmp2
      double complex ztmp,zmul,ztmp2
     
      double precision ctmp(3)

      allocate(tloc(nd,0:nterms,-nterms:nterms))


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
     1          mexp(idim,j,jbox,2)*zmul
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
     1          mexp(idim,j,jbox,2)*zmul
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
     1          mexp(idim,j,jbox,1)*zmul
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
     1         mexp(idim,j,jbox,1)*zmul
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
        
        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 2
      jbox = ichild(2,ibox)
      if(jbox.gt.0) then
        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*xs(1,i)
            mexpdownphys(idim,i) = (mexpdall(idim,i) + 
     1          mexpd1234(idim,i))*xs(-1,i)
          enddo
        enddo
 
        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif
  
c      add contributions due to child 3
      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpdall(idim,i) + 
     1          mexpd1234(idim,i))*ys(-1,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
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
     1         mexpd1234(idim,i))*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 5
      jbox = ichild(5,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1.0d0/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1           mexpu5678(idim,i))*zs(1,i)
            mexpdownphys(idim,i) = mexpdall(idim,i)*rtmp
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 6
      jbox = ichild(6,ibox)

      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpuall(idim,i)+
     1         mexpu5678(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
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

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
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

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,1,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(iaddr(2,jbox)),nterms)

      endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processnsexp(nd,ibox,ilev,nboxes,centers,ichild,
     1           rscale,bs,nterms,iaddr,rmlexp,rlams,whts,nlams,
     2           nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nnall,nall,
     3           nn1256,n1256,nn12,n12,nn56,n56,
     4           nsall,sall,ns3478,s3478,ns34,s34,ns78,s78,mexpup,
     5           mexpdown,mexpupphys,mexpdownphys,
     6           mexpnall,mexpn3478,mexpn34,mexpn78,mexpsall,
     7           mexps1256,mexps12,mexps56,rdplus,
     8           xs,ys,zs,fexpback,rlsc,rscpow,
     9           pgboxwexp,cntlist4,list4,nlist4s,ilist4,mnlist4)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) nd
      integer(8) ibox,ilev,nboxes,nterms,nlams,nthmax
      integer(8) nphysical(nlams),nfourier(nlams)
      integer(8) iaddr(2,nboxes)
      integer(8) ichild(8,nboxes)
      integer(8) nexptot,nexptotp,nmax
      integer(8) nnall,nsall,nn1256,ns3478,nn12,nn56,ns34,ns78
      integer(8) nall(*),sall(*),n1256(*),s3478(*)
      integer(8) n12(*),n56(*),s34(*),s78(*)
      double precision rscale,bs
      double complex zk2
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:)
      double complex, allocatable :: tloc2(:,:,:)
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision rdplus(0:nterms,0:nterms,-nterms:nterms)
      double precision rmlexp(*),centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nd,nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpnall(nd,nexptotp),mexpsall(nd,nexptotp)
      double complex mexps1256(nd,nexptotp),mexpn3478(nd,nexptotp)
      double complex mexps12(nd,nexptotp),mexps56(nd,nexptotp)
      double complex mexpn34(nd,nexptotp),mexpn78(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)
      integer(8) cntlist4,list4(*),nlist4s(*),ilist4(*),mnlist4
      integer(8) nlist4
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)

c      temp variables
      integer(8) jbox,ctr,ii,jj,i,ix,iy,iz,j,idim,kbox
      double complex ztmp,zmul,ztmp2
      double precision rtmp,rtmp2
    
      double precision ctmp(3)
      allocate(tloc(nd,0:nterms,-nterms:nterms))
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
     1           mexp(idim,j,jbox,4)*zmul
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
     1          mexp(idim,j,jbox,4)*zmul
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
            mexps12(idim,j) = mexps12(idim,j)+mexp(idim,j,jbox,4)*zmul
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
            mexps56(idim,j) = mexps56(idim,j)+mexp(idim,j,jbox,4)*zmul
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
     1         mexp(idim,j,jbox,3)*zmul
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
     1         mexp(idim,j,jbox,3)*zmul
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
            mexpn34(idim,j) = mexpn34(idim,j)+mexp(idim,j,jbox,3)*zmul
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
            mexpn78(idim,j) = mexpn78(idim,j)+mexp(idim,j,jbox,3)*zmul
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
            mexpdownphys(idim,i) = mexpsall(idim,i)+mexps1256(idim,i)+ 
     1         mexps12(idim,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 2
      jbox = ichild(2,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpsall(idim,i) + 
     1         mexps1256(idim,i) + mexps12(idim,i))*ys(-1,i)      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)


      endif
  
c      add contributions due to child 3
      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+mexpn34(idim,i)+
     1          mexpn3478(idim,i))*zs(1,i)
            mexpdownphys(idim,i) = mexpsall(idim,i)*rtmp
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 4
      jbox = ichild(4,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+mexpn34(idim,i)+
     1         mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
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

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
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
     1          mexps1256(idim,i) + mexps56(idim,i))*ztmp2      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

c      add contributions due to child 7
 
      jbox = ichild(7,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+mexpn78(idim,i)+
     1          mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

c      add contributions due to child 8
      jbox = ichild(8,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpnall(idim,i)+mexpn78(idim,i)+
     1         mexpn3478(idim,i))*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,2,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)


        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processewexp(nd,ibox,ilev,nboxes,centers,ichild,
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
     9           xs,ys,zs,fexpback,rlsc,rscpow,
     6           pgboxwexp,cntlist4,list4,nlist4s,ilist4,mnlist4)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) nd
      integer(8) ibox,ilev,nboxes,nterms,nlams,nthmax
      integer(8) nphysical(nlams),nfourier(nlams)
      integer(8) iaddr(2,nboxes)
      integer(8) ichild(8,nboxes)
      integer(8) nexptot,nexptotp,nmax
      integer(8) neall,nwall,ne1357,nw2468,ne13,ne57,nw24,nw68
      integer(8) ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
      integer(8) eall(*),wall(*),e1357(*),w2468(*)
      integer(8) e13(*),e57(*),w24(*),w68(*)
      integer(8) e1(*),e3(*),e5(*),e7(*),w2(*),w4(*),w6(*),w8(*)
      double precision rscale,bs
      double complex zk2
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:),tloc2(:,:,:)
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision rdminus(0:nterms,0:nterms,-nterms:nterms)
      double precision rmlexp(*),centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpeall(nd,nexptotp),mexpwall(nd,nexptotp)
      double complex mexpw1357(nd,nexptotp),mexpe2468(nd,nexptotp)
      double complex mexpw13(nd,nexptotp),mexpw57(nd,nexptotp)
      double complex mexpe24(nd,nexptotp),mexpe68(nd,nexptotp)
      double complex mexpw1(nd,nexptotp),mexpw3(nd,nexptotp)
      double complex mexpw5(nd,nexptotp),mexpw7(nd,nexptotp)
      double complex mexpe2(nd,nexptotp),mexpe4(nd,nexptotp)
      double complex mexpe6(nd,nexptotp),mexpe8(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)
      integer(8) cntlist4,list4(*),nlist4s(*),ilist4(*),mnlist4
      integer(8) nlist4
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)

c      temp variables
      integer(8) jbox,ctr,ii,jj,i,ix,iy,iz,j,l,idim,kbox
      double complex ztmp,zmul,ztmp2
      double precision rtmp,rtmp2
     
      double precision ctmp(3)

      allocate(tloc(nd,0:nterms,-nterms:nterms))
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
     1         mexp(idim,j,jbox,6)*zmul
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
     1         mexp(idim,j,jbox,6)*zmul
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
            mexpw13(idim,j) = mexpw13(idim,j)+mexp(idim,j,jbox,6)*zmul
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
            mexpw57(idim,j) = mexpw57(idim,j)+mexp(idim,j,jbox,6)*zmul
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
            mexpw1(idim,j) = mexpw1(idim,j) + mexp(idim,j,jbox,6)*zmul
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
             mexpw3(idim,j) = mexpw3(idim,j) + mexp(idim,j,jbox,6)*zmul
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
             mexpw5(idim,j) = mexpw5(idim,j) + mexp(idim,j,jbox,6)*zmul
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
            mexpw7(idim,j) = mexpw7(idim,j) + mexp(idim,j,jbox,6)*zmul
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
     1         mexp(idim,j,jbox,5)*zmul
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
     1          mexp(idim,j,jbox,5)*zmul
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
            mexpe24(idim,j) = mexpe24(idim,j)+mexp(idim,j,jbox,5)*zmul
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
            mexpe68(idim,j) = mexpe68(idim,j)+mexp(idim,j,jbox,5)*zmul
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
            mexpe2(idim,j) = mexpe2(idim,j) + mexp(idim,j,jbox,5)*zmul
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
            mexpe4(idim,j) = mexpe4(idim,j) + mexp(idim,j,jbox,5)*zmul
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
            mexpe6(idim,j) = mexpe6(idim,j) + mexp(idim,j,jbox,5)*zmul
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
            mexpe8(idim,j) = mexpe8(idim,j) + mexp(idim,j,jbox,5)*zmul
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
     1         mexpw13(idim,i)+mexpw1(idim,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 2
      jbox = ichild(2,ibox)

      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1         mexpe24(idim,i)+mexpe2(idim,i))*zs(1,i)      
            mexpdownphys(idim,i) = mexpwall(idim,i)*rtmp
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif
  
c      add contributions due to child 3
      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = (mexpwall(idim,i)+mexpw1357(idim,i)+
     1         mexpw13(idim,i)+mexpw3(idim,i))*ys(-1,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
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
     1         mexpe24(idim,i)+mexpe4(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 5
      jbox = ichild(5,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*xs(-1,i)
            mexpdownphys(idim,i) = (mexpwall(idim,i)+mexpw1357(idim,i)+
     1             mexpw57(idim,i)+mexpw5(idim,i))*xs(1,i)
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)
      endif

c      add contributions due to child 6
      jbox = ichild(6,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(-1,i)*zs(1,i)
          ztmp2 = xs(1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1           mexpe68(idim,i)+mexpe6(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
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
     1          mexpw57(idim,i)+mexpw7(idim,i))*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

c      add contributions due to child 8
      jbox = ichild(8,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(-1,i)*ys(1,i)*zs(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = (mexpeall(idim,i)+mexpe2468(idim,i)+
     1         mexpe68(idim,i)+mexpe8(idim,i))*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        nlist4=nlist4s(jbox)
        do i=1,nlist4
          kbox=ilist4((jbox-1)*mnlist4+i)
          call l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1         list4(kbox),nexptotp,xs,ys,zs,
     2         centers(1,kbox),centers(1,jbox),bs,3,cntlist4)
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(iaddr(2,jbox)),nterms)

      endif

      return
      end
c
c
c--------------------------------------------------------------------
      subroutine processlist3udexp(nd,ibox,nboxes,centers,
     1           bs,nterms,nexptotp,mexp,nuall,uall,
     3           ndall,dall,mexpuall,mexpdall,
     5           xs,ys,zs)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) idim,nd
      integer(8) ibox,nboxes,nterms,nlams,nthmax
      integer(8) nexptot,nexptotp
      integer(8) nuall,ndall
      integer(8) uall(*),dall(*)
      double precision bs
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision centers(3,nboxes)
      double complex mexpuall(nd,nexptotp),mexpdall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)

c      temp variables
      integer(8) jbox,i,ix,iy,iz,j
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
C        print *,"ulist j: ",jbox
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
C        print *,"dlist j: ",jbox
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

      subroutine processlist3nsexp(nd,ibox,nboxes,centers,
     1           bs,nterms,nexptotp,mexp,nnall,nall,
     3           nsall,sall,mexpnall,mexpsall,
     5           xs,ys,zs)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) nd
      integer(8) ibox,nboxes,nterms,nlams,nthmax
      integer(8) nexptotp
      integer(8) nnall,nsall
      integer(8) nall(*),sall(*)
      double precision bs
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision centers(3,*)
      double complex mexpnall(nd,nexptotp),mexpsall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)

c      temp variables
      integer(8) jbox,i,ix,iy,iz,j,idim
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

      subroutine processlist3ewexp(nd,ibox,nboxes,centers,
     1           bs,nterms,nexptotp,mexp,neall,eall,nwall,wall,
     4           mexpeall,mexpwall,xs,ys,zs)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) nd
      integer(8) ibox,nboxes,nterms,nlams,nthmax
      integer(8) nexptotp
      integer(8) neall,nwall
      integer(8) eall(*),wall(*)
      double precision bs 
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision centers(3,*)
      double complex mexpeall(nd,nexptotp),mexpwall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)

c      temp variables
      integer(8) jbox,i,ix,iy,iz,j,l,idim
      double complex ztmp,zmul,ztmp2
      double precision rtmp
     
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

      subroutine lpw_ud_eval_p(nd,center,boxsize,ntarg,targ,nlam,rlams,
     1   whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot)
      implicit none
      integer(8) nd
      real *8 center(3),boxsize,targ(3,ntarg),rlams(nlam),pot(nd,ntarg)
      real *8 whts(nlam)
      integer(8) ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:)
      integer(8) itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr
      complex *16 rz
      real *8, allocatable :: rexp(:),rexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(rexp(nlam),rexpinv(nlam),cc(nphmax))


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize


        do i=1,nlam
          rr = exp(-z*rlams(i))
          rexp(i) = rr*whts(i) 
          rexpinv(i) = whts(i)/rr
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            cc(iphys) = exp(ima*rlams(il)*(x*cos(alpha) + y*sin(alpha)))
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz = 
     1            (mexpupphys(idim,ii)*rexp(il)*cc(iphys) + 
     2             mexpdownphys(idim,ii)*rexpinv(il)*
     3             conjg(cc(iphys)))*hh
              pot(idim,itarg) = pot(idim,itarg) + real(rz) 
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

      subroutine lpw_ns_eval_p(nd,center,boxsize,ntarg,targ,nlam,rlams,
     1   whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot)
      implicit none
      integer(8) nd
      real *8 center(3),boxsize,targ(3,ntarg),rlams(nlam),pot(nd,ntarg)
      real *8 whts(nlam)
      integer(8) ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:)
      integer(8) itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr
      complex *16 rz
      real *8, allocatable :: rexp(:),rexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(rexp(nlam),rexpinv(nlam),cc(nphmax))


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

        do i=1,nlam
          rr = exp(-y*rlams(i))
          rexp(i) = rr*whts(i) 
          rexpinv(i) = whts(i)/rr
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            cc(iphys) = exp(ima*rlams(il)*(z*cos(alpha) + x*sin(alpha)))
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz = (mexpupphys(idim,ii)*rexp(il)*cc(iphys) + 
     2             mexpdownphys(idim,ii)*rexpinv(il)*
     3             conjg(cc(iphys)))*hh
              pot(idim,itarg) = pot(idim,itarg) + real(rz)
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

      subroutine lpw_ew_eval_p(nd,center,boxsize,ntarg,targ,nlam,rlams,
     1   whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot)
      implicit none
      integer(8) nd
      real *8 center(3),boxsize,targ(3,ntarg),rlams(nlam),pot(nd,ntarg)
      real *8 whts(nlam)
      integer(8) ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:)
      integer(8) itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr
      complex *16 rz
      real *8, allocatable :: rexp(:),rexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(rexp(nlam),rexpinv(nlam),cc(nphmax))


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

        do i=1,nlam
          rr = exp(-x*rlams(i))
          rexp(i) = rr*whts(i) 
          rexpinv(i) = whts(i)/rr
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            cc(iphys) = exp(ima*rlams(il)*(-z*cos(alpha)+y*sin(alpha)))
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz = (mexpupphys(idim,ii)*rexp(il)*cc(iphys) + 
     2             mexpdownphys(idim,ii)*rexpinv(il)*
     3             conjg(cc(iphys)))*hh
              pot(idim,itarg) = pot(idim,itarg) + real(rz)
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

      subroutine lpw_ud_eval_g(nd,center,boxsize,ntarg,targ,nlam,rlams,
     1   whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot,grad)
      implicit none
      integer(8) nd
      real *8 center(3),boxsize,targ(3,ntarg),rlams(nlam),pot(nd,ntarg)
      real *8 grad(nd,3,ntarg)
      real *8 whts(nlam)
      integer(8) ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),crc(:),crs(:)
      integer(8) itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr,binv
      complex *16 rz,rz1,rz2
      real *8, allocatable :: rexp(:),rexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(rexp(nlam),rexpinv(nlam),cc(nphmax))
      allocate(crc(nphmax),crs(nphmax))

      binv = 1.0d0/boxsize


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

        do i=1,nlam
          rr = exp(-z*rlams(i))
          rexp(i) = rr*whts(i) 
          rexpinv(i) = whts(i)/rr
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            crc(iphys) = rlams(il)*cos(alpha)*ima
            crs(iphys) = rlams(il)*sin(alpha)*ima
            cc(iphys) = exp(crc(iphys)*x + crs(iphys)*y)
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz1 = mexpupphys(idim,ii)*rexp(il)*cc(iphys)*hh
              rz2 = mexpdownphys(idim,ii)*rexpinv(il)*conjg(cc(iphys))*
     1              hh
              rz = rz1 + rz2
              pot(idim,itarg) = pot(idim,itarg) + real(rz)
              grad(idim,1,itarg) = grad(idim,1,itarg) + 
     1            real((rz1-rz2)*crc(iphys))*binv
              grad(idim,2,itarg) = grad(idim,2,itarg) + 
     1            real((rz1-rz2)*crs(iphys))*binv
              grad(idim,3,itarg) = grad(idim,3,itarg) - 
     1            real(rz1-rz2)*binv*rlams(il)
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

      subroutine lpw_ns_eval_g(nd,center,boxsize,ntarg,targ,nlam,rlams,
     1   whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot,grad)
      implicit none
      integer(8) nd
      real *8 center(3),boxsize,targ(3,ntarg),rlams(nlam),pot(nd,ntarg)
      real *8 grad(nd,3,ntarg)
      real *8 whts(nlam)
      integer(8) ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),crc(:),crs(:)
      integer(8) itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr,binv
      complex *16 rz,rz1,rz2
      real *8, allocatable :: rexp(:),rexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(rexp(nlam),rexpinv(nlam),cc(nphmax))
      allocate(crc(nphmax),crs(nphmax))

      binv = 1.0d0/boxsize


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

        do i=1,nlam
          rr = exp(-y*rlams(i))
          rexp(i) = rr*whts(i) 
          rexpinv(i) = whts(i)/rr
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            crc(iphys) = rlams(il)*cos(alpha)*ima
            crs(iphys) = rlams(il)*sin(alpha)*ima
            cc(iphys) = exp(crc(iphys)*z + crs(iphys)*x)
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz1 = mexpupphys(idim,ii)*rexp(il)*cc(iphys)*hh
              rz2 = mexpdownphys(idim,ii)*rexpinv(il)*conjg(cc(iphys))*
     1              hh
              rz = rz1 + rz2
              pot(idim,itarg) = pot(idim,itarg) + real(rz)
              grad(idim,1,itarg) = grad(idim,1,itarg) + 
     1            real((rz1-rz2)*crs(iphys))*binv
              grad(idim,2,itarg) = grad(idim,2,itarg) - 
     1            real(rz1-rz2)*binv*rlams(il)
              grad(idim,3,itarg) = grad(idim,3,itarg) + 
     1            real((rz1-rz2)*crc(iphys))*binv

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

      subroutine lpw_ew_eval_g(nd,center,boxsize,ntarg,targ,nlam,rlams,
     1   whts,nphys,nexptotp,nphmax,mexpupphys,mexpdownphys,pot,grad)
      implicit none
      integer(8) nd
      real *8 center(3),boxsize,targ(3,ntarg),rlams(nlam),pot(nd,ntarg)
      real *8 grad(nd,3,ntarg)
      real *8 whts(nlam)
      integer(8) ntarg,nlam,nphys(nlam),nexptotp,nphmax
      complex *16 mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      complex *16 ima
      complex *16, allocatable :: cc(:),crc(:),crs(:)
      integer(8) itarg,i,j,k,l,il,ii,iphys,istart,idim
      real *8 pi2inv,rexp1,alpha,pi2,x,y,z
      real *8 h,hh,rr,binv
      complex *16 rz,rz1,rz2
      real *8, allocatable :: rexp(:),rexpinv(:)
      data pi2inv/0.15915494309189535d0/
      data pi2/6.283185307179586d0/
      data ima/(0.0d0,1.0d0)/

      allocate(rexp(nlam),rexpinv(nlam),cc(nphmax))
      allocate(crc(nphmax),crs(nphmax))

      binv = 1.0d0/boxsize


      do itarg=1,ntarg
        x = (targ(1,itarg) - center(1))/boxsize
        y = (targ(2,itarg) - center(2))/boxsize
        z = (targ(3,itarg) - center(3))/boxsize

c
c         note: using vectorized notation to ensure simd instruction use
c
        rexpinv = 1.0d0/rexp

        do i=1,nlam
          rr = exp(-x*rlams(i))
          rexp(i) = rr*whts(i) 
          rexpinv(i) = whts(i)/rr
        enddo


        istart = 0
        do il=1,nlam
          h = pi2/nphys(il)
          hh = 1.0d0/nphys(il)

          do iphys = 1,nphys(il)
            alpha = (iphys-1)*h
            crc(iphys) = rlams(il)*cos(alpha)*ima
            crs(iphys) = rlams(il)*sin(alpha)*ima
            cc(iphys) = exp(-crc(iphys)*z + crs(iphys)*y)
          enddo
          do iphys = 1,nphys(il)
            ii = istart + iphys
            do idim=1,nd
              rz1 = mexpupphys(idim,ii)*rexp(il)*cc(iphys)*hh
              rz2 = mexpdownphys(idim,ii)*rexpinv(il)*conjg(cc(iphys))*
     1              hh
              rz = rz1 + rz2
              pot(idim,itarg) = pot(idim,itarg) + real(rz)
              grad(idim,1,itarg) = grad(idim,1,itarg) - 
     1            real(rz1-rz2)*binv*rlams(il)
              grad(idim,2,itarg) = grad(idim,2,itarg) + 
     1            real((rz1-rz2)*crs(iphys))*binv
              grad(idim,3,itarg) = grad(idim,3,itarg) - 
     1            real((rz1-rz2)*crc(iphys))*binv

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
c--------------------------------------------------------------------      
c
c
c--------------------------------------------------------------------
      subroutine processgboxudexp(nd,mexpugbox,mexpdgbox,jbox,
     1           nexptotp,mexpuall,mexpdall,
     2           xs,ys,zs)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) idim,nd
      integer(8) jbox,i
      integer(8) nexptotp
      double complex mexpugbox(nd,nexptotp)
      double complex mexpdgbox(nd,nexptotp)
      double complex mexpuall(nd,nexptotp)
      double complex mexpdall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rtmp
      double complex ztmp,ztmp2

c
cc       move all ghost box contributions to the child 1
c

c      add contributions due to child 1
      if(jbox.eq.1) then
        do i=1,nexptotp
          do idim=1,nd
            mexpuall(idim,i) = mexpuall(idim,i) + mexpugbox(idim,i)
            mexpdall(idim,i) = mexpdall(idim,i) + mexpdgbox(idim,i)
          enddo
        enddo
      endif
      
c      add contributions due to child 2
      if(jbox.eq.2) then
        do i=1,nexptotp
          do idim=1,nd
            mexpuall(idim,i) = mexpuall(idim,i) + 
     1                         mexpugbox(idim,i)*xs(-1,i)
            mexpdall(idim,i) = mexpdall(idim,i) + 
     1                         mexpdgbox(idim,i)*xs(1,i)
          enddo
        enddo
      endif
  
c      add contributions due to child 3
      if(jbox.eq.3) then
        do i=1,nexptotp
          do idim=1,nd
            mexpuall(idim,i) = mexpuall(idim,i) +
     1                         mexpugbox(idim,i)*ys(-1,i)
            mexpdall(idim,i) = mexpdall(idim,i) + 
     1                         mexpdgbox(idim,i)*ys(1,i)
          enddo
        enddo
      endif

c      add contributions due to child 4
      if(jbox.eq.4) then
        do i=1,nexptotp
          ztmp = ys(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)
          do idim=1,nd
            mexpuall(idim,i) = mexpuall(idim,i) + 
     1                         mexpugbox(idim,i)*ztmp2
            mexpdall(idim,i) = mexpdall(idim,i) + 
     1                         mexpdgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 5
      if(jbox.eq.5) then
        do i=1,nexptotp
          rtmp = 1.0d0/zs(1,i)
          do idim=1,nd
            mexpuall(idim,i) = mexpuall(idim,i) + 
     1                         mexpugbox(idim,i)*rtmp
            mexpdall(idim,i) = mexpdall(idim,i) + 
     1                         mexpdgbox(idim,i)*zs(1,i)
          enddo
        enddo
      endif

c      add contributions due to child 6
      if(jbox.eq.6) then
        do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpuall(idim,i) = mexpuall(idim,i) + 
     1                         mexpugbox(idim,i)*ztmp2
            mexpdall(idim,i) = mexpdall(idim,i) + 
     1                         mexpdgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 7
      if(jbox.eq.7) then
        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpuall(idim,i) = mexpuall(idim,i) + 
     1                         mexpugbox(idim,i)*ztmp2
            mexpdall(idim,i) = mexpdall(idim,i) + 
     1                         mexpdgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 8
      if(jbox.eq.8) then
        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)*xs(1,i)
          ztmp2 = xs(-1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpuall(idim,i) = mexpuall(idim,i) + 
     1                         mexpugbox(idim,i)*ztmp2
            mexpdall(idim,i) = mexpdall(idim,i) + 
     1                         mexpdgbox(idim,i)*ztmp
          enddo
        enddo
      endif

      return
      end
c--------------------------------------------------------------------      
c
c
c--------------------------------------------------------------------
      subroutine processgboxnsexp(nd,mexpngbox,mexpsgbox,jbox,
     1           nexptotp,mexpnall,mexpsall,
     2           xs,ys,zs)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) idim,nd
      integer(8) jbox,i
      integer(8) nexptotp
      double complex mexpngbox(nd,nexptotp)
      double complex mexpsgbox(nd,nexptotp)
      double complex mexpnall(nd,nexptotp)
      double complex mexpsall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rtmp
      double complex ztmp,ztmp2

c
cc       move all ghost box contributions to the child 1
c

c      add contributions due to child 1
      if(jbox.eq.1) then
        do i=1,nexptotp
          do idim=1,nd
            mexpnall(idim,i) = mexpnall(idim,i) + mexpngbox(idim,i)
            mexpsall(idim,i) = mexpsall(idim,i) + mexpsgbox(idim,i)
          enddo
        enddo
      endif
      
c      add contributions due to child 2
      if(jbox.eq.2) then
        do i=1,nexptotp
          do idim=1,nd
            mexpnall(idim,i) = mexpnall(idim,i) + 
     1                         mexpngbox(idim,i)*ys(-1,i)
            mexpsall(idim,i) = mexpsall(idim,i) + 
     1                         mexpsgbox(idim,i)*ys(1,i)
          enddo
        enddo
      endif
  
c      add contributions due to child 3
      if(jbox.eq.3) then
        do i=1,nexptotp
          rtmp = 1/zs(1,i)
          do idim=1,nd
            mexpnall(idim,i) = mexpnall(idim,i) +
     1                         mexpngbox(idim,i)*rtmp
            mexpsall(idim,i) = mexpsall(idim,i) + 
     1                         mexpsgbox(idim,i)*zs(1,i)
          enddo
        enddo
      endif

c      add contributions due to child 4
      if(jbox.eq.4) then
        do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpnall(idim,i) = mexpnall(idim,i) + 
     1                         mexpngbox(idim,i)*ztmp2
            mexpsall(idim,i) = mexpsall(idim,i) + 
     1                         mexpsgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 5
      if(jbox.eq.5) then
        do i=1,nexptotp
          do idim=1,nd
            mexpnall(idim,i) = mexpnall(idim,i) + 
     1                         mexpngbox(idim,i)*xs(-1,i)
            mexpsall(idim,i) = mexpsall(idim,i) + 
     1                         mexpsgbox(idim,i)*xs(1,i)
          enddo
        enddo
      endif

c      add contributions due to child 6
      if(jbox.eq.6) then
        do i=1,nexptotp
          ztmp = ys(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)
          do idim=1,nd
            mexpnall(idim,i) = mexpnall(idim,i) + 
     1                         mexpngbox(idim,i)*ztmp2
            mexpsall(idim,i) = mexpsall(idim,i) + 
     1                         mexpsgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 7
      if(jbox.eq.7) then
        do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpnall(idim,i) = mexpnall(idim,i) + 
     1                         mexpngbox(idim,i)*ztmp2
            mexpsall(idim,i) = mexpsall(idim,i) + 
     1                         mexpsgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 8
      if(jbox.eq.8) then
        do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpnall(idim,i) = mexpnall(idim,i) + 
     1                         mexpngbox(idim,i)*ztmp2
            mexpsall(idim,i) = mexpsall(idim,i) + 
     1                         mexpsgbox(idim,i)*ztmp
          enddo
        enddo
      endif

      return
      end
c--------------------------------------------------------------------      
c
c
c--------------------------------------------------------------------
      subroutine processgboxewexp(nd,mexpegbox,mexpwgbox,jbox,
     1           nexptotp,mexpeall,mexpwall,
     2           xs,ys,zs)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) idim,nd
      integer(8) jbox,i
      integer(8) nexptotp
      double complex mexpegbox(nd,nexptotp)
      double complex mexpwgbox(nd,nexptotp)
      double complex mexpeall(nd,nexptotp)
      double complex mexpwall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rtmp
      double complex ztmp,ztmp2

c
cc       move all ghost box contributions to the child 1
c

c      add contributions due to child 1
      if(jbox.eq.1) then
        do i=1,nexptotp
          do idim=1,nd
            mexpeall(idim,i) = mexpeall(idim,i) + mexpegbox(idim,i)
            mexpwall(idim,i) = mexpwall(idim,i) + mexpwgbox(idim,i)
          enddo
        enddo
      endif
      
c      add contributions due to child 2
      if(jbox.eq.2) then
        do i=1,nexptotp
          rtmp = 1/zs(1,i)
          do idim=1,nd
            mexpeall(idim,i) = mexpeall(idim,i) + 
     1                         mexpegbox(idim,i)*rtmp
            mexpwall(idim,i) = mexpwall(idim,i) + 
     1                         mexpwgbox(idim,i)*zs(1,i)
          enddo
        enddo
      endif
  
c      add contributions due to child 3
      if(jbox.eq.3) then
        do i=1,nexptotp
          do idim=1,nd
            mexpeall(idim,i) = mexpeall(idim,i) +
     1                         mexpegbox(idim,i)*ys(-1,i)
            mexpwall(idim,i) = mexpwall(idim,i) + 
     1                         mexpwgbox(idim,i)*ys(1,i)
          enddo
        enddo
      endif

c      add contributions due to child 4
      if(jbox.eq.4) then
        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpeall(idim,i) = mexpeall(idim,i) + 
     1                         mexpegbox(idim,i)*ztmp2
            mexpwall(idim,i) = mexpwall(idim,i) + 
     1                         mexpwgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 5
      if(jbox.eq.5) then
        do i=1,nexptotp
          do idim=1,nd
            mexpeall(idim,i) = mexpeall(idim,i) + 
     1                         mexpegbox(idim,i)*xs(1,i)
            mexpwall(idim,i) = mexpwall(idim,i) + 
     1                         mexpwgbox(idim,i)*xs(-1,i)
          enddo
        enddo
      endif

c      add contributions due to child 6
      if(jbox.eq.6) then
        do i=1,nexptotp
          ztmp = xs(-1,i)*zs(1,i)
          ztmp2 = xs(1,i)/zs(1,i)
          do idim=1,nd
            mexpeall(idim,i) = mexpeall(idim,i) + 
     1                         mexpegbox(idim,i)*ztmp2
            mexpwall(idim,i) = mexpwall(idim,i) + 
     1                         mexpwgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 7
      if(jbox.eq.7) then
        do i=1,nexptotp
          ztmp = xs(-1,i)*ys(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)
          do idim=1,nd
            mexpeall(idim,i) = mexpeall(idim,i) + 
     1                         mexpegbox(idim,i)*ztmp2
            mexpwall(idim,i) = mexpwall(idim,i) + 
     1                         mexpwgbox(idim,i)*ztmp
          enddo
        enddo
      endif

c      add contributions due to child 8
      if(jbox.eq.8) then
        do i=1,nexptotp
          ztmp = xs(-1,i)*ys(1,i)*zs(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpeall(idim,i) = mexpeall(idim,i) + 
     1                         mexpegbox(idim,i)*ztmp2
            mexpwall(idim,i) = mexpwall(idim,i) + 
     1                         mexpwgbox(idim,i)*ztmp
          enddo
        enddo
      endif

      return
      end
c--------------------------------------------------------------------   c
c
c--------------------------------------------------------------------
      subroutine l3dlist4pw(ilev,nd,nexptotp,nexptot,nterms,nmax,
     1           nlams,nlege,nthmax,nlevels,
     2           ifcharge,ifdipole,list4,itree,laddr,ipointer,
     3           nfourier,nphysical,rdminus,rdplus,rlsc,
     4           rscales,boxsize,zshift,sourcesort,chargesort,
     5           dipvecsort,centers,xshift,yshift,fexpe,fexpo,
     6           mexpf1,mexpf2,tmp,mptemp,wlege,rlams,rscpow,
     7           pgboxwexp,cntlist4)
c--------------------------------------------------------------------
c-------------------------------------------------------------------
      implicit none
ccc   input/output variables
      integer(8) ilev
      integer(8) nd
      integer(8) nexptotp,nexptot
      integer(8) nterms,nmax,nlams,nlege,nthmax
      integer(8) nlevels,cntlist4
      integer(8) ifcharge,ifdipole
      integer(8) list4(*),itree(*),laddr(2,0:nlevels)
      integer(8) ipointer(32)
      integer(8) nfourier(*)
      integer(8) nphysical(*)
      double precision rscales
      double precision boxsize
      double precision zshift(5,nexptotp)
      double precision mptemp((nmax+1)*(2*nmax+1)*2*nd)
      double precision wlege(*)
      double precision rlams(*)
      double precision rscpow(*)
      double precision sourcesort(3,*)
      double precision chargesort(nd,*)
      double precision dipvecsort(nd,3,*)
      double precision centers(3,*)
      double precision rdminus(0:nmax,0:nmax,-nmax:nmax)
      double precision rdplus(0:nmax,0:nmax,-nmax:nmax)
      double precision rlsc(0:nmax,0:nmax,nlams)
      double complex xshift(-5:5,nexptotp),yshift(-5:5,nexptotp)
      double complex fexpe(*),fexpo(*)
      double complex mexpf1(nd,nexptot),mexpf2(nd,nexptot)
      double complex tmp(nd,0:nmax,-nmax:nmax)
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)
ccc   scoped function variables
      integer(8) ibox,jbox,i,idim,nlist3,j
      integer(8) istart,iend,npts
      integer(8) jstart,jend,npts0
      integer(8) gboxfl(2,8)
      integer(8), allocatable :: gboxind(:)
      integer(8) itmp
      double precision time1,time2,omp_get_wtime
      double precision gboxsubcenters(3,8)
      double precision, allocatable ::  gboxsort(:,:)
      double precision, allocatable ::  gboxcgsort(:,:)
      double precision, allocatable ::  gboxdpsort(:,:,:)
      double complex, allocatable :: gboxmexp(:,:)
      double complex, allocatable :: gboxwexp(:,:,:,:)


c
c     count number of boxes are in list4 of this level
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      pgboxwexp=0d0
c     form mexp for all list4 type box at first ghost box center
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i)
C$OMP$PRIVATE(gboxind,gboxsort,gboxfl,gboxsubcenters)
C$OMP$PRIVATE(gboxwexp,gboxmexp,gboxcgsort,gboxdpsort)
C$OMP$PRIVATE(mexpf1,mexpf2,tmp,mptemp,itmp)
      do ibox=laddr(1,ilev),laddr(2,ilev)
        if(list4(ibox).gt.0) then
          istart=itree(ipointer(10)+ibox-1)
          iend=itree(ipointer(11)+ibox-1)
          npts = iend-istart+1
          if(npts.gt.0) then
            allocate(gboxind(npts))
            allocate(gboxsort(3,npts))
            allocate(gboxmexp(nd*(nterms+1)*(2*nterms+1),8))
            allocate(gboxwexp(nd,nexptotp,6,8))
            call subdividebox(sourcesort(1,istart),npts,
     1           centers(1,ibox),boxsize,
     2           gboxind,gboxfl,gboxsubcenters)
            itmp = 3
            call dreorderf(itmp,npts,sourcesort(1,istart),
     1           gboxsort,gboxind)
            if(ifcharge.eq.1) then
              allocate(gboxcgsort(nd,npts))
              call dreorderf(nd,npts,chargesort(1,istart),
     1             gboxcgsort,gboxind)
            endif
            if(ifdipole.eq.1) then
              allocate(gboxdpsort(nd,3,npts))
              call dreorderf(3*nd,npts,dipvecsort(1,1,istart),
     1             gboxdpsort,gboxind)
            endif
cccccccccccccc  bad code, note gboxmexp is an array not scalar
            gboxmexp=0
            gboxwexp=0
            do i=1,8
              if(gboxfl(1,i).gt.0) then
                jstart=gboxfl(1,i)
                jend=gboxfl(2,i)
                npts0=jend-jstart+1
                jbox=list4(ibox)
                if(npts0.gt.0) then
                  if(ifcharge.eq.1.and.ifdipole.eq.0) then
                    call l3dformmpc(nd,rscales,
     1                   gboxsort(1,jstart),
     2                   gboxcgsort(1,jstart),
     3                   npts0,gboxsubcenters(1,i),nterms,
     4                   gboxmexp(1,i),wlege,nlege)          
                  endif
                  if(ifcharge.eq.0.and.ifdipole.eq.1) then
                    call l3dformmpd(nd,rscales,
     1                   gboxsort(1,jstart),
     2                   gboxdpsort(1,1,jstart),
     3                   npts0,gboxsubcenters(1,i),nterms,
     4                   gboxmexp(1,i),wlege,nlege)          
                  endif
                  if(ifcharge.eq.1.and.ifdipole.eq.1) then
                    call l3dformmpcd(nd,rscales,
     1                   gboxsort(1,jstart),
     2                   gboxcgsort(1,jstart),
     3                   gboxdpsort(1,1,jstart),
     4                   npts0,gboxsubcenters(1,i),nterms,
     5                   gboxmexp(1,i),wlege,nlege)          
                  endif
ccc    convert to plane wave
                  call mpscale(nd,nterms,gboxmexp(1,i),
     1                 rscpow,tmp)
c
cc                process up down for current box
c
                  call mpoletoexp(nd,tmp,nterms,nlams,nfourier,
     1                 nexptot,mexpf1,mexpf2,rlsc)

                  call ftophys(nd,mexpf1,nlams,rlams,nfourier,nphysical,
     1                 nthmax,gboxwexp(1,1,1,i),fexpe,fexpo)

                  call ftophys(nd,mexpf2,nlams,rlams,nfourier,nphysical,
     1                 nthmax,gboxwexp(1,1,2,i),fexpe,fexpo)

                  call processgboxudexp(nd,gboxwexp(1,1,1,i),
     1                 gboxwexp(1,1,2,i),i,nexptotp,
     2                 pgboxwexp(1,1,jbox,1),
     3                 pgboxwexp(1,1,jbox,2),
     4                 xshift,yshift,zshift)
c
cc                process north-south for current box
c
                  call rotztoy(nd,nterms,tmp,mptemp,rdminus)
                  call mpoletoexp(nd,mptemp,nterms,nlams,nfourier,
     1                 nexptot,mexpf1,mexpf2,rlsc)

                  call ftophys(nd,mexpf1,nlams,rlams,nfourier,nphysical,
     1                 nthmax,gboxwexp(1,1,3,i),fexpe,fexpo)

                  call ftophys(nd,mexpf2,nlams,rlams,nfourier,nphysical,
     1                 nthmax,gboxwexp(1,1,4,i),fexpe,fexpo)

                  call processgboxnsexp(nd,gboxwexp(1,1,3,i),
     1                 gboxwexp(1,1,4,i),i,nexptotp,
     2                 pgboxwexp(1,1,jbox,3),
     3                 pgboxwexp(1,1,jbox,4),
     4                 xshift,yshift,zshift)

c
cc                process east-west for current box

                  call rotztox(nd,nterms,tmp,mptemp,rdplus)
                  call mpoletoexp(nd,mptemp,nterms,nlams,nfourier,
     1                 nexptot,mexpf1,mexpf2,rlsc)

                  call ftophys(nd,mexpf1,nlams,rlams,nfourier,nphysical,
     1                 nthmax,gboxwexp(1,1,5,i),fexpe,fexpo)

                  call ftophys(nd,mexpf2,nlams,rlams,nfourier,nphysical,
     1                 nthmax,gboxwexp(1,1,6,i),fexpe,fexpo)
                
                  call processgboxewexp(nd,gboxwexp(1,1,5,i),
     1                 gboxwexp(1,1,6,i),i,nexptotp,
     2                 pgboxwexp(1,1,jbox,5),
     3                 pgboxwexp(1,1,jbox,6),
     4                 xshift,yshift,zshift)
                endif
              endif
            enddo
            deallocate(gboxind,gboxsort)
            if(ifcharge.eq.1) then
              deallocate(gboxcgsort)
            endif
            if(ifdipole.eq.1) then
              deallocate(gboxdpsort)
            endif
            deallocate(gboxmexp)
            deallocate(gboxwexp)
          endif
        endif
      enddo
C$OMP END PARALLEL DO
      call cpu_time(time2)
C$    time2=omp_get_wtime()

      return
      end
c--------------------------------------------------------------------      
c--------------------------------------------------------------------   c
c
c--------------------------------------------------------------------
      subroutine l3dlist4shift(nd,mexpupphys,mexpdownphys,pgboxwexp,
     1           jbox,nexptotp,xs,ys,zs,
     2           censrc,centrg,boxsize,dirtype,cntlist4)
c--------------------------------------------------------------------
c-------------------------------------------------------------------
      implicit none
ccc   input/output variables
      integer(8) nd
      integer(8) nexptotp
      integer(8) jbox
      integer(8) cntlist4
      integer dirtype
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision boxsize
      double precision censrc(3),centrg(3)
ccc   scoped function variables
      integer(8) dir
      integer(8) i,ix,iy,iz,idim
      double complex zmul
      double precision rtmp
      double precision ctmp(3)


      call getlist4pwdir(dir,censrc,centrg,boxsize)

      ctmp(1) = censrc(1) - boxsize/2.0d0
      ctmp(2) = censrc(2) - boxsize/2.0d0
      ctmp(3) = censrc(3) - boxsize/2.0d0

      if(dir.eq.1.and.dirtype.eq.1) then
C        print *,"dir:",dir,"dirtype:",dirtype
        ix = 1.05d0*(ctmp(1)-centrg(1))/boxsize
        iy = 1.05d0*(ctmp(2)-centrg(2))/boxsize
        iz = 1.05d0*(ctmp(3)-centrg(3))/boxsize
        do i=1,nexptotp
          zmul = zs(iz,i)*xs(ix,i)*ys(iy,i)
          do idim=1,nd
            mexpdownphys(idim,i) = mexpdownphys(idim,i) + 
     1          pgboxwexp(idim,i,jbox,2)*zmul
          enddo
        enddo
      else if(dir.eq.2.and.dirtype.eq.1) then
C        print *,"dir:",dir,"dirtype:",dirtype
        ix = 1.05d0*(ctmp(1)-centrg(1))/boxsize
        iy = 1.05d0*(ctmp(2)-centrg(2))/boxsize
        iz = 1.05d0*(ctmp(3)-centrg(3))/boxsize
        do i=1,nexptotp
          zmul = zs(-iz,i)*xs(-ix,i)*ys(-iy,i)
          do idim=1,nd
            mexpupphys(idim,i) = mexpupphys(idim,i) + 
     1          pgboxwexp(idim,i,jbox,1)*zmul
          enddo
        enddo
      else if(dir.eq.3.and.dirtype.eq.2) then
C        print *,"dir:",dir,"dirtype:",dirtype
        ix = 1.05d0*(ctmp(1)-centrg(1))/boxsize
        iy = 1.05d0*(ctmp(2)-centrg(2))/boxsize
        iz = 1.05d0*(ctmp(3)-centrg(3))/boxsize
        do i=1,nexptotp
          zmul = zs(iy,i)*xs(iz,i)*ys(ix,i)
          do idim=1,nd
            mexpdownphys(idim,i) = mexpdownphys(idim,i) + 
     1          pgboxwexp(idim,i,jbox,4)*zmul
          enddo
        enddo
      else if(dir.eq.4.and.dirtype.eq.2) then
C        print *,"dir:",dir,"dirtype:",dirtype
        ix = 1.05d0*(ctmp(1)-centrg(1))/boxsize
        iy = 1.05d0*(ctmp(2)-centrg(2))/boxsize
        iz = 1.05d0*(ctmp(3)-centrg(3))/boxsize
        do i=1,nexptotp
          zmul = zs(-iy,i)*xs(-iz,i)*ys(-ix,i)
          do idim=1,nd
            mexpupphys(idim,i) = mexpupphys(idim,i) + 
     1          pgboxwexp(idim,i,jbox,3)*zmul
          enddo
        enddo
      else if(dir.eq.5.and.dirtype.eq.3) then
C        print *,"dir:",dir,"dirtype:",dirtype
        ix = 1.05d0*(ctmp(1)-centrg(1))/boxsize
        iy = 1.05d0*(ctmp(2)-centrg(2))/boxsize
        iz = 1.05d0*(ctmp(3)-centrg(3))/boxsize
        do i=1,nexptotp
          zmul = zs(ix,i)*xs(-iz,i)*ys(iy,i)
          do idim=1,nd
            mexpdownphys(idim,i) = mexpdownphys(idim,i) + 
     1          pgboxwexp(idim,i,jbox,6)*zmul
          enddo
        enddo
      else if(dir.eq.6.and.dirtype.eq.3) then
C        print *,"dir:",dir,"dirtype:",dirtype
        ix = 1.05d0*(ctmp(1)-centrg(1))/boxsize
        iy = 1.05d0*(ctmp(2)-centrg(2))/boxsize
        iz = 1.05d0*(ctmp(3)-centrg(3))/boxsize
        do i=1,nexptotp
          zmul = zs(-ix,i)*xs(iz,i)*ys(-iy,i)
          do idim=1,nd
            mexpupphys(idim,i) = mexpupphys(idim,i) + 
     1          pgboxwexp(idim,i,jbox,5)*zmul
          enddo
        enddo
      else
C        print *,"dir:",dir
      endif
       
      return
      end
c--------------------------------------------------------------------      
c
c
c--------------------------------------------------------------------
      subroutine processlist3udexplong(nd,ibox,nboxes,centers,
     1           rscale,nterms,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nuall,uall,
     3           ndall,dall,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpuall,mexpdall,
     5           xs,ys,zs,fexpback,rlsc,rscpow)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) idim,nd
      integer(8) ibox,nboxes,nterms,nlams,nthmax
      integer(8) nphysical(nlams),nfourier(nlams)
      integer(8) nexptot,nexptotp
      integer(8) nuall,ndall
      integer(8) uall(*),dall(*)
      double precision rscale
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:)  
      double complex mexp(nd,nexptotp,nboxes,6)
      double complex rmlexp(nd*(nterms+1)*(2*nterms+1),8)
      double precision centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nd,nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpuall(nd,nexptotp),mexpdall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)

c      temp variables
      integer(8) jbox,i,ix,iy,iz,j
      double precision rtmp
      double complex ztmp,zmul,ztmp2
     
      double precision ctmp(3)

      allocate(tloc(nd,0:nterms,-nterms:nterms))


      do i=1,nexptotp
        do idim=1,nd
          mexpuall(idim,i) = 0
          mexpdall(idim,i) = 0
        enddo
      enddo
      
   
      ctmp(1) = centers(1,ibox) - rscale/2.0d0
      ctmp(2) = centers(2,ibox) - rscale/2.0d0
      ctmp(3) = centers(3,ibox) - rscale/2.0d0
  
      
      do i=1,nuall
        jbox = uall(i)
C        print *,"ulist j: ",jbox
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
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
C        print *,"dlist j: ",jbox
        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

        do j=1,nexptotp
          zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
          do idim=1,nd
            mexpuall(idim,j) = mexpuall(idim,j) + 
     1                         mexp(idim,j,jbox,1)*zmul
          enddo
        enddo
      enddo
  
C      print *,mexpuall
C      print *,mexpdall


c
cc       move contributions to the children
c

      jbox=1
c      add contributions due to child 1
C      jbox = ichild(1,ibox)
      if(jbox.gt.0) then
        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)
            mexpdownphys(idim,i) = mexpdall(idim,i)
          enddo
        enddo

       call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

C        print *,tloc
        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

c
c         NOTE: fix rscpow to be 1/rscpow
c
        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(1,1),nterms)
      endif
      
c      add contributions due to child 2
C      jbox = ichild(2,ibox)
      if(jbox.gt.0) then
        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*xs(1,i)
            mexpdownphys(idim,i) = mexpdall(idim,i)*xs(-1,i)
          enddo
        enddo
 
        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(1,2),nterms)

      endif
  
c      add contributions due to child 3
C      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = mexpdall(idim,i)*ys(-1,i)
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(1,3),nterms)

      endif

c      add contributions due to child 4
C      jbox = ichild(4,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(1,4),nterms)

      endif

c      add contributions due to child 5
C      jbox = ichild(5,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1.0d0/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpuall(idim,i)*zs(1,i)
            mexpdownphys(idim,i) = mexpdall(idim,i)*rtmp
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(1,5),nterms)

      endif

c      add contributions due to child 6
C      jbox = ichild(6,ibox)

      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpuall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(1,6),nterms)


      endif

c      add contributions due to child 7
C      jbox = ichild(7,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpuall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(1,7),nterms)

      endif

c      add contributions due to child 8
C      jbox = ichild(8,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)*xs(1,i)
          ztmp2 = xs(-1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpuall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpdall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call mpscale(nd,nterms,tloc,rscpow,tloc)
        call mpadd(nd,tloc,rmlexp(1,8),nterms)

      endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processlist3nsexplong(nd,ibox,nboxes,centers,
     1           rscale,nterms,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nnall,nall,
     3           nsall,sall,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpnall,mexpsall,
     5           rdplus,xs,ys,zs,fexpback,rlsc,rscpow)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) nd
      integer(8) ibox,nboxes,nterms,nlams,nthmax
      integer(8) nphysical(nlams),nfourier(nlams)
      integer(8) nexptot,nexptotp
      integer(8) nnall,nsall
      integer(8) nall(*),sall(*)
      double precision rscale
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:)
      double complex, allocatable :: tloc2(:,:,:)
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision rdplus(0:nterms,0:nterms,-nterms:nterms)
      double complex rmlexp(nd*(nterms+1)*(2*nterms+1),8)
      double precision centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nd,nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpnall(nd,nexptotp),mexpsall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)

c      temp variables
      integer(8) jbox,i,ix,iy,iz,j,idim
      double complex ztmp,zmul,ztmp2
      double precision rtmp
    
      double precision ctmp(3)
      allocate(tloc(nd,0:nterms,-nterms:nterms))
      allocate(tloc2(nd,0:nterms,-nterms:nterms))


      do i=1,nexptotp
        do idim=1,nd
          mexpnall(idim,i) = 0
          mexpsall(idim,i) = 0
        enddo
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
           do idim=1,nd
             mexpsall(idim,j) = mexpsall(idim,j) + 
     1                          mexp(idim,j,jbox,4)*zmul
           enddo
        enddo

      enddo

      do i=1,nsall
        jbox = sall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
        do j=1,nexptotp
          zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
          do idim=1,nd
            mexpnall(idim,j) = mexpnall(idim,j) + 
     1                         mexp(idim,j,jbox,3)*zmul
          enddo
        enddo
      enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
C      jbox = ichild(1,ibox)
      jbox=1
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)
            mexpdownphys(idim,i) = mexpsall(idim,i)
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,1),nterms)

      endif

c      add contributions due to child 2
C      jbox = ichild(2,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = mexpsall(idim,i)*ys(-1,i)      
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,2),nterms)


      endif
  
c      add contributions due to child 3
C      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpnall(idim,i)*zs(1,i)
            mexpdownphys(idim,i) = mexpsall(idim,i)*rtmp
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,3),nterms)

      endif

c      add contributions due to child 4
C      jbox = ichild(4,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpnall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,4),nterms)

      endif

c      add contributions due to child 5
C      jbox = ichild(5,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*xs(1,i)
            mexpdownphys(idim,i) = mexpsall(idim,i)*xs(-1,i)      
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,5),nterms)

      endif

c      add contributions due to child 6
C      jbox = ichild(6,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,6),nterms)
      endif

c      add contributions due to child 7
C      jbox = ichild(7,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          ztmp2 = xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpnall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,7),nterms)
      endif

c      add contributions due to child 8
C      jbox = ichild(8,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)*xs(1,i)
          ztmp2 = ys(-1,i)*xs(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpnall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpsall(idim,i)*ztmp2      
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nd,nterms,tloc,tloc2,rdplus)


        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,8),nterms)
      endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processlist3ewexplong(nd,ibox,nboxes,centers,
     1           rscale,nterms,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,neall,eall,
     3           nwall,wall,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpeall,mexpwall,
     5           rdminus,xs,ys,zs,fexpback,rlsc,rscpow)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
      implicit none
      integer(8) nd
      integer(8) ibox,nboxes,nterms,nlams,nthmax
      integer(8) nphysical(nlams),nfourier(nlams)
      integer(8) nexptot,nexptotp
      integer(8) neall,nwall
      integer(8) eall(*),wall(*)
      double precision rscale
      double precision rlams(*),whts(*)
      double complex, allocatable :: tloc(:,:,:),tloc2(:,:,:)
      double complex mexp(nd,nexptotp,nboxes,6)
      double precision rdminus(0:nterms,0:nterms,-nterms:nterms)
      double complex rmlexp(nd*(nterms+1)*(2*nterms+1),8)
      double precision centers(3,*)
      double complex mexpup(nd,nexptot),mexpdown(nexptot)
      double complex mexpupphys(nd,nexptotp),mexpdownphys(nd,nexptotp)
      double complex mexpeall(nd,nexptotp),mexpwall(nd,nexptotp)
      double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
      double complex fexpback(*)

c      temp variables
      integer(8) jbox,i,ix,iy,iz,j,l,idim
      double complex ztmp,zmul,ztmp2
      double precision rtmp
     
      double precision ctmp(3)

      allocate(tloc(nd,0:nterms,-nterms:nterms))
      allocate(tloc2(nd,0:nterms,-nterms:nterms))


      do i=1,nexptotp
        do idim=1,nd
          mexpeall(idim,i) = 0
          mexpwall(idim,i) = 0
        enddo
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
          do idim=1,nd
            mexpwall(idim,j) = mexpwall(idim,j) + 
     1                         mexp(idim,j,jbox,6)*zmul
          enddo
        enddo
      enddo

      do i=1,nwall
        jbox = wall(i)

        ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
        iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
        iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         
        do j=1,nexptotp
          zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
          do idim=1,nd
            mexpeall(idim,j) = mexpeall(idim,j) + 
     1                         mexp(idim,j,jbox,5)*zmul
          enddo
        enddo
      enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
C      jbox = ichild(1,ibox)
      jbox=1
      if(jbox.gt.0) then
        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)
            mexpdownphys(idim,i) = mexpwall(idim,i)
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,1),nterms)

      endif

c      add contributions due to child 2
C      jbox = ichild(2,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          rtmp = 1/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpeall(idim,i)*zs(1,i)      
            mexpdownphys(idim,i) = mexpwall(idim,i)*rtmp
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,2),nterms)

      endif
  
c      add contributions due to child 3
C      jbox = ichild(3,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*ys(1,i)
            mexpdownphys(idim,i) = mexpwall(idim,i)*ys(-1,i)
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,3),nterms)

      endif

c      add contributions due to child 4
C      jbox = ichild(4,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          ztmp2 = ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpeall(idim,i)*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,4),nterms)

      endif

c      add contributions due to child 5
C      jbox = ichild(5,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*xs(-1,i)
            mexpdownphys(idim,i) = mexpwall(idim,i)*xs(1,i)
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,5),nterms)
      endif

c      add contributions due to child 6
C      jbox = ichild(6,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(-1,i)*zs(1,i)
          ztmp2 = xs(1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  =  mexpeall(idim,i)*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,6),nterms)
      endif

c      add contributions due to child 7
C      jbox = ichild(7,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(-1,i)*ys(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*ztmp
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,7),nterms)

      endif

c      add contributions due to child 8
C      jbox = ichild(8,ibox)
      if(jbox.gt.0) then

        do i=1,nexptotp
          ztmp = xs(-1,i)*ys(1,i)*zs(1,i)
          ztmp2 = xs(1,i)*ys(-1,i)/zs(1,i)
          do idim=1,nd
            mexpupphys(idim,i)  = mexpeall(idim,i)*ztmp      
            mexpdownphys(idim,i) = mexpwall(idim,i)*ztmp2
          enddo
        enddo

        call phystof(nd,mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
        call phystof(nd,mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal(nd,tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nd,nterms,tloc,tloc2,rdminus)

        call mpscale(nd,nterms,tloc2,rscpow,tloc2)
        call mpadd(nd,tloc2,rmlexp(1,8),nterms)

      endif

      return
      end
c--------------------------------------------------------------------      
