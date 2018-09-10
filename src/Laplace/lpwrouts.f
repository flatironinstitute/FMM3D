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
      do 100 i = 1,100
	     facts(i) = facts(i-1)*dsqrt(i+0.0d0)
100   continue
c
      do 1000 nl = 1,nlambs
c
c     compute powers of lambda_nl
c
         rlampow(0) = 1.0d0
         rmul = rlams(nl)
         do 200 j = 1,nterms
            rlampow(j) = rlampow(j-1)*rmul
200      continue
         do 600 j = 0,nterms
            do 400 k = 0,j
               rlsc(j,k,nl) = rlampow(j)/(facts(j-k)*facts(j+k))
400         continue
600      continue
1000  continue
      return
      end
c-------------------------------------------------------------
      subroutine mkexps(rlams,nlambs,numphys,nexptotp,xs,ys,zs)
      implicit double precision (a-h,o-z)
      double complex ima
      double complex xs(-5:5,nexptotp)
      double complex ys(-5:5,nexptotp)
      double precision zs(5,nexptotp)
      double precision     rlams(nlambs),u
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
      do 400 nl = 1,nlambs
         hu=2*pi/numphys(nl)
         do 200 mth = 1,numphys(nl)
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
200      continue
         ntot = ntot+numphys(nl)
400   continue
      return
      end
c***********************************************************************
      subroutine mkfexp(nlambs,numfour,numphys,fexpe,fexpo,fexpback)
      implicit double precision (a-h,o-z)
      double complex ima
      double complex fexpe(1)
      double complex fexpo(1)
      double complex fexpback(1)
      integer  nlambs,numphys(nlambs),numfour(nlambs)
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
      do 600 i=1,nlambs
	     nalpha = numphys(i)
         halpha=2*pi/nalpha
         do 400 j=1,nalpha
            alpha=(j-1)*halpha
	        do 200 mm = 2,numfour(i),2
               fexpe(nexte)  = cdexp(ima*(mm-1)*alpha)
	           nexte = nexte + 1
200         continue
	        do 300 mm = 3,numfour(i),2
               fexpo(nexto)  = cdexp(ima*(mm-1)*alpha)
	           nexto = nexto + 1
300         continue
400      continue
600   continue
      next = 1
      do 1600 i=1,nlambs
	     nalpha = numphys(i)
         halpha=2*pi/nalpha
	     do 1400 mm = 2,numfour(i)
            do 1200 j=1,nalpha
               alpha=(j-1)*halpha
               fexpback(next)  = cdexp(-ima*(mm-1)*alpha)
	           next = next + 1
1200        continue
1400     continue
1600  continue
      return
      end

c---------------------------------------------------
      subroutine mpoletoexp(mpole,nterms,nlambs,numtets,nexptot,
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
c     mpole       in: double complex (0:nterms, -nterms:nterms)
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
c     rlsc        in: double precision(0:nterms, 0:nterms,nlambs)
c                 scaled discretization points in the \lambda
c                 integral
c
c     OUTPUT 
c     mexpupf     out: double complex (nexptot)
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
c     mexpdownf   out: double complex (nexptot)
c                 Fourier coefficients of the function 
c                 mexpdown(\lambda,\alpha) for successive
c                 discrete \lambda values
c---------------------------------------------------------------

      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot
      double complex mpole(0:nterms,-nterms:nterms)
      double complex mexpupf(*)
      double complex mexpdownf(*)
      double complex zeyep,ztmp1,ztmp2
      double precision rlsc(0:nterms,0:nterms,nlambs)

c     Temp variables
      double precision sgn
      integer ntot,ncurrent,nl,mth,nm

      ntot = 0
      do nl=1,nlambs
         sgn = -1.0d0
         zeyep = 1.0d0
         do mth = 0,numtets(nl)-1
            ncurrent = ntot + mth + 1
            ztmp1 = 0.0d0
            ztmp2 = 0.0d0
            sgn = -sgn
            do nm = mth,nterms,2
               ztmp1 = ztmp1 + rlsc(nm,mth,nl)*mpole(nm,mth)
            enddo

            do nm=mth+1,nterms,2
               ztmp2 = ztmp2 + rlsc(nm,mth,nl)*mpole(nm,mth)
            enddo
            mexpupf(ncurrent) = (ztmp1+ztmp2)*zeyep
            mexpdownf(ncurrent) = sgn*(ztmp1-ztmp2)*zeyep
            zeyep = zeyep*dcmplx(0.0d0,1.0d0)
         enddo
         ntot = ntot + numtets(nl)
      enddo

      return
      end

c -----------------------------------------------------------------
      subroutine exptolocal(local,nterms,rlambs,whts,nlambs,numtets,
     1                     nthmax,nexptot,lexp1f,lexp2f,scale,rlsc)
c-----------------------------------------------------------------
c     INPUT arguments
c     nterms           in: integer
c                      Order of local expansion
c
c     rlambs           in: double precision(nlambs)
c                      discretization points in the \lambda integral
c
c     whts             in: double precision(nlambs)
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
c     lexp1f(nexptot)  double complex(nexptot)
c                      Fourier coefficients of the function 
c                      lexp1 for discrete \lambda values.
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
c     local(0:nterms,-nterms:nterms): output local expansion of order
c                                     nterms
        
      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot,nthmax
      integer ncurrent,ntot,nl
      double complex local(0:nterms,-nterms:nterms)
      double complex lexp1f(nexptot),lexp2f(nexptot)
      double complex zeye(0:nterms)
      double precision rlambs(nlambs), rlambpow(0:nterms) ,whts(nlambs)
      double precision rmul,rlsc(0:nterms,0:nterms,nlambs)
      double precision scale, rscale(0:nterms)
      double complex ima
    
c     Temporary variables
      integer i, nm, mth, j, mmax

      data ima/(0.0d0,1.0d0)/


      zeye(0) = 1.0d0
      do i=1,nterms
         zeye(i) = zeye(i-1)*ima
      enddo

      rscale(0) = 1
      do nm=0,nterms
         if(nm.gt.0) rscale(nm) = rscale(nm-1)*scale
         do mth = -nterms,nterms
            local(nm,mth) = 0.0d0
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
               local(nm,mth) = local(nm,mth)+rlsc(nm,mth,nl)*
     1          (lexp1f(ncurrent)+lexp2f(ncurrent))*whts(nl)
            enddo
         enddo
         do nm=1,nterms,2
            mmax = numtets(nl) - 1
            if(mmax.gt.nm) mmax = nm
            rmul = rlambpow(nm)
            do mth =0,mmax
               ncurrent = ntot+mth
               local(nm,mth) = local(nm,mth)-rlsc(nm,mth,nl)*
     1          (lexp1f(ncurrent)-lexp2f(ncurrent))*whts(nl)
            enddo
         enddo
         ntot = ntot + numtets(nl)
      enddo

      do nm=0,nterms
         do mth = 0,nm
            local(nm,mth) = local(nm,mth)*zeye(mth)
            if(mth.gt.0) local(nm,-mth) = dconjg(local(nm,mth))
         enddo
      enddo

      return
      end
c------------------------------------------------
      subroutine phystof(mexpf,nlambs,numfour,numphys,
     1                      mexpphys,fexpback)
      implicit double precision (a-h,o-z)
      double complex mexpf(1)
      double complex mexpphys(1),ima
      double complex fexpback(1)
      double precision     alphas(0:100)
      integer  nlambs,numfour(nlambs),numphys(nlambs),nthmax
      data ima/(0.0d0,1.0d0)/
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
c     nthmax =      max_j numfour(j)
c
c     on output:
c
c     mexpf(*):     fourier coefficients of the function 
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
        do 800 mm = 2,numfour(i)
           mexpf(nftot+mm) = 0.0d0
           do 600 ival=1,nalpha
              mexpf(nftot+mm) = mexpf(nftot+mm) +
     1          fexpback(next)*mexpphys(nptot+ival)
                next = next+1
600        continue
           mexpf(nftot+mm) = mexpf(nftot+mm)/nalpha
800     continue
        nftot = nftot+numfour(i)
        nptot = nptot+numphys(i)
 2000 continue
      return
      end
c
c------------------------------------------------

c
      subroutine ftophys(mexpf,nlambs,rlams,numfour,numphys,
     1                      nthmax,mexpphys,fexpe,fexpo)
      implicit double precision (a-h,o-z)
      double complex mexpf(1)
      double complex mexpphys(1),ima,ctmp
      double complex fexpe(1)
      double complex fexpo(1)
      double precision     rlams(nlambs)
      double precision     alphas(0:200)
      integer  nlambs,numfour(nlambs),numphys(nlambs),nthmax
      data ima/(0.0d0,1.0d0)/
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
      nexte = 1
      nexto = 1
      do 2000 i=1,nlambs
        do 1200 ival=1,numphys(i)
           mexpphys(nptot+ival) = mexpf(nftot+1)
           do 200 mm = 2,numfour(i),2
              rt1 = dimag(fexpe(nexte))*dreal(mexpf(nftot+mm))
              rt2 = dreal(fexpe(nexte))*dimag(mexpf(nftot+mm))
	      rtmp = 2*(rt1+rt2)
              nexte = nexte + 1
              mexpphys(nptot+ival) = mexpphys(nptot+ival) +
     1                dcmplx(0.0d0,rtmp)
200        continue
           do 400 mm = 3,numfour(i),2
              rt1 = dreal(fexpo(nexto))*dreal(mexpf(nftot+mm))
              rt2 = dimag(fexpo(nexto))*dimag(mexpf(nftot+mm))
	      rtmp = 2*(rt1-rt2)
              nexto = nexto + 1
              mexpphys(nptot+ival) = mexpphys(nptot+ival) +
     1                dcmplx(rtmp,0.0d0)
400        continue
 1200   continue
        nftot = nftot+numfour(i)
        nptot = nptot+numphys(i)
 2000 continue
      return
      end
c
c
c--------------------------------------------------------------------
      subroutine processudexp(ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nuall,uall,
     3           nu1234,u1234,ndall,dall,nd5678,d5678,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpuall,mexpu5678,mexpdall,
     5           mexpd1234,xs,ys,zs,fexpback,rlsc,rscpow)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer nuall,ndall,nu1234,nd5678
       integer uall(*),dall(*),u1234(*),d5678(*)
       double precision rscale
       double precision rlams(*),whts(*)
       double complex tloc(0:nterms,-nterms:nterms)
       double complex mexp(nexptotp,nboxes,6)
       double precision rmlexp(*),centers(3,*)
       double complex mexpup(nexptot),mexpdown(nexptot)
       double complex mexpupphys(nexptotp),mexpdownphys(nexptotp)
       double complex mexpuall(nexptotp),mexpdall(nexptotp)
       double complex mexpd1234(nexptotp),mexpu5678(nexptotp)
       double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       double precision zs(5,nexptotp)
       double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
       double complex fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j
       double complex ztmp,zmul
     
       double precision ctmp(3)


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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

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


        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo


        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processnsexp(ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nnall,nall,
     3           nn1256,n1256,nn12,n12,nn56,n56,
     4           nsall,sall,ns3478,s3478,ns34,s34,ns78,s78,mexpup,
     5           mexpdown,mexpupphys,mexpdownphys,
     6           mexpnall,mexpn3478,mexpn34,mexpn78,mexpsall,
     7           mexps1256,mexps12,mexps56,rdplus,
     8           xs,ys,zs,fexpback,rlsc,rscpow)
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
       double precision rscale
       double complex zk2
       double complex rlams(*),whts(*)
       double complex tloc(0:nterms,-nterms:nterms)
       double complex tloc2(0:nterms,-nterms:nterms)
       double complex mexp(nexptotp,nboxes,6)
       double precision rdplus(0:nterms,0:nterms,-nterms:nterms)
       double precision rmlexp(*),centers(3,*)
       double complex mexpup(nexptot),mexpdown(nexptot)
       double complex mexpupphys(nexptotp),mexpdownphys(nexptotp)
       double complex mexpnall(nexptotp),mexpsall(nexptotp)
       double complex mexps1256(nexptotp),mexpn3478(nexptotp)
       double complex mexps12(nexptotp),mexps56(nexptotp)
       double complex mexpn34(nexptotp),mexpn78(nexptotp)
       double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       double precision zs(5,nexptotp)
       double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
       double complex fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j
       double complex ztmp,zmul
     
       double precision ctmp(3)


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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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


        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processewexp(ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,neall,eall,
     3           ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,ne3,e3,ne5,e5,
     4           ne7,e7,nwall,wall,nw2468,w2468,nw24,w24,nw68,w68,
     5           nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6           mexpup,mexpdown,mexpupphys,mexpdownphys,
     7           mexpeall,mexpe2468,mexpe24,mexpe68,mexpe2,mexpe4,
     8           mexpe6,mexpe8,mexpwall,mexpw1357,mexpw13,mexpw57,
     9           mexpw1,mexpw3,mexpw5,mexpw7,rdminus,
     9           xs,ys,zs,fexpback,rlsc,rscpow)
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
       double precision rscale
       double complex zk2
       double complex rlams(*),whts(*)
       double complex tloc(0:nterms,-nterms:nterms)
       double complex tloc2(0:nterms,-nterms:nterms)
       double complex mexp(nexptotp,nboxes,6)
       double precision rdminus(0:nterms,0:nterms,-nterms:nterms)
       double precision rmlexp(*),centers(3,*)
       double complex mexpup(nexptot),mexpdown(nexptot)
       double complex mexpupphys(nexptotp),mexpdownphys(nexptotp)
       double complex mexpeall(nexptotp),mexpwall(nexptotp)
       double complex mexpw1357(nexptotp),mexpe2468(nexptotp)
       double complex mexpw13(nexptotp),mexpw57(nexptotp)
       double complex mexpe24(nexptotp),mexpe68(nexptotp)
       double complex mexpw1(nexptotp),mexpw3(nexptotp)
       double complex mexpw5(nexptotp),mexpw7(nexptotp)
       double complex mexpe2(nexptotp),mexpe4(nexptotp)
       double complex mexpe6(nexptotp),mexpe8(nexptotp)
       double complex xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       double precision zs(5,nexptotp)
       double precision rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
       double complex fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j,l
       double complex ztmp,zmul
     
       double precision ctmp(3)


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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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


        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

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

        call exptolocal(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      
