c
c        This file contains a set of subroutines for the handling 
c        of Legendre expansions. It contains 19 subroutines that are 
c        user-callable. Following is a brief description of these 
c        subroutines.
c
c   legeexps - constructs Legendre nodes, and  corresponding Gaussian
c        weights. Also constructs the matrix v converting the 
c         coefficients of a legendre expansion into its values at 
c         the n Gaussian nodes, and its inverse u, converting the
c         values of a function at n Gaussian nodes into the
c         coefficients of the corresponding Legendre series.
c
c   legepol - evaluates a single Legendre polynomial (together
c         with its derivative) at the user-provided point
c
c   legepols - evaluates a bunch of Legendre polynomials
c         at the user-provided point
c   legepls2 - an accelerated version of legepols, evaluating a 
c         Legendre polynomials at the user-provided point; maximum
c         order of the polynomials to be evaluated is 290
c
c   legepolders - evaluates a bunch of Legendre polynomials
c         at the user-provided point, and the derivatives of the 
c         said polynomials 
c   legeinmt - for the user-specified n, constructs the matrices of 
c        spectral indefinite integration differentiation on the n 
c        Gaussian nodes on the interval [-1,1]. 
c
c   legeinte - computes the indefinite integral of the legendre 
c        expansion polin getting the expansion polout
c
c   legediff -  differentiates the legendre expansion polin getting 
c        the expansion polout
c
c   legefder - computes the value and the derivative of a Legendre 
c        expansion at point X in interval [-1,1]; this subroutine 
c        is not designed to be very efficient, but it does not
c        use any exdternally supplied arrays
c
c   legefde2 - the same as legefder, except it is desigmed to be 
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed
c   
c   legeexev - computes the value of a Legendre expansion with 
c        at point X in interval [-1,1]; same as legefder, but does
c        not compute the derivative of the expansion
c
c   legeexe2 - the same as legeexev, except it is desigmed to be 
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed
c
c   lematrin - constructs the matrix interpolating functions from 
c        the n-point Gaussian grid on the interval [-1,1] to an 
c        arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c   
c   levecin - constructs the coefficients of the standard 
c        interpolation formula connecting the values of a 
c        function at n Gaussian nodes on the interval [a,b] with
c        its value at the point x \in R^1
c
c   legeodev - evaluates at the point x a Legendre expansion
c        having only odd-numbered elements; this is a fairly 
c        efficient code, using external arrays that are 
c        precomputed
c
c   legeevev - evaluates at the point x a Legendre expansion
c        having only even-numbered elements; this is a fairly 
c        efficient code, using external arrays that are 
c        precomputed
c
c   legepeven - evaluates even-numbered Legendre polynomials 
c        of the argument x; this is a fairly efficient code, 
c        using external arrays that are precomputed
c
c   legepodd - evaluates odd-numbered Legendre polynomials 
c        of the argument x; this is a fairly efficient code, 
c        using external arrays that are precomputed
c
C   legefdeq - computes the value and the derivative of a
c        Legendre Q-expansion with coefficients coefs
C     at point X in interval (-1,1); please note that this is
c     the evil twin of the subroutine legefder, evaluating the
c     proper (P-function) Legendre expansion; this subroutine 
c        is not designed to be very efficient, but it does not
c        use any exdternally supplied arrays
c
c   legeq - calculates the values and derivatives of a bunch 
c        of Legendre Q-functions at the user-specified point 
c        x on the interval (-1,1)
c
c   legeqs - calculates the value and the derivative of a single
c        Legendre Q-function at the user-specified point 
c        x on the interval (-1,1)
c
c   legecfde - computes the value and the derivative of a Legendre 
c        expansion with complex coefficients at point X in interval 
c        [-1,1]; this subroutine is not designed to be very efficient, 
c        but it does not use any exdternally supplied arrays. This is
c        a complex version of the subroutine legefder.
c
c   legecfd2 - the same as legecfde, except it is designed to be 
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed. This is a complex version of the 
c        subroutine legefde2.
c
c   legecva2 - the same as legecfd2, except it is does not evaluate
c        the derivative of the function
c
c
        subroutine legeexps(itype,n,x,u,v,whts)
        implicit double precision (a-h,o-z)
        dimension x(*),whts(*),u(n,n),v(n,n)
c
c         this subroutine constructs the gaussiaqn nodes 
c         on the interval [-1,1], and the weights for the 
c         corresponding order n quadrature. it also constructs
c         the matrix v converting the coefficients
c         of a legendre expansion into its values at the n
c         gaussian nodes, and its inverse u, converting the
c         values of a function at n gaussian nodes into the
c         coefficients of the corresponding legendre series.
c         no attempt has been made to make this code efficient, 
c         but its speed is normally sufficient, and it is 
c         mercifully short.
c
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of gaussian nodes and weights to be generated
c  
c                 output parameters:
c
c  x - the order n gaussian nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its 
c         legendre expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term legendre expansion into its values at
c         n legendre nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only 
c         if itype .ge. 1
c
c       . . . construct the nodes and the weights of the n-point gaussian 
c             quadrature
c
        ifwhts=0
        if(itype. gt. 0) ifwhts=1
        call legewhts(n,x,whts,ifwhts)
c
c       construct the matrix of values of the legendre polynomials
c       at these nodes        
c
        if(itype .ne. 2) return
        do 1400 i=1,n
c
        call legepols(x(i),n-1,u(1,i) )
 1400 continue
c
        do 1800 i=1,n
        do 1600 j=1,n
        v(i,j)=u(j,i)
 1600 continue
 1800 continue
c
c       now, v converts coefficients of a legendre expansion
c       into its values at the gaussian nodes. construct its 
c       inverse u, converting the values of a function at 
c       gaussian nodes into the coefficients of a legendre 
c       expansion of that function
c
        do 2800 i=1,n
        d=1
        d=d*(2*i-1)/2
        do 2600 j=1,n
        u(i,j)=v(j,i)*whts(j)*d
 2600 continue
 2800 continue
        return
        end
c
c
c
c
c
        subroutine legewhts_old(n,ts,whts,ifwhts)
        implicit double precision (a-h,o-z)
        dimension ts(*),whts(*)
c
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on 
c        the interval [-1,1]
c
c                input parameters:
c
c  n - the number of nodes in the quadrature
c
c                output parameters:
c
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n) 
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c
c         use newton to find all roots of the legendre polynomial
c
        ts(n/2+1)=0
        do 2000 i=1,n/2
c
        xk=ts(i)
        ifout=0
        deltold=1
        do 1400 k=1,10
        call legepol(xk,n,pol,der)
        delta=-pol/der
ccccc         call prin2('delta=*',delta,1)
        xk=xk+delta
        if(abs(delta) .lt. eps) ifout=ifout+1
c
cccc        call prin2('delta=*',delta,1)

        
        if(ifout .eq. 3) goto 1600
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c
c       now, use the explicit integral formulae 
c       to obtain the weights
c
        if(ifwhts .eq. 0) return
        a=-1
        b=1
        do 2200 i=1,n/2+1
        call prodend(a,ts,n,i,fm)
        call prodend(b,ts,n,i,fp)
        whts(i)=fp-fm
        whts(n-i+1)=whts(i)
 2200 continue
        return
        end
c
c
c
c
c
        subroutine legewhts(n,ts,whts,ifwhts)
        implicit double precision (a-h,o-z)
        dimension ts(*),whts(*)
c
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on 
c        the interval [-1,1]
c
c                input parameters:
c
c  n - the number of nodes in the quadrature
c
c                output parameters:
c
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c
        eps=1.0d-14
        done = 1
        pi=datan(done)*4
        h=pi/(2.0d0*n)


        do i=1,n
          t=(2*i-1)*h
          ts(n-i+1)=cos(t)
        enddo

c
c         use newton to find all roots of the legendre polynomial
c
        ts(n/2+1)=0
        do 2000 i=1,n/2
c
        xk=ts(i)
        ifout=0
        deltold=1
        do 1400 k=1,10
        call legepol_sum(xk,n,pol,der,sum)
        delta=-pol/der
        xk=xk+delta
        if(abs(delta) .lt. eps) ifout=ifout+1
c
        if(ifout .eq. 3) goto 1600
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c     
c        construct the weights via the orthogonality relation
c
        if(ifwhts .eq. 0) return
c
        do 2400 i=1,(n+1)/2
        call legepol_sum(ts(i),n,pol,der,sum)
        whts(i)=1/sum
        whts(n-i+1)=whts(i)
 2400 continue
c
        return
        end
c
c
c
c
c
        subroutine legepol_sum(x,n,pol,der,sum)
        implicit double precision (a-h,o-z)
c
        done=1
        sum=0 
c
        pkm1=1
        pk=x
        sum=sum+pkm1**2 /2
        sum=sum+pk**2 *(1+done/2)
c
        pk=1
        pkp1=x
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200

        sum=0 
c
        pol=1
        der=0
        sum=sum+pol**2 /2
        if(n .eq. 0) return
c
        pol=x
        der=1
        sum=sum+pol**2*(1+done/2)
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        sum=sum+pkp1**2*(k+1+done/2)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c
c
c
c
c
        subroutine legepol(x,n,pol,der)
        implicit double precision (a-h,o-z)
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c
        pol=x
        der=1
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c
c
c
c
c
        subroutine prodend(x,xs,n,i,f)
        implicit double precision (a-h,o-z)
        dimension xs(*)
c
c      evaluate the product
c
        f=1
        dlarge=1.0d20
        dsmall=f/dlarge
c
        large=0
        do 2000 j=1,n
        dd=dabs(f)
        if( dd .gt. dsmall) goto 1200
         f=f*10000
         large=large-1
 1200 continue
c
        if( dd .lt. dlarge) goto 1400
        f=f/10000
        large=large+1
 1400 continue
        if(j .eq. i) goto 2000
        f=f*(x-xs(j))/(xs(i)-xs(j))
 2000 continue
        d10000=10000
        f=f*d10000**large
        f=f**2*(x-xs(i))
        return
        end
c
c
c
c
c
        subroutine legepols(x,n,pols)
        implicit double precision (a-h,o-z)
        dimension pols(*)
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=x
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c
        return
        end

c
c
c
c
c
      SUBROUTINE legepolders(X,VALs,ders,N)
      IMPLICIT double precision (A-H,O-Z)
      double precision vals(*),ders(*)
C
C     This subroutine computes the values and the derivatives
c     of n+1 first Legendre polynomials at the point x 
C     in interval [-1,1].
c
c                input parameters:
c
C     X = evaluation point
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VALs = computed values of Legendre polynomials
C     ders = computed values of the derivatives
C
C


        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        vals(1)=1
        ders(1)=0
c
        vals(2)=x
        ders(2)=1
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c
        derj=derj/j

        vals(j+1)=pj
        ders(j+1)=derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
        subroutine legepls2(x,n,pols)
        implicit double precision (a-h,o-z)
        dimension pols(*),pjcoefs1(2000),pjcoefs2(300)
        save
        data ifcalled/0/
c
c        if need be - initialize the arrays pjcoefs1, pjcoefs2
c
        if(ifcalled .eq. 1) goto 1100
c
        done=1
        ninit=290
        do 1050 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1050 continue
c 
        ifcalled=1
 1100 continue

        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=x
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=x*pk*pjcoefs1(k+1)+pkm1*pjcoefs2(k+1)
        pols(k+2)=pkp1
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legeinmt(n,ainte,adiff,x,whts,endinter,
     1      itype,w)
        implicit double precision (a-h,o-z)
        dimension ainte(*),w(*),x(*),whts(*),adiff(*),endinter(*)
c
c
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes 
c        on the interval [-1,1]. Actually, this is omnly a 
c        memory management routine. All the actual work is done
c        by the subroutine legeinm0 (see)
c
c                           input parameters:
c
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION: 
c       itype=1 means that only the matrix ainte will 
c               be constructed
c       itype=2 means that only the matrix adiff will 
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c
c                           output paramaters:
c
c  ainte - the matrix of spectral indefinite integration on 
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on 
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c  endinter - the interpolation coefficients converting the 
c          values of a function at n Gaussian nodes into its
c          value at 1 (the right end of the interval)
c
c                           work arrays:
c
c  w - must be 3* n**2 + 2*n +50 *8 locations long
c
c        . . . allocate memory for the construction of the integrating
c              matrix
c
        ipolin=1
        lpolin=n+5
c
        ipolout=ipolin+lpolin
        lpolout=n+5
c
        iu=ipolout+lpolout
        lu=n**2+1
c
        iv=iu+lu
        lv=n**2+1
c
        iw=iv+lv
        lw=n**2+1
c
        ltot=iw+lw
c
c        construct the integrating matrix
c
        call legeinm0(n,ainte,adiff,w(ipolin),w(ipolout),
     1      x,whts,w(iu),w(iv),w(iw),itype,endinter)
c
        return
        end
c
c
c
c
c
        subroutine legeinm0(n,ainte,adiff,polin,polout,
     1      x,whts,u,v,w,itype,endinter)
        implicit double precision (a-h,o-z)
        dimension ainte(n,n),u(n,n),v(n,n),w(n,n),
     1      endinter(*),x(n),whts(n),polin(n),polout(n),
     2      adiff(n,n)
c
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes 
c        on the interval [-1,1]
c
c                           input parameters:
c
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION: 
c       itype=1 means that only the matrix ainte will 
c               be constructed
c       itype=2 means that only the matrix adiff will 
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c
c                           output paramaters:
c
c  ainte - the matrix of spectral indefinite integration on 
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on 
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c
c                           work arrays:
c
c  polin, polout - must be n+3 double precision locations each
c
c  u, v, w - must be n**2+1 double precision locations each
c
c        . . . construct the matrices of the forward and inverse 
c              Legendre transforms
c
        itype2=2
        call legeexps(itype2,n,x,u,v,whts)
c
cccc         call prin2('after legeexps, u=*',u,n*n)
c
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the indefinite integral of that function
c
        if(itype. eq. 2) goto 2000
c
        do 1600 i=1,n
c
        do 1200 j=1,n+2
        polin(j)=0
 1200 continue
c
        polin(i)=1
c
        call legeinte(polin,n,polout)
c
        do 1400 j=1,n
        ainte(j,i)=polout(j)
 1400 continue
c
 1600 continue
c
cccc         call prin2('ainte initially is*',ainte,n*n)
c
c        multiply the three, obtaining the integrating matrix
c
        call matmul(ainte,u,w,n)
        call matmul(v,w,ainte,n)
c
 2000 continue
c
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the derivative of that function
c
        if(itype. eq. 1) goto 3000
c
        do 2600 i=1,n
c
        do 2200 j=1,n+2
        polin(j)=0
 2200 continue
c
        polin(i)=1
c
        call legediff(polin,n,polout)
c
        do 2400 j=1,n
        adiff(j,i)=polout(j)
cccc        ainte(i,j)=polout(j)
 2400 continue
c
 2600 continue
c
cccc         call prin2('adiff initially is*',adiff,n*n)
c
c        multiply the three, obtaining the integrating matrix
c
        call matmul(adiff,u,w,n)
        call matmul(v,w,adiff,n)
c
 3000 continue
c
c        construct the vector of interpolation coefficients
c        converting the values of a polynomial at the Gaussian
c        nodes into its value at the right end of the interval
c
        do 3400 i=1,n
c
        d=0
        do 3200 j=1,n
        d=d+u(j,i)
 3200 continue
        endinter(i)=d
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine legeinte(polin,n,polout)
        implicit double precision (a-h,o-z)
        dimension polin(*),polout(*)
c
c       this subroutine computes the indefinite integral of the 
c       legendre expansion polin getting the expansion polout
c
c
c                       input parameters:
c
c  polin - the legendre expansion to be integrated
c  n - the order of the expansion polin 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c
c                       output parameters:
c
c  polout - the legendre expansion of the integral of the function 
c         represented by the expansion polin
c
        do 1200 i=1,n+2
        polout(i)=0
 1200 continue
c
        do 2000 k=2,n+1
        j=k-1
c
cccc        polout(k+1)=polin(k)/(2*j+1)+polout(k+1)
        polout(k+1)=polin(k)/(2*j+1)
        polout(k-1)=-polin(k)/(2*j+1)+polout(k-1)
c
 2000 continue
c
        polout(2)=polin(1)+polout(2)
c
        dd=0
        sss=-1
        do 2200 k=2,n+1
c
        dd=dd+polout(k)*sss
        sss=-sss
 2200 continue
c
ccc        call prin2('dd=*',dd,1)
        polout(1)=-dd
c
        return
        end
c
c
c
c
c
        subroutine legediff(polin,n,polout)
        implicit double precision (a-h,o-z)
        dimension polin(*),polout(*)
c
c       this subroutine differentiates the legendre 
c       expansion polin getting the expansion polout
c
c
c                       input parameters:
c
c  polin - the legendre expansion to be differentiated
c  n - the order of the expansion polin 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c
c                       output parameters:
c
c  polout - the legendre expansion of the derivative of the function 
c         represented by the expansion polin
c
        do 1200 k=1,n+1
        polout(k)=0
 1200 continue
c
        pk=polin(n+1)
        pkm1=polin(n)
        pkm2=0
        do 2000 k=n+1,2,-1
c
        j=k-1
c         
        polout(k-1)=pk*(2*j-1)
        if(k .ge. 3) pkm2=polin(k-2)+pk
c
        pk=pkm1
        pkm1=pkm2
c
 2000 continue
         return
         end
c
c
c
c
c
      SUBROUTINE legeFDER(X,VAL,der,PEXP,N)
      IMPLICIT double precision (A-H,O-Z)
      double precision PEXP(*)
C
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients PEXP
C     at point X in interval [-1,1].
c
c                input parameters:
c
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
C


        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c
        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c
      RETURN
      END


c
c
c
c
c
      SUBROUTINE legeFDE2(X,VAL,der,PEXP,N,
     1    pjcoefs1,pjcoefs2,ninit)
      IMPLICIT double precision (A-H,O-Z)
      double precision PEXP(*),pjcoefs1(*),pjcoefs2(*)
c
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients PEXP
C     at point X in interval [-1,1].
c
c                input parameters:
c
C  X - evaluation point
C  PEXP - expansion coefficients
C  N  - order of expansion 
c  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call 
c      on a previous call to this subroutine. Please note that this
c      is only an input parameter if the parameter ninit (see below) 
c      has been set to 0; otherwise, these are output parameters
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c       elements of each of the arrays pjcoefs1, pjcoefs2. On the first 
c       call to this subroutine, ninit should be set to the maximum 
c       order n for which this subroutine might have to be called; 
c       on subsequent calls, ninit should be set to 0. PLEASE NOTE 
c       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEXE2. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C  VAL - computed value
C  der - computed value of the derivative
C
C
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 1600 J = 2,N
c
cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2


        val=val+pexp(j+1)*pj
c
cccc        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
        derj=pjcoefs1(j)*(pjm1+x*derjm1)+pjcoefs2(j)*derjm2

ccc         call prin2('derj=*',derj,1)


cccc        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 1600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
      SUBROUTINE legeexev(X,VAL,PEXP,N)
      IMPLICIT double precision (A-H,O-Z)
      double precision PEXP(*)
C
C     This subroutine computes the value o a Legendre
c     expansion with coefficients PEXP at point X in interval [-1,1]
C
c                input parameters:
c
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C
        done=1
        pjm2=1
        pjm1=x
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c
        RETURN
        END
c
c
c
c
c
      SUBROUTINE legeexe2(X,VAL,PEXP,N,
     1      pjcoefs1,pjcoefs2,ninit)
      IMPLICIT double precision (A-H,O-Z)
      double precision PEXP(*),pjcoefs1(*),pjcoefs2(*)
c
C     This subroutine computes the value o a Legendre
c     expansion with coefficients PEXP at point X in interval [-1,1]
C
c                input parameters:
c
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C
        done=1
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2

cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c
        RETURN
        END
c
c
c
c
c
        subroutine lematrin(n,m,xs,amatrint,ts,w)
        implicit double precision (a-h,o-z)
        dimension amatrint(m,n),xs(*),w(*),ts(*)
c
c
c        This subroutine constructs the matrix interpolating
c        functions from the n-point Gaussian grid on the interval [-1,1]
c        to an arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c
c                 Input parameters:
c
c  n - the number of interpolation nodes
c  m - the number of nodes to which the functions will be interpolated
c  xs - the points at which the function is to be interpolated
c
c                  Output parameters:
c
c  amatrint - the m \times n matrix conerting the values of a function
c        at the n Legendre nodes into its values at m user-specified
c        (arbitrary) nodes
c  ts - the n Gaussian nodes on the interval [-1,1]
c
c                  Work arrays:
c
c  w - must be at least 2*n**2+n + 100 double precision locations long 
c

        icoefs=1
        lcoefs=n+2
c
        iu=icoefs+lcoefs
        lu=n**2+10
c
        iv=iu+lu
c
        ifinit=1
        do 2000 i=1,m
c
        call levecin(n,xs(i),ts,w(iu),w(iv),w(icoefs),ifinit)
c
        do 1400 j=1,n
        amatrint(i,j)=w(j)
 1400 continue
c
        ifinit=0
 2000 continue
c
        return
        end

c
c
c
c
c
        subroutine levecin(n,x,ts,u,v,coefs,ifinit)
        implicit double precision (a-h,o-z)
        dimension u(n,n),v(n,n),ts(*),coefs(*)
c
c        This subroutine constructs the coefficients of the 
c        standard interpolation formula connecting the values of a 
c        function at n Gaussian nodes on the interval [a,b] with
c        its value at the point x \in R^1
c
c                 Input parameters:
c
c  n - the number of interpolation nodes
c  x - the points at which the function is to be interpolated
c  ts - the n Gaussian nodes on the interval [-1,1]; please note that
c        it is an input parameter only if the parameter ifinit (see 
c        below) has been set to 1; otherwise, it is an output parameter
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its 
c        legendre expansion; please note that
c        it is an input parameter only if the parameter ifinit (see 
c        below) has been set to 1; otherwise, it is an output parameter
c  ifinit - an integer parameter telling the subroutine whether it should
c        initialize the Legendre expander; 
c     ifinit=1 will cause the subroutine to perform the initialization
c     ifinit=0 will cause the subroutine to  skip the initialization
c
c                  Output parameters:
c
c  coefs - the interpolation coefficients
c
c                 Work arrays: 
c
c  v - must be at least n*n double precision locations long
c
c       . . . construct the n Gausian nodes on the interval [-1,1]; 
c             also the corresponding Gaussian expansion-evaluation 
c             matrices
c
        itype=2
        if(ifinit .ne.0) call legeexps(itype,n,ts,u,v,coefs)
c
c       evaluate the n Legendre polynomials at the point where the
c       functions will have to be interpolated
c
        call legepols(x,n+1,v)
c
c       apply the interpolation matrix to the ector of values 
c       of polynomials from the right 
c
        call lematvec(u,v,coefs,n)
        return
        end
c
c
c
c
c
        subroutine lematvec(a,x,y,n)
        implicit double precision (a-h,o-z)
        dimension a(n,n),x(n),y(n)
c
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(j,i)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c
c
c
c
c
        subroutine matmul(a,b,c,n)
        implicit double precision (a-h,o-z)
        dimension a(n,n),b(n,n),c(n,n)
c
        do 2000 i=1,n
        do 1800 j=1,n
        d=0
        do 1600 k=1,n
        d=d+a(i,k)*b(k,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
c
c
c
c
        entry matmua(a,b,c,n)
ccc          call prin2('in matmua, a=*',a,n**2)
ccc          call prin2('in matmua, b=*',b,n**2)
        do 3000 i=1,n
        do 2800 j=1,n
        d=0
        do 2600 k=1,n
        d=d+a(i,k)*b(j,k)
 2600 continue
        c(i,j)=d
 2800 continue
 3000 continue
ccc          call prin2('exiting, c=*',c,n**2)
        return
        end

c
c
c
c
c
        subroutine legeodev(x,nn,coefs,val,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit double precision (a-h,o-z)
        dimension coepnm1(*),coepnp1(*),
     1            coexpnp1(*),coefs(*)
c
c
c       This subroutine evaluates at the point x a Legendre expansion
c       having only odd-numbered elements
c
c                  Input parameters:
c
c  x - point on the interval [-1,1] at which the Legendre expansion 
c       is to be evaluated
c  nn - order of the expansion to be evaluated
c  coefs - odd-numbered coefficients of the Legendre expansion
c       to be evaluated at the point x (nn/2+2 of them things)
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEPODD. IF these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 double precision elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  val - the value at the point x of the Legendre expansion with 
c       coefficients coefs (see above) 
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 double precision elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters
c 
c       
        if(ninit .eq. 0) goto 1400
        done=1
        n=0
        i=0
c
        do 1200 nnn=2,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pi=x
        pip1=x*(2.5d0*x22-1.5d0)
c
        val=coefs(1)*pi+coefs(2)*pip1

        do 2000 i=1,nn/2-2
c
        pip2 = coepnm1(i)*pi +
     1      (coepnp1(i)+coexpnp1(i)*x22)*pip1
c
        val=val+coefs(i+2)*pip2
c
        pi=pip1
        pip1=pip2


 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legeevev(x,nn,coefs,val,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit double precision (a-h,o-z)
        dimension coepnm1(*),coepnp1(*),
     1            coexpnp1(*),coefs(*)
c
c
c       This subroutine evaluates at the point x a Legendre expansion
c       having only even-numbered elements
c
c                  Input parameters:
c
c  x - point on the interval [-1,1] at which the Legendre expansion 
c       is to be evaluated
c  nn - order of the expansion to be evaluated
c  coefs - even-numbered coefficients of the Legendre expansion
c       to be evaluated at the point x (nn/2+2 of them things)
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEPEVEN. IF these aqrrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 double precision elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  val - the value at the point x of the Legendre expansion with 
c       coefficients coefs (see above) 
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 double precision elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters
c       
        if(ninit .eq. 0) goto 1400
c
        done=1
        n=-1
        i=0
        do 1200 nnn=1,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pi=1
        pip1=1.5d0*x22-0.5d0
c
        val=coefs(1)+coefs(2)*pip1
c
c       n is greater than 2. conduct recursion
c
        do 2000 i=1,nn/2-2
c
        pip2 = coepnm1(i)*pi +
     1      (coepnp1(i)+coexpnp1(i)*x22) *pip1
        val=val+coefs(i+2)*pip2
c
        pi=pip1
        pip1=pip2
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legepeven(x,nn,pols,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit double precision (a-h,o-z)
        dimension pols(*),coepnm1(*),coepnp1(*),
     1            coexpnp1(*)
c
c       This subroutine evaluates even-numbered Legendre polynomials 
c       of the argument x, up to order nn+1
c
c                  Input parameters:
c
c  x - the argument for which the Legendre polynomials are 
c       to be evaluated
c  nn - the maximum order for which the Legendre polynomials are
c       to be evaluated
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine ill initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 double precision elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c       PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEVEV. IF these aqrrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c 
c                  Output parameters:
c
c  pols - even-numbered Legendre polynomials of the input parameter x
c         (nn/2+2 of them things)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 double precision elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEVEV. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c 
c       
        if(ninit .eq. 0) goto 1400
c
        done=1
        n=-1
        i=0
        do 1200 nnn=1,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pols(1)=1
        pols(2)=1.5d0*x22-0.5d0
c
c       n is greater than 2. conduct recursion
c
        do 2000 i=1,nn/2
c
        pols(i+2) = coepnm1(i)*pols(i) +
     1      (coepnp1(i)+coexpnp1(i)*x22) *pols(i+1)
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legepodd(x,nn,pols,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit double precision (a-h,o-z)
        dimension pols(*),coepnm1(*),coepnp1(*),
     1            coexpnp1(*)
c
c       This subroutine evaluates odd-numbered Legendre polynomials 
c       of the argument x, up to order nn+1
c
c                  Input parameters:
c
c  x - the argument for which the Legendre polynomials are 
c       to be evaluated
c  nn - the maximum order for which the Legendre polynomials are
c       to be evaluated
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 double precision elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  pols - the odd-numbered Legendre polynomials of the input parameter x
c         (nn/2+2 of them things)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_1(x),
c       pols(2) = P_3(x), pols(3) = P_5 (x), etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 double precision elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS 
c       SUBROUTINE ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES 
c       SUSED BY THE UBROUTINE LEGEODEV. IF these arrays have been 
c       initialized by one of these two subroutines, they do not need 
c       to be initialized by the other one.
c       
        if(ninit .eq. 0) goto 1400
        done=1
        n=0
        i=0
c
        do 1200 nnn=2,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pols(1)=x
        pols(2)=x*(2.5d0*x22-1.5d0)
c
        do 2000 i=1,nn/2
c
        pols(i+2) = coepnm1(i)*pols(i) +
     1      (coepnp1(i)+coexpnp1(i)*x22)*pols(i+1)
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legefdeq(x,val,der,coefs,n)
        implicit double precision (a-h,o-z)
        dimension coefs(*)
C
C     This subroutine computes the value and the derivative
c     of a Legendre Q-expansion with coefficients coefs
C     at point X in interval (-1,1); please note that this is
c     the evil twin of the subroutine legefder, evaluating the
c     proper (P-function) Legendre expansion
c
c                input parameters:
c
C  X = evaluation point
C  coefs = expansion coefficients
C  N  = order of expansion 
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
c
        val=0
        der=0
c
        d= log( (1+x) /(1-x) ) /2
        pkm1=d
        pk=d*x-1
c
        pk=d
        pkp1=d*x-1

        derk=(1/(1+x)+1/(1-x)) /2
        derkp1=d + derk *x 
c
        val=coefs(1)*pk+coefs(2)*pkp1
        der=coefs(1)*derk+coefs(2)*derkp1
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
c
        if(n .eq. 0) return
c
        return
 1200 continue
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
c
        derkm1=derk
        derk=derkp1
c
        derkp1= ( (2*k+1)*pk+(2*k+1)*x*derk - k*derkm1 )/(k+1)
c
        val=val+coefs(k+2)*pkp1
        der=der+coefs(k+2)*derkp1
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legeqs(x,n,pols,ders)
        implicit double precision (a-h,o-z)
        dimension pols(*),ders(*)
c
c       This subroutine calculates the values and derivatives of 
c       a bunch of Legendre Q-functions at the user-specified point 
c       x on the interval (-1,1)
c
c                     Input parameters:
c
c  x - the point on the interval [-1,1] where the Q-functions and 
c       their derivatives are to be evaluated
c  n - the highest order for which the functions are to be evaluated
c  
c                     Output parameters:
c
c  pols - the values of the Q-functions (the evil twins of the 
c       Legeendre polynomials) at the point x (n+1 of them things)
c  ders - the derivatives of the Q-functions (the evil twins of the 
c       Legeendre polynomials) at the point x (n+1 of them things)
c  
c
        d= log( (1+x) /(1-x) ) /2
        pkm1=d
        pk=d*x-1
c
        pk=d
        pkp1=d*x-1

        derk=(1/(1+x)+1/(1-x)) /2
        derkp1=d + derk *x 
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=pk
        ders(1)=derk
        if(n .eq. 0) return
c
        pols(2)=pkp1
        ders(2)=derkp1
        return
 1200 continue
c
        pols(1)=pk
        pols(2)=pkp1
c
c       n is greater than 2. conduct recursion
c
        ders(1)=derk
        ders(2)=derkp1
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
c
        derkm1=derk
        derk=derkp1
c
        derkp1= ( (2*k+1)*pk+(2*k+1)*x*derk - k*derkm1 )/(k+1)
        ders(k+2)=derkp1
 2000 continue
c
        return
        end
c
c
c
c
c


        subroutine legeq(x,n,pol,der)
        implicit double precision (a-h,o-z)
c
c       This subroutine calculates the value and derivative of 
c       a Legendre Q-function at the user-specified point 
c       x on the interval (-1,1)
c
c
c                     Input parameters:
c
c  x - the point on the interval [-1,1] where the Q-functions and 
c       their derivatives are to be evaluated
c  n - the order for which the function is to be evaluated
c  
c                     Output parameters:
c
c  pol - the value of the n-th Q-function (the evil twin of the 
c       Legeendre polynomial) at the point x 
c  ders - the derivatives of the Q-function at the point x 
c  
c
        d= log( (1+x) /(1-x) ) /2
        pk=d
        pkp1=d*x-1
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pol=d

        der=(1/(1+x)+1/(1-x)) /2

        if(n .eq. 0) return
c
        pol=pkp1
        der=d + der *x 
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end

c
c
c
c
c
      SUBROUTINE legecFDE(X,VAL,der,PEXP,N)
      IMPLICIT double precision (A-H,O-Z)
      double complex PEXP(*),val,der
C
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with complex coefficients PEXP
C     at point X in interval [-1,1].
c
c                input parameters:
c
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
C
        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c
        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
      SUBROUTINE legecFD2(X,VAL,der,PEXP,N,
     1    pjcoefs1,pjcoefs2,ninit)
      IMPLICIT double precision (A-H,O-Z)
      double precision pjcoefs1(*),pjcoefs2(*)
      double complex PEXP(*),val,der
c
C     This subroutine computes the value and the derivative
c     of a Legendre expansion with complex coefficients PEXP
C     at point X in interval [-1,1].
c
c                input parameters:
c
C  X - evaluation point
C  PEXP - expansion coefficients
C  N  - order of expansion 
c  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call 
c      on a previous call to this subroutine. Please note that this
c      is only an input parameter if the parameter ninit (see below) 
c      has been set to 0; otherwise, these are output parameters
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c       elements of each of the arrays pjcoefs1, pjcoefs2. On the first 
c       call to this subroutine, ninit should be set to the maximum 
c       order n for which this subroutine might have to be called; 
c       on subsequent calls, ninit should be set to 0. PLEASE NOTE 
c       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEXE2. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C  VAL - computed value
C  der - computed value of the derivative
C
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 1600 J = 2,N
c
cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2


        val=val+pexp(j+1)*pj
c
cccc        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
        derj=pjcoefs1(j)*(pjm1+x*derjm1)+pjcoefs2(j)*derjm2

ccc         call prin2('derj=*',derj,1)


cccc        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 1600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
      SUBROUTINE legecva2(X,VAL,PEXP,N,
     1    pjcoefs1,pjcoefs2,ninit)
      IMPLICIT double precision (A-H,O-Z)
      double precision pjcoefs1(*),pjcoefs2(*)
      double complex PEXP(*),val
c
C     This subroutine computes the value of a Legendre expansion 
c     with complex coefficients PEXP at point X in interval [-1,1].
c
c                input parameters:
c
C  X - evaluation point
C  PEXP - expansion coefficients
C  N  - order of expansion 
c  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call 
c      on a previous call to this subroutine. Please note that this
c      is only an input parameter if the parameter ninit (see below) 
c      has been set to 0; otherwise, these are output parameters
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c       elements of each of the arrays pjcoefs1, pjcoefs2. On the first 
c       call to this subroutine, ninit should be set to the maximum 
c       order n for which this subroutine might have to be called; 
c       on subsequent calls, ninit should be set to 0. PLEASE NOTE 
c       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEXE2. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C  VAL - computed value
C
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
c
        DO 1600 J = 2,N
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2
        val=val+pexp(j+1)*pj
c
        pjm2=pjm1
        pjm1=pj
 1600   CONTINUE
c
      RETURN
      END


