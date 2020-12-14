c
c  This file contains a suite of evaluation codes for associated 
c  Legendre functions with various scalings, arguent types, etc.
c  Following is a brief description of the user-callable subroutines.
c  (FORTRAN 77 VERSION).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Simple routines for real argument.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       ylgndr - evaluate normalized Legendre functions
c
c       ylgndr2 - evaluate normalized Legendre functions 
c                 and their derivatives
c
c       ylgndr2s - evaluate normalized Legendre functions and their
c            derivatives with scaling. 
c            For m>0, the Ynm(x) values are scaled by 1/sqrt(1-x^2),
c            For m>0, the Ynm(x) derivatives with respect to x are 
c            scaled by sqrt(1-x^2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Simple routines for real argument (not scaled by sqrt(2*n+1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       ylgndru - evaluate normalized Legendre functions
c
c       ylgndru2 - evaluate normalized Legendre functions 
c                 and their derivatives
c
c       ylgndru2s - evaluate normalized Legendre functions and their
c            derivatives with scaling. 
c            For m>0, the Ynm(x) values are scaled by 1/sqrt(1-x^2),
c            For m>0, the Ynm(x) derivatives with respect to x are 
c            scaled by sqrt(1-x^2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Simple routines for real argument
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       ylgndr2sm - evaluate normalized Legendre functions and their
c            derivatives with scaling. 
c            For m>0, the Ynm(x) values are scaled by 1/sqrt(1-x^2)**m,
c            For m>0, the Ynm(x) derivatives with respect to x are 
c            scaled by 1/sqrt(1-x^2)**(m-2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Accelerated codes with precompution for real argument.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       ylgndrini - precomputation routine for recursion coefficients.
c                   Used by "fast" routines ylgndrf, ylgndr2f, ylgndr2sf.
c
c       ylgndrf - faster evaluation of normalized Legendre functions. 
c                 Requires prior call to ylgndrini.
c
c       ylgndr2f - faster evaluation of normalized Legendre functions 
c                  and their derivatives. 
c                  Requires prior call to ylgndrini.
c
c       ylgndr2sf - faster evaluation of  normalized Legendre functions 
c                   and their derivatives (with scaling as in ylgndr2s).
c                   Requires prior call to ylgndrini.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Accelerated codes with precompution for real argument.       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       ylgndrfex - faster evaluation of normalized Legendre functions. 
c                 Requires prior call to ylgndrini.
c                 In this routine, internal scaling by 10^(-id) is used
c                 to handle recurrence underflow near the poles.
c
c       ylgndr2fex - faster evaluation of normalized Legendre functions 
c                  and their derivatives. 
c                  Requires prior call to ylgndrini.
c                 In this routine, internal scaling by 10^(-id) is used
c                 to handle recurrence underflow near the poles.
c
c       ylgndr2sfex - faster evaluation of  normalized Legendre functions 
c                   and their derivatives (with scaling as in ylgndr2s).
c                   Requires prior call to ylgndrini.
c                 In this routine, internal scaling by 10^(-id) is used
c                 to handle recurrence underflow near the poles.
c
c   ylgndrfe, ylgndr2fe, ylgndr2sfe are the corresponding low-level
c   routines that do not post-multiply function values and derivatives 
c   by the scaling factors 10^(id). 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Simple codes for complex argument.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       zylgndr - evaluate normalized Legendre functions 
c                 of a complex argument
c
c       zylgndr2 - evaluate normalized Legendre functions 
c                  of a complex argument and their derivatives 
c
c       zylgndr2s - evaluate normalized Legendre functions 
c                  of a complex argument and their derivatives 
c                  (with scaling as in ylgndr2s)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Accelerated codes with precompution for complex argument.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       zylgndrf - faster evaluation of normalized Legendre functions 
c                  of a complex argument.
c                  Requires prior call to ylgndrini.
c
c       zylgndr2f - faster evaluation of normalized Legendre functions 
c                   of a complex argument and their derivatives.
c                   Requires prior call to ylgndrini.
c
c       zylgndr2sf - faster evaluation of normalized Legendre functions 
c                    of a complex argument and their derivatives 
c                   (with scaling as in ylgndr2s).
c                   Requires prior call to ylgndrini.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Miscellaneous, special purpose codes.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       zylgndrbr - evaluate normalized Legendre functions of a 
c            complex argument with modified branch cut.
c
c       zylgndrsc - evaluate normalized Legendre functions of a complex
c            argument with scaling to prevent overflow for large values
c            of z.
c
c       ylgndrpm, ylgndrpm_opt - Given Y_nm(x), return Y_nm(-x)
c
c       ylgndr2pm, ylgndrpm_opt - Given Y_nm(x) and Y'_nm(x), 
c                                 return Y_nm(-x) and Y'_nm(-x)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Truncated evaluation routines: computation carried out for
c       orders up to specified parameter rather than max degree.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       ylgndr2s_trunc - evaluate normalized Legendre functions and their
c            derivatives (with scaling as in ylgndr2s), 
c            BUT recursion carried out only to m = m2 rather than nmax.
c
c       ylgndrf_trunc - faster evaluation of normalized Legendre functions. 
c            Same as ylgndrf, BUT recursion carried out only to m = m2 
c            rather than nmax.
c            Requires prior call to ylgndrini.
c
c       ylgndr2f_trunc - faster evaluation of  normalized Legendre 
c            functions and their derivatives. Same as ylgndr2f, 
c            BUT recursion carried out only to m = m2 rather than nmax.
c
c       ylgndr2sf_trunc - faster evaluation of normalized Legendre 
c            functions and their derivatives (with scaling). 
c            Same as ylgndr2sf, BUT recursion carried out only to 
c            m = m2 rather than nmax.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Accelerated codes with precompution for real argument: 
c       precomputation is done only once for orders up to nmax, 
c       the routines are able to revert to simple routines for
c       higher orders
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       ylgndrfwini - precomputation routine for recursion coefficients.
c                   Used by "fast" routines ylgndrfw, ylgndr2fw, ylgndr2sfw.
c
c       ylgndrfw - faster evaluation of normalized Legendre functions. 
c                 Requires prior call to ylgndrfwini.
c
c       ylgndr2fw - faster evaluation of normalized Legendre functions 
c                  and their derivatives. 
c                  Requires prior call to ylgndrfwini.
c
c       ylgndr2sfw - faster evaluation of  normalized Legendre functions 
c                   and their derivatives (with scaling as in ylgndr2s).
c                   Requires prior call to ylgndrfwini.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       (not scaled by sqrt(2*n+1))
c
c       ylgndrufw - faster evaluation of normalized Legendre functions. 
c                 Requires prior call to ylgndrfwini.
c
c       ylgndru2fw - faster evaluation of normalized Legendre functions 
c                  and their derivatives. 
c                  Requires prior call to ylgndrfwini.
c
c       ylgndru2sfw - faster evaluation of  normalized Legendre functions 
c                   and their derivatives (with scaling as in ylgndr2s).
c                   Requires prior call to ylgndrfwini.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     CODES
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ylgndr(nmax, x, y)
      implicit none
c
c     Evaluate normalized Legendre functions defined as:
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                 must be non-negative
c     x                    -1 <= x <= 1
c     y(0:nmax,0:nmax)     resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x)
c     for 0 <= n <= nmax  and  0 <= m <= n.  Other elements of y
c     will contain undefined values.
c
cf2py intent(in) nmax
cf2py intent(in) x
cf2py intent(out) y

      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), u
c
c
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndr2(nmax, x, y, d)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x
cf2py intent(out) y,d
c
      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du
c
      u=-sqrt((1-x)*(1+x))
      du=x/sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x/u**2
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndr2s(nmax, x, y, d)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     For m>0, 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x
cf2py intent(out) y,d
c
      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
      do n=m+2, nmax
        y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+(1-x**2)*y(m,m))*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+(1-x**2)*y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
      subroutine ylgndru(nmax, x, y)
      implicit none
c
c     Evaluate normalized Legendre functions defined as:
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                 must be non-negative
c     x                    -1 <= x <= 1
c     y(0:nmax,0:nmax)     resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x)
c     for 0 <= n <= nmax  and  0 <= m <= n.  Other elements of y
c     will contain undefined values.
c
cf2py intent(in) nmax
cf2py intent(in) x
cf2py intent(out) y

      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), u
c
c
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c      do n=0, nmax
c	 do m=0, n
c	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
c         enddo
c      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndru2(nmax, x, y, d)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x
cf2py intent(out) y,d
c
      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du
c
      u=-sqrt((1-x)*(1+x))
      du=x/sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x/u**2
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c      do n=0, nmax
c	 do m=0, n
c	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
c	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
c         enddo
c      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndru2s(nmax, x, y, d)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     For m>0, 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x
cf2py intent(out) y,d
c
      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, u2
c
      u=-sqrt((1-x)*(1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
      do n=m+2, nmax
        y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+u2*y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c      do n=0, nmax
c	 do m=0, n
c	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
c	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
c         enddo
c      enddo
      return
      end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      subroutine ylgndr2sm(nmax, x, y, d)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     For m>0, 
c          the functions are scaled by 1/sqrt(1-x**2)**m
c          the derivatives are scaled by 1/sqrt(1-x**2)**(m-2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x
cf2py intent(out) y,d
c
      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, u2
c
      u=-sqrt((1-x))*sqrt((1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
      do n=m+2, nmax
        y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+u2*y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       faster version for real argument 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine ylgndrini(nmax, rat1, rat2)
      implicit none
c
c     Precompute the recurrence coefficients for the fast
c     evaluation of normalized Legendre functions and their derivatives
c    
c     Parameters:
c       nmax                      must be non-negative
c       rat1(0:nmax,0:nmax)       recurrence coefficient
c       rat2(0:nmax,0:nmax)       recurrence coefficient
c
cf2py intent(in) nmax
cf2py intent(out) rat1, rat2
c
      integer(8) nmax, m, n
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      rat1(0,0)=1
      rat2(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  rat1(m,m)=sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  rat2(m,m)=1
	 if (m.lt.nmax)  rat1(m+1,m)=sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  rat2(m+1,m)=1
	 do n=m+2, nmax
	    rat1(n,m)=(2*n-1)
            rat2(n,m)=sqrt((n+m-1.0d0)*(n-m-1.0d0))
	    rat1(n,m)=rat1(n,m)/sqrt(dble(n-m)*(n+m))
	    rat2(n,m)=rat2(n,m)/sqrt(dble(n-m)*(n+m))
         enddo
      enddo
c
c      do m=0, nmax
c	 do n=m, nmax
c            rat1(m,n)=rat1(n,m)
c            rat2(m,n)=rat2(n,m)
c         enddo
c      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndrf(nmax, x, y, rat1, rat2)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nmax, x, rat1, rat2
cf2py intent(out) y
c
      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), u
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2f(nmax, x, y, d, rat1, rat2)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x, rat1, rat2
cf2py intent(out) y, d
c
      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
c
      u=-sqrt((1-x)*(1+x))
      du=x/sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x/u**2
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2sf(nmax, x, y, d, rat1, rat2)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x, rat1, rat2
cf2py intent(out) y, d
c
      integer(8) nmax, m, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, u2
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
      do n=m+2, nmax
        y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndru2sf(nmax, x, y, d, rat1, rat2)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) =  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx =  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x, rat1, rat2
cf2py intent(out) y, d
c
      integer(8) nmax, n, m
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, u2
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
      do n=m+2, nmax
        y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
c      do n=0, nmax
c	 do m=0, n
c	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
c	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
c        enddo
c      enddo
      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       complex valued Legendre functions
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine zylgndr(nmax, z, y)
      implicit none
c
c     Evaluate normalized Legendre function for complex argument
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     z                     double complex
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(z)
c     for 0 <= n <= nmax  and  0 <= m <= n.  Other elements of y
c     will contain undefined values.
c
cf2py intent(in) nmax, z
cf2py intent(out) y
c
      integer(8) nmax
      double complex z, y(0:nmax,0:nmax)
c
      integer(8) m,n
      double complex u
c
      u=-sqrt(1-z*z)
      y(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*z*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine zylgndr2(nmax, z, y, d)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c       for complex argument
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     z                     double complex
c     y(0:nmax,0:nmax)      resulting function values, double complex
c     d(0:nmax,0:nmax)      resulting derivative values, double complex
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention is valid for the derivative.
c
cf2py intent(in) nmax, z
cf2py intent(out) y, d
c
      integer(8) nmax, m, n
      double complex z, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du
      u=-sqrt(1-z*z)
      du=-z/u
      y(0,0)=1
      d(0,0)=0
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*z/u**2
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  d(m+1,m)=(z*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*z*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(z*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine zylgndr2s(nmax, z, y, d)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(z) with m>0 
c          the functions are scaled by 1/sqrt(1-z**2)
c          the derivatives are scaled by sqrt(1-z**2)
c
c
c      Ynm(z) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(z)
c
c      d Ynm(z) / dz = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(z) / dz
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     z                     double complex
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(z) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, z
cf2py intent(out) y, d
c
      integer(8) nmax, m, n
      double complex z, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u
      double complex ztmp

      u=-sqrt(1-z*z)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*sqrt(2*m+1.0d0)
      if (m.lt.nmax)  d(m+1,m)=(z*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
      do n=m+2, nmax
        y(n,m)=((2*n-1)*z*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        d(n,m)=((2*n-1)*(z*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*z
c
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(z*d(m,m)+(1-z**2)*y(m,m))*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*z*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(z*d(n-1,m)+(1-z**2)*y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine zylgndrf(nmax, z, y, rat1, rat2)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(z) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(z)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     z                     double complex
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(z) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nmax, z, rat1, rat2
cf2py intent(out) y
c
      integer(8) nmax, m, n
      double complex z, y(0:nmax,0:nmax), u
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-z*z)
      y(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*z*y(n-1,m)-rat2(n,m)*y(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine zylgndr2f(nmax, z, y, d, rat1, rat2)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(z) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(z)
c
c      d Ynm(z) / dz = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(z) / dz
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     z                     double complex
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(z) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, z, rat1, rat2
cf2py intent(out) y, d
c
      integer(8) nmax, m, n
      double complex z, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-z*z)
      du=z/sqrt(1-z*z)
      y(0,0)=1
      d(0,0)=0
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=(d(m-1,m-1)*u+du)*rat1(m,m)
ccc	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*z/u**2
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  d(m+1,m)=(z*d(m,m)+y(m,m))*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*z*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(z*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine zylgndr2sf(nmax, z, y, d, rat1, rat2)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(z) with m>0 
c          the functions are scaled by 1/sqrt(1-z**2)
c          the derivatives are scaled by sqrt(1-z**2)
c
c
c      Ynm(z) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(z)
c
c      d Ynm(z) / dz = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(z) / dz
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     z                     double complex
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(z) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, z, rat1, rat2
cf2py intent(out) y, d
c
      integer(8) nmax, m, n
      double complex z, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt(1-z*z)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*rat1(m+1,m)
      if (m.lt.nmax)  d(m+1,m)=(z*d(m,m)+y(m,m))*rat1(m+1,m)
      do n=m+2, nmax
        y(n,m)=rat1(n,m)*z*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(z*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*z
c
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(z*d(m,m)+(1-z**2)*y(m,m))*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*z*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(z*d(n-1,m)+(1-z**2)*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       complex valued Legendre functions, modified branch cut
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine zylgndrbr(nmax, z, y)
      implicit none
c
c     Evaluate normalized Legendre function for complex argument
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     z                     double complex
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(z)
c     for 0 <= n <= nmax  and  0 <= m <= n.  Other elements of y
c     will contain undefined values.
c
c       Modified branch cut at (0,+i) 
c
cf2py intent(in) nmax, z
cf2py intent(out) y
c
      integer(8) nmax
      double complex z, y(0:nmax,0:nmax)
c
      integer(8) m,n
      double complex u
c
      u=-sqrt(1-z*z)
c
c     branch cut at (0,+i), select the lower branch 
c     of complex square root
c
      if( imag(1-z*z) .gt. 0 .and. real(1-z*z) .lt. 0) u=+sqrt(1-z*z)
ccc      call prin2('in zylgndrbr, u=*', -u, 2)
ccc      call prin2('in zylgndrbr, 1-z^2=*', 1-z*z, 2)
c
      y(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*z*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       complex valued Legendre functions (scaled)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      subroutine zylgndrsc(nmax, z,scale, ysc)
      implicit none
c
c     Evaluate scaled versions of zylgndr.  zylgndr is the complex version
c     of the normalized Legendre function.  Scaling removes the possible
c     conditioning errors from zylgndr evaluated at large arguments.
c
c     y_nm = scale**(-n) * ysc_nm
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                      must be non-negative
c     z (double complex)            double complex argument
c     scale (real*8)
c          looks like scale wants to be O(1/z) and scale<1.
c     ysc(0:nmax,0:nmax)        resulting function values
c
cf2py intent(in) nmax, z, scale
cf2py intent(out) ysc
c
      integer(8) nmax
      double complex z, ysc(0:nmax,0:nmax)
      real*8 scale
c
      integer(8) m,n
      double complex u,ztmp
c
      ztmp = 1-z*z
      u=-sqrt(ztmp)
      if(abs(imag(z)).le.1.0d-16.and.abs(real(z)).gt.1) then
        if(imag(u).lt.0) u = dconjg(u)
      endif
      ysc(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  then
            ysc(m,m)=ysc(m-1,m-1)*scale*u*sqrt((2*m-1.0d0)/(2*m))
c           call prinf('m=*',m,1)
c           call prin2('ysc(m,m)=*',ysc(m,m),2)
         endif
	 if (m.lt.nmax)  then
            ysc(m+1,m)=z*scale*ysc(m,m)*sqrt(2*m+1.0d0)
         endif
	 do n=m+2, nmax
	    ysc(n,m)=((2*n-1)*scale*z*ysc(n-1,m) - 
     1           sqrt((n+m-1.0d0)*(n-m-1.0d0))*scale**2*ysc(n-2,m))
     2           /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    ysc(n,m)=ysc(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       truncated recurrences for Legendre functions:
c    
c       Ynm(theta) for n = 0,nmax  but m = -m2,...,m2  with m2 < nmax.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
      subroutine ylgndr2s_trunc(nmax, m2, x, y, d)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, m2, x
cf2py intent(out) y, d
c
      integer(8) nmax, m, m2, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
      do n=m+2, nmax
        y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, m2
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+(1-x**2)*y(m,m))*sqrt(2*m+1.0d0)
	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+(1-x**2)*y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, min(n,m2)
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
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
      subroutine ylgndrf_trunc(nmax, m2, x, y, rat1, rat2)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., min(m2,n).
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nmax, m2, x, rat1, rat2
cf2py intent(out) y
c
      integer(8) nmax, m, m2, n
      double precision x, y(0:nmax,0:nmax), u
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      do m=0, m2
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, min(n,m2)
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2f_trunc(nmax, m2, x, y, d, rat1, rat2)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c    
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c    
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c    
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., min(m2,n).
c    
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c    
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, m2, x, rat1, rat2
cf2py intent(out) y, d
c
      integer(8) nmax, m, m2, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      du=x/sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
      do m=0, m2
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x/u**2
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, min(n,m2)
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2sf_trunc(nmax, m2, x, y, d, rat1, rat2)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., min(n,m2).
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, m2, x, rat1, rat2
cf2py intent(out) y, d
c
      integer(8) nmax, m, m2, n
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, u2
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
c
      u=-sqrt((1-x)*(1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
      do n=m+2, nmax
        y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, m2
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	 do m=0, min(n,m2)
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Symmetries for Legendre functions
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine ylgndrpm(nterms,y)
        implicit none
        integer(8) n,m,nterms
        double precision y(0:nterms,0:nterms)
c
c       Given Y_nm(x), return Y_nm(-x)
c
        do n=0,nterms
           do m=0,n
              if( mod(n+m,2) .eq. 1 ) y(n,m)=-y(n,m)
           enddo       
        enddo
c
        return
        end
c
c
c
c
c
        subroutine ylgndr2pm(nterms,y,d)
        implicit none
        integer(8) nterms,n,m
        double precision y(0:nterms,0:nterms)
        double precision d(0:nterms,0:nterms)
c
c       Given Y_nm(x), return Y_nm(-x)
c       Given Y'_nm(x), return Y'_nm(-x)
c
cf2py intent(in) nterms, y
cf2py intent(out) d
c
        do n=0,nterms
           do m=0,n
              if( mod(n+m,2) .eq. 1 ) y(n,m)=-y(n,m)
              if( mod(n+m,2) .eq. 0 ) d(n,m)=-d(n,m)
           enddo       
        enddo
c
        return
        end
c
c
c
c
c
        subroutine ylgndrpm_opt(nterms,y)
        implicit none
        integer(8) nterms,n,m
        double precision y(0:nterms,0:nterms)
c
c       Given Y_nm(x), return Y_nm(-x)
c
        do n=0,nterms,2
           do m=1,n,2
              y(n,m)=-y(n,m)
           enddo       
        enddo
c
        do n=1,nterms,2
           do m=0,n,2
              y(n,m)=-y(n,m)
           enddo       
        enddo
c
        return
        end
c
c
c
c
c
        subroutine ylgndr2pm_opt(nterms,y,d)
        implicit none
        integer(8) nterms,n,m
        double precision y(0:nterms,0:nterms)
        double precision d(0:nterms,0:nterms)
c
c       Given Y_nm(x), return Y_nm(-x)
c       Given Y'_nm(x), return Y'_nm(-x)
c
cf2py intent(in) nterms, y
cf2py intent(out) d
c
        do n=0,nterms,2
           do m=0,n,2
              d(n,m)=-d(n,m)
           enddo       
           do m=1,n,2
              y(n,m)=-y(n,m)
           enddo       
        enddo
c
        do n=1,nterms,2
           do m=0,n,2
              y(n,m)=-y(n,m)
           enddo       
           do m=1,n,2
              d(n,m)=-d(n,m)
           enddo       
        enddo
c
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       fast version of real valued Legendre functions, with storage
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine ylgndrfwini(nmax, w, lw, lused)
      implicit none
c
c     Precompute the recurrence coefficients for the fast
c     evaluation of normalized Legendre functions and their derivatives
c    
c     Parameters:
c       nmax             must be non-negative
c       w                  contains rat1 and rat2 arrays
c
      integer(8) nmax,irat1,irat2,lw,lused
      double precision w(*)
c
cf2py intent(in) nmax, lw
cf2py intent(out) w, lused
c
      irat1=1
      irat2=1+(nmax+1)**2
      lused=2*(nmax+1)**2
      if( lused .gt. lw ) return
      
      call ylgndrini(nmax, w(irat1), w(irat2))
      return
      end
c
c
c
c
c
      subroutine ylgndrfw(nterms, x, y, w, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)          resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nterms, x, w, nmax
cf2py intent(out) y
c
      integer(8) nterms,nmax,irat1,irat2
      double precision x, y(0:nterms,0:nterms), w(*)
c
      if( nterms .le. nmax ) then
        irat1=1
        irat2=1+(nmax+1)**2
        call ylgndrfw0(nterms, x, y, w(irat1), w(irat2), nmax)
      else
        call ylgndr(nterms, x, y)
      endif
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2fw(nterms, x, y, d, w, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, w, nmax
cf2py intent(out) y, d
c
      integer(8) nterms,nmax,irat1,irat2
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision w(*)    
c
      if( nterms .le. nmax ) then
        irat1=1
        irat2=1+(nmax+1)**2
        call ylgndr2fw0(nterms, x, y, d, w(irat1), w(irat2), nmax)
      else
        call ylgndr2(nterms, x, y, d)
      endif
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2sfw(nterms, x, y, d, w, nmax)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, w, nmax
cf2py intent(out) y, d
c
      integer(8) nterms,nmax,irat1,irat2
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision w(*)
c
      if( nterms .le. nmax ) then
        irat1=1
        irat2=1+(nmax+1)**2
        call ylgndr2sfw0(nterms, x, y, d, w(irat1), w(irat2), nmax)
      else
        call ylgndr2s(nterms, x, y, d)
      endif
c
      return
      end
c
c
c
c
c
      subroutine ylgndrfw0(nterms, x, y, rat1, rat2, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y
c
      integer(8) nterms,nmax,m,n
      double precision x, y(0:nterms,0:nterms), u
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      do m=0, nterms
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.lt.nterms)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 do n=m+2, nterms
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
         enddo
      enddo
c
      do n=0, nterms
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2fw0(nterms, x, y, d, rat1, rat2, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y, d
c
      integer(8) nterms, nmax, m, n
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision u, du
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      du=x/sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
      do m=0, nterms
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x/u**2
	 if (m.lt.nterms)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nterms)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
	 do n=m+2, nterms
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      do n=0, nterms
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2sfw0(nterms, x, y, d, rat1, rat2, nmax)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y, d
c
      integer(8) nterms, nmax, m, n
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision u, u2
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nterms)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      if (m.lt.nterms)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
      do n=m+2, nterms
        y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nterms
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nterms)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nterms)  
ccc     $      d(m+1,m)=(x*d(m,m)+(1-x**2)*y(m,m))*rat1(m+1,m)
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
	 do n=m+2, nterms
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
ccc	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+(1-x**2)*y(n-1,m))-
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
         enddo
      enddo
      do n=0, nterms
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndrufw(nterms, x, y, w, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)          resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nterms, x, w, nmax
cf2py intent(out) y
c
      integer(8) nterms,nmax,irat1,irat2
      double precision x, y(0:nterms,0:nterms), w(*)
c
      if( nterms .le. nmax ) then
        irat1=1
        irat2=1+(nmax+1)**2
        call ylgndrufw0(nterms, x, y, w(irat1), w(irat2), nmax)
      else
        call ylgndru(nterms, x, y)
      endif
c
      return
      end
c
c
c
c
c
      subroutine ylgndru2fw(nterms, x, y, d, w, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, w, nmax
cf2py intent(out) y, d
c
      integer(8) nterms,nmax,irat1,irat2
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision w(*)
c
      if( nterms .le. nmax ) then
        irat1=1
        irat2=1+(nmax+1)**2
        call ylgndru2fw0(nterms, x, y, d, w(irat1), w(irat2), nmax)
      else
        call ylgndru2(nterms, x, y, d)
      endif
c
      return
      end
c
c
c
c
c
      subroutine ylgndru2sfw(nterms, x, y, d, w, nmax)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, w, nmax
cf2py intent(out) y, d
c
      integer(8) nterms,nmax,irat1,irat2
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision w(*)
c
      if( nterms .le. nmax ) then
        irat1=1
        irat2=1+(nmax+1)**2
        call ylgndru2sfw0(nterms, x, y, d, w(irat1), w(irat2), nmax)
      else
        call ylgndru2s(nterms, x, y, d)
      endif
c
      return
      end
c
c
c
c
c
      subroutine ylgndrufw0_old(nterms, x, y, rat1, rat2, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y
c
      integer(8) nterms, nmax, n, m
      double precision x, y(0:nterms,0:nterms), u
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      do 10 m=0, nterms
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.lt.nterms)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 do 20 n=m+2, nterms
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
 20      continue
 10   continue
c
      return
      end
c
c
c
c
c
      subroutine ylgndrufw0(nterms, x, y, rat1, rat2, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y
c
      integer(8) nterms, nmax, n, m
      double precision x, y(0:nterms,0:nterms), u
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)

      y(0,0)=1
      if( nterms .eq. 0 ) return

      y(1,0)=x*y(0,0)*rat1(1,0)

      u=-sqrt((1-x)*(1+x))
      do m=1, nterms-1
	 y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      enddo

      y(nterms,nterms)=y(nterms-1,nterms-1)*u*rat1(nterms,nterms)

c      do m=0, nterms-1
c	 do n=m+2, nterms
c	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
c         enddo
c      enddo

      do n=2, nterms
         do m=0, n-2
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndru2fw0_old(nterms, x, y, d, rat1, rat2, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y, d
c
      integer(8) nterms, nmax, n, m
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision u, du
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      du=x/sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
      do 10 m=0, nterms
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x/u**2
	 if (m.lt.nterms)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nterms)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
	 do 20 n=m+2, nterms
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
 20      continue
 10   continue
c
      return
      end
c
c
c
c
c
      subroutine ylgndru2fw0(nterms, x, y, d, rat1, rat2, nmax)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y, d
c
      integer(8) nterms, nmax, n, m
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision u, u2
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)

      y(0,0)=1
      d(0,0)=0
      if( nterms .eq. 0 ) return

      y(1,0)=x*y(0,0)*rat1(1,0)
      d(1,0)=(x*d(0,0)+y(0,0))*rat1(1,0)

      u=-sqrt((1-x)*(1+x))
      u2=(1-x)*(1+x)
      do m=1, nterms-1
         y(m,m)=y(m-1,m-1)*u*rat1(m,m)
         d(m,m)=y(m,m)*(-m)*x/u2
         y(m+1,m)=x*y(m,m)*rat1(m+1,m)
         d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
      enddo

      y(nterms,nterms)=y(nterms-1,nterms-1)*u*rat1(nterms,nterms)
      d(nterms,nterms)=y(nterms,nterms)*(-nterms*x)

      do n=2, nterms
         do m=0, n-2
         y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
         d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndru2sfw0_old(nterms, x, y, d, rat1, rat2, nmax)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y, d
c
      integer(8) nterms, nmax, n, m
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision u, u2
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nterms)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      if (m.lt.nterms)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
      do 120 n=m+2, nterms
        y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
120   continue
c
c       ... then, evaluate scaled associated Legendre functions
c
      do 210 m=1, nterms
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
c
	 if (m.lt.nterms)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nterms)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
	 do 220 n=m+2, nterms
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
220      continue
210   continue
c
      return
      end
c
c
c

      subroutine ylgndru2sfw0(nterms, x, y, d, rat1, rat2, nmax)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c
c      Ynm(x) = sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nterms
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nterms                      must be non-negative
c     x                         -1 <= x <= 1
c     y(0:nterms,0:nterms)      resulting function values
c     d(0:nterms,0:nterms)      resulting derivative values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nterms and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nterms, x, rat1, rat2, nmax
cf2py intent(out) y, d
c
      integer(8) nterms, nmax, n, m
      double precision x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
      double precision u, u2
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)

      y(0,0)=1
      d(0,0)=0
      if( nterms .eq. 0 ) return

      y(1,0)=x*y(0,0)*rat1(1,0)
      d(1,0)=(x*d(0,0)+y(0,0))*rat1(1,0)

      u=-sqrt((1-x)*(1+x))
      u2=(1-x)*(1+x)
      do m=1, nterms-1
         if( m .eq. 1 ) y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
         if( m .gt. 1 ) y(m,m)=y(m-1,m-1)*u*rat1(m,m)
         d(m,m)=y(m,m)*(-m)*x
         y(m+1,m)=x*y(m,m)*rat1(m+1,m)
         d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
      enddo

      y(nterms,nterms)=y(nterms-1,nterms-1)*u*rat1(nterms,nterms)
      d(nterms,nterms)=y(nterms,nterms)*(-nterms*x)

      do n=2, nterms
         m=0
         y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
         d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
         do m=1, n-2
         y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
         d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      subroutine ylgndr2smi(nmax, x, y, d, id)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     For m>0, 
c          the functions are scaled by 1/sqrt(1-x**2)**m
c          the derivatives are scaled by 1/sqrt(1-x**2)**(m-2)
c
c      scaled by 10^(-id) to handle underflow/overflow
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c     id(0:nmax,0:nmax)     external scaling factor 10^(-id)
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x, id
cf2py intent(out) y,d
c
      integer(8) nmax, m, n, id(0:nmax,0:nmax)
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, u2
c
      u=-sqrt((1-x))*sqrt((1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*sqrt(2*m+1.0d0)
      id(0,m)=0
      id(1,m)=0
      do n=m+2, nmax
        y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        d(n,m)=((2*n-1)*(x*d(n-1,m)+y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
        id(n,m)=id(n-1,m)
        if( abs(y(n,m)) .gt. 1d+50) then
        y(n,m)=y(n,m)/1d50
        d(n,m)=d(n,m)/1d50
        id(n,m)=id(n,m)+50
        y(n-1,m)=y(n-1,m)/1d50
        d(n-1,m)=d(n-1,m)/1d50
        id(n-1,m)=id(n-1,m)+50
        endif
        if( abs(y(n,m)) .lt. 1d-50) then
        y(n,m)=y(n,m)*1d50
        d(n,m)=d(n,m)*1d50
        id(n,m)=id(n,m)-50
        y(n-1,m)=y(n-1,m)*1d50
        d(n-1,m)=d(n-1,m)*1d50
        id(n-1,m)=id(n-1,m)-50
        endif
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*(-1)*sqrt((2*m-1.0d0)/(2*m))
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
	 if (m.gt.0)  id(m,m)=id(m-1,m-1)
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*sqrt(2*m+1.0d0)
	 if (m.lt.nmax)  id(m+1,m)=id(m,m)

        if( abs(y(m,m)) .gt. 1d+50) then
        y(m,m)=y(m,m)/1d50
        d(m,m)=d(m,m)/1d50
        id(m,m)=id(m,m)+50
        if (m.lt.nmax) then
        y(m+1,m)=y(m+1,m)/1d50
        d(m+1,m)=d(m+1,m)/1d50
        id(m+1,m)=id(m+1,m)+50
        endif
        endif
c
        if( abs(y(m,m)) .lt. 1d-50) then
        y(m,m)=y(m,m)*1d50
        d(m,m)=d(m,m)*1d50
        id(m,m)=id(m,m)-50
        if (m.lt.nmax) then
        y(m+1,m)=y(m+1,m)*1d50
        d(m+1,m)=d(m+1,m)*1d50
        id(m+1,m)=id(m+1,m)-50
        endif
        endif
c

	 do n=m+2, nmax
	    y(n,m)=((2*n-1)*x*y(n-1,m) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*y(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
	    d(n,m)=((2*n-1)*(x*d(n-1,m)+u2*y(n-1,m)) - 
     1               sqrt((n+m-1.0d0)*(n-m-1.0d0))*d(n-2,m))
     2               /sqrt((n-m+0.0d0)*(n+m))
            id(n,m)=id(n-1,m)

        if( abs(y(n,m)) .gt. 1d+50) then
        y(n,m)=y(n,m)/1d50
        d(n,m)=d(n,m)/1d50
        id(n,m)=id(n,m)+50
        y(n-1,m)=y(n-1,m)/1d50
        d(n-1,m)=d(n-1,m)/1d50
        id(n-1,m)=id(n-1,m)+50
        endif
        if( abs(y(n,m)) .lt. 1d-50) then
        y(n,m)=y(n,m)*1d50
        d(n,m)=d(n,m)*1d50
        id(n,m)=id(n,m)-50
        y(n-1,m)=y(n-1,m)*1d50
        d(n-1,m)=d(n-1,m)*1d50
        id(n-1,m)=id(n-1,m)-50
        endif

         enddo
      enddo
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       real valued Legendre functions, with scaling factor
c       to handle underflow near the poles
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
      subroutine ylgndrfe(nmax, x, y, rat1, rat2, id)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      scaled by 10^(-id) to handle underflow/overflow
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     id(0:nmax,0:nmax)     external scaling factor 10^(-id)
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nmax, x, rat1, rat2, id
cf2py intent(out) y
c
      integer(8) nmax, m, n, id(0:nmax,0:nmax)
      double precision x, y(0:nmax,0:nmax), u
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      y(0,0)=1
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  id(m,m)=id(m-1,m-1)
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  id(m+1,m)=id(m,m)

        if( abs(y(m,m)) .gt. 1d+50) then
        y(m,m)=y(m,m)/1d50
        id(m,m)=id(m,m)+50
        if (m.lt.nmax) then
        y(m+1,m)=y(m+1,m)/1d50
        id(m+1,m)=id(m+1,m)+50
        endif
        endif
c
        if( abs(y(m,m)) .lt. 1d-50) then
        y(m,m)=y(m,m)*1d50
        id(m,m)=id(m,m)-50
        if (m.lt.nmax) then
        y(m+1,m)=y(m+1,m)*1d50
        id(m+1,m)=id(m+1,m)-50
        endif
        endif

	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
            id(n,m)=id(n-1,m)

        if( abs(y(n,m)) .gt. 1d+50) then
        y(n,m)=y(n,m)/1d50
        id(n,m)=id(n,m)+50
        y(n-1,m)=y(n-1,m)/1d50
        id(n-1,m)=id(n-1,m)+50
        endif
        if( abs(y(n,m)) .lt. 1d-50) then
        y(n,m)=y(n,m)*1d50
        id(n,m)=id(n,m)-50
        y(n-1,m)=y(n-1,m)*1d50
        id(n-1,m)=id(n-1,m)-50
        endif

         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2fe(nmax, x, y, d, rat1, rat2, id)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      scaled by 10^(-id) to handle underflow/overflow
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c     id(0:nmax,0:nmax)     external scaling factor 10^(-id)
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x, rat1, rat2, id
cf2py intent(out) y, d
c
      integer(8) nmax, m, n, id(0:nmax,0:nmax)
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
c
      u=-sqrt((1-x)*(1+x))
      du=x/sqrt((1-x)*(1+x))
      y(0,0)=1
      d(0,0)=0
      id(0,0)=0
      do m=0, nmax
	 if (m.gt.0)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x/u**2
	 if (m.gt.0)  id(m,m)=id(m-1,m-1)
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
	 if (m.lt.nmax)  id(m+1,m)=id(m,m)

        if( abs(y(m,m)) .gt. 1d+50) then
        y(m,m)=y(m,m)/1d50
        d(m,m)=d(m,m)/1d50
        id(m,m)=id(m,m)+50
        if (m.lt.nmax) then
        y(m+1,m)=y(m+1,m)/1d50
        d(m+1,m)=d(m+1,m)/1d50
        id(m+1,m)=id(m+1,m)+50
        endif
        endif
c
        if( abs(y(m,m)) .lt. 1d-50) then
        y(m,m)=y(m,m)*1d50
        d(m,m)=d(m,m)*1d50
        id(m,m)=id(m,m)-50
        if (m.lt.nmax) then
        y(m+1,m)=y(m+1,m)*1d50
        d(m+1,m)=d(m+1,m)*1d50
        id(m+1,m)=id(m+1,m)-50
        endif
        endif

	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
            id(n,m)=id(n-1,m)

        if( abs(y(n,m)) .gt. 1d+50) then
        y(n,m)=y(n,m)/1d50
        d(n,m)=d(n,m)/1d50
        id(n,m)=id(n,m)+50
        y(n-1,m)=y(n-1,m)/1d50
        d(n-1,m)=d(n-1,m)/1d50
        id(n-1,m)=id(n-1,m)+50
        endif
        if( abs(y(n,m)) .lt. 1d-50) then
        y(n,m)=y(n,m)*1d50
        d(n,m)=d(n,m)*1d50
        id(n,m)=id(n,m)-50
        y(n-1,m)=y(n-1,m)*1d50
        d(n-1,m)=d(n-1,m)*1d50
        id(n-1,m)=id(n-1,m)-50
        endif

         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2sfe(nmax, x, y, d, rat1, rat2, id)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c      scaled by 10^(-id) to handle underflow/overflow
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c     id(0:nmax,0:nmax)     external scaling factor 10^(-id)
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x, rat1, rat2, id
cf2py intent(out) y, d
c
      integer(8) nmax, m, n, id(0:nmax,0:nmax)
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, u2
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
      u=-sqrt((1-x)*(1+x))
      u2 = (1-x)*(1+x)
      y(0,0)=1
      d(0,0)=0
      id(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
      if (m.lt.nmax)  d(m+1,m)=(x*d(m,m)+y(m,m))*rat1(m+1,m)
      if (m.lt.nmax)  id(m+1,m)=id(m,m)

      do n=m+2, nmax
        y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(x*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
        id(n,m)=id(n-1,m)

        if( abs(y(n,m)) .gt. 1d+50) then
        y(n,m)=y(n,m)/1d50
        d(n,m)=d(n,m)/1d50
        id(n,m)=id(n,m)+50
        y(n-1,m)=y(n-1,m)/1d50
        d(n-1,m)=d(n-1,m)/1d50
        id(n-1,m)=id(n-1,m)+50
        endif
        if( abs(y(n,m)) .lt. 1d-50) then
        y(n,m)=y(n,m)*1d50
        d(n,m)=d(n,m)*1d50
        id(n,m)=id(n,m)-50
        y(n-1,m)=y(n-1,m)*1d50
        d(n-1,m)=d(n-1,m)*1d50
        id(n-1,m)=id(n-1,m)-50
        endif

      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
	 if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
	 if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
	 if (m.gt.0)  d(m,m)=y(m,m)*(-m)*x
	 if (m.gt.0)  id(m,m)=id(m-1,m-1)
c
	 if (m.lt.nmax)  y(m+1,m)=x*y(m,m)*rat1(m+1,m)
	 if (m.lt.nmax)  
     $      d(m+1,m)=(x*d(m,m)+u2*y(m,m))*rat1(m+1,m)
	 if (m.lt.nmax)  id(m+1,m)=id(m,m)

        if( abs(y(m,m)) .gt. 1d+50) then
        y(m,m)=y(m,m)/1d50
        d(m,m)=d(m,m)/1d50
        id(m,m)=id(m,m)+50
        if (m.lt.nmax) then
        y(m+1,m)=y(m+1,m)/1d50
        d(m+1,m)=d(m+1,m)/1d50
        id(m+1,m)=id(m+1,m)+50
        endif
        endif
c
        if( abs(y(m,m)) .lt. 1d-50) then
        y(m,m)=y(m,m)*1d50
        d(m,m)=d(m,m)*1d50
        id(m,m)=id(m,m)-50
        if (m.lt.nmax) then
        y(m+1,m)=y(m+1,m)*1d50
        d(m+1,m)=d(m+1,m)*1d50
        id(m+1,m)=id(m+1,m)-50
        endif
        endif


	 do n=m+2, nmax
	    y(n,m)=rat1(n,m)*x*y(n-1,m)-rat2(n,m)*y(n-2,m)
	    d(n,m)=rat1(n,m)*(x*d(n-1,m)+u2*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
            id(n,m)=id(n-1,m)

        if( abs(y(n,m)) .gt. 1d+50) then
        y(n,m)=y(n,m)/1d50
        d(n,m)=d(n,m)/1d50
        id(n,m)=id(n,m)+50
        y(n-1,m)=y(n-1,m)/1d50
        d(n-1,m)=d(n-1,m)/1d50
        id(n-1,m)=id(n-1,m)+50
        endif
        if( abs(y(n,m)) .lt. 1d-50) then
        y(n,m)=y(n,m)*1d50
        d(n,m)=d(n,m)*1d50
        id(n,m)=id(n,m)-50
        y(n-1,m)=y(n-1,m)*1d50
        d(n-1,m)=d(n-1,m)*1d50
        id(n-1,m)=id(n-1,m)-50
        endif

         enddo
      enddo
c
      do n=0, nmax
	 do m=0, n
	    y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
	    d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end
c
c
c
c
c
      subroutine ylgndrfex(nmax, x, y, rat1, rat2, id)
      implicit none
c
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     id(0:nmax,0:nmax)     internal scaling factor 10^(-id)
c
c     In this routine, intermediary values are scaled 
c     by 10^(-id) to handle recurrence underflow/overflow
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. 
c
cf2py intent(in) nmax, x, rat1, rat2, id
cf2py intent(out) y
c
      integer(8) nmax, m, n, id(0:nmax,0:nmax)
      double precision x, y(0:nmax,0:nmax), u, sc
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)

      call ylgndrfe(nmax, x, y, rat1, rat2, id)

      do m=0, nmax
         do n=m, nmax
            sc = 10d0**(id(n,m))
	    y(n,m)=y(n,m)*sc
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2fex(nmax, x, y, d, rat1, rat2, id)
      implicit none
c
c     Evaluate normalized Legendre functions and their derivatives
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c     id(0:nmax,0:nmax)     internal scaling factor 10^(-id)
c
c     In this routine, intermediary values are scaled 
c     by 10^(-id) to handle recurrence underflow/overflow
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x, rat1, rat2, id
cf2py intent(out) y, d
c
      integer(8) nmax, m, n, id(0:nmax,0:nmax)
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, du, sc
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
c
      call ylgndr2fe(nmax, x, y, d, rat1, rat2, id)

      do m=0, nmax
         do n=m, nmax
            sc = 10d0**(id(n,m))
	    y(n,m)=y(n,m)*sc
	    d(n,m)=d(n,m)*sc 
         enddo
      enddo
c
      return
      end
c
c
c
c
c
      subroutine ylgndr2sfex(nmax, x, y, d, rat1, rat2, id)
      implicit none
c
c     Evaluate scaled normalized Legendre functions and their derivatives
c
c     Only for Ynm(x) with m>0 
c          the functions are scaled by 1/sqrt(1-x**2)
c          the derivatives are scaled by sqrt(1-x**2)
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c      d Ynm(x) / dx = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) d Pnm(x) / dx
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c     d(0:nmax,0:nmax)      resulting derivative values
c     id(0:nmax,0:nmax)     internal scaling factor 10^(-id)
c
c     In this routine, intermediary values are scaled 
c     by 10^(-id) to handle recurrence underflow/overflow
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values. The same convention for the derivatives.
c
cf2py intent(in) nmax, x, rat1, rat2, id
cf2py intent(out) y, d
c
      integer(8) nmax, m, n, id(0:nmax,0:nmax)
      double precision x, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u, u2, sc
      double precision rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)

      call ylgndr2sfe(nmax, x, y, d, rat1, rat2, id)

      do m=0, nmax
         do n=m, nmax
            sc = 10d0**(id(n,m))
	    y(n,m)=y(n,m)*sc
	    d(n,m)=d(n,m)*sc
         enddo
      enddo
c
      return
      end
c
c
c
c
c
