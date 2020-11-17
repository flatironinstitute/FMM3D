c
c      read plane wave quadratures from file
c

      subroutine lreadall(eps,zk,nquad,cxs,cws,nfour,nphys,lw,ier)
c
cc        this subroutine computes lambda quadrature nodes and weights
c         for given helmholtz parameter and the number of terms
c         per lambda required in fourier and physical domains
c
c         input:
c         eps:   precision requested 
c                  min precision possible is 0.5d-12
c
c
c         zk -      Helmholtz parameter
c                   currently code works for zk in
c                    [0.02,pi/2]+[0,0.02]i
c  
c         lw -      length of various arrays sent to the routine
c
c         OUTPUT
c         nquad -    number of discretization nodes for lambda
c                    quadrature
c
c         cxs(nquad) -  complex *16, lambda quadrature nodes, 
c                       cxs actually stores \sqrt(\lambda^2 - k^2)
c
c         cws(nquad) - complex *16, lambda quadrature weights
c
c         nfour(nquad) - integer(8), number of fourier modes required
c                        per lambda
c
c         nphys(nquad) - integer(8), number of physical modes required
c                        per lambda
c
c         ier - error code
c               ier = 0, no error, successful return of all parameters
c               ier = 16 if lw<nquad
c
c
       
        implicit real *8 (a-h,o-z)
        complex *16 cxs(*),cws(*),zk
        real *8 cx(100),cy(100),eps
        integer(8) nfour(*),nphys(*)

        done = 1
        pi = atan(done)*4



c
        nnx=11
        nny=11
c
        cx(1)=0.00d0
        cx(2)=0.02d0
        cx(3)=0.04d0
        cx(4)=0.10d0
        cx(5)=0.20d0
        cx(6)=0.40d0
        cx(7)=1.00d0
        cx(8)=pi/2
        cx(9)=pi
        cx(10)=2*pi
        cx(11)=3*pi
        cx(12)=4*pi
        cx(13)=5*pi
        cx(14)=6*pi
        cx(15)=7*pi
        cx(16)=8*pi
        cx(17)=9*pi
        cx(18)=10*pi
        cx(19)=11*pi
        cx(20)=12*pi
        cx(21)=13*pi
        cx(22)=14*pi
        cx(23)=15*pi
        cx(24)=16*pi
c
        cy(1)=0.00d0
        cy(2)=0.02d0
        cy(3)=0.04d0
        cy(4)=0.10d0
        cy(5)=0.20d0
        cy(6)=0.40d0
        cy(7)=1.00d0
        cy(8)=pi/2
        cy(9)=pi
        cy(10)=2*pi
        cy(11)=3*pi
        cy(12)=4*pi
        cy(13)=5*pi
        cy(14)=6*pi
        cy(15)=7*pi
        cy(16)=8*pi
        cy(17)=9*pi
        cy(18)=10*pi
        cy(19)=11*pi
        cy(20)=12*pi
        cy(21)=13*pi
        cy(22)=14*pi
        cy(23)=15*pi
        cy(24)=16*pi


        ier = 0

        call hwts3getd(ier,zk,nnx,nny,cx,cy,iquad,ix,iy)

        iprec = 1
        if(eps.lt.0.5d-3) iprec = 2
        if(eps.lt.0.5d-6) iprec = 3
        if(eps.lt.0.5d-9) iprec = 4

        call hwts3(ier,iprec,zk,cxs,cws,nquad)
        
        if(nquad.gt.lw) then

        ier = 16
        call prin2('not enough memory for quad nodes=*',ier,1)
        stop

        endif

        call numphysfour(iprec,iquad,nquad,nphys,nfour)

        return
        end
 
 
