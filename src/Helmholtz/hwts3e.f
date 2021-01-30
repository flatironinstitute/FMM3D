ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the exponential representation routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine hwts3e(ier,eps,rk,cxs,cws,n)
        implicit real *8 (a-h,o-z)
        integer(8) ier,n,idomain
        complex *16 cxs(*),cws(*)
        complex *16 rk,ima
        real *8, allocatable :: xs(:),ws(:)
        complex *16, allocatable :: xs1(:),ws1(:)
        data ima/(0.0d0,1.0d0)/
c
c       
c       This subroutine returns the nodes and weights for 
c       the exponential representation 
c
c       exp(ima*rk*sqrt(u^2+v^2)/sqrt(u^2+v^2) = 
c       \int_0^{+\infty} exp(-sqrt(t^2-rk^2)*t) J_0(tv) t/sqrt(t^2-rk^2) dt
c
c       u in [1,4], v in [0, 4*sqrt(2)]
c
c       All tables are valid in 
c
c       rk range [0, 16*pi] \times [0, 12*pi]
c
c
c
c       Input parameter:
c       
c       iprec - precision required:
c           iprec = 0 - 2 digits
c           iprec = 1 - 3 digits
c           iprec = 2 - 6 digits
c           iprec = 3 - 9 digits
c           iprec = 4 - 12 digits
c
c       Output parameters:
c
c       n - the number of nodes
c       cxs - the nodes: complex *16 (n)
c       cws - the weights: complex *16 (n)
c

        lxs = 10000
        allocate(xs(lxs),ws(lxs),xs1(lxs),ws1(lxs))

        n = 0
c
        ier = 0
c       
        call hwts3dgetd(ier,rk,idomain)
c
ccc        call prinf('ier=*',ier,1)
ccc        write(*,*) 'idomain=*', idomain
ccc        call prinf('idomain=*',idomain,1)
c

        iprec = 0
        if(eps.lt.0.5d-2) iprec = 1
        if(eps.lt.0.5d-3) iprec = 2
        if(eps.lt.0.5d-6) iprec = 3
        if(eps.lt.0.5d-9) iprec = 4
        

        if( iprec .eq. 0 ) errmax = .1d-2
        if( iprec .eq. 1 ) errmax = .1d-3
        if( iprec .eq. 2 ) errmax = .1d-6
        if( iprec .eq. 3 ) errmax = .1d-9
        if( iprec .eq. 4 ) errmax = .1d-12
        
        if( iprec .eq. 0 ) call hwts3p0(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 1 ) call hwts3p1(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 2 ) call hwts3p2(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 3 ) call hwts3p3(idomain,xs,ws,n,iquadtype,err)
        if( iprec .eq. 4 ) call hwts3p4(idomain,xs,ws,n,iquadtype,err)

 1100   continue
        
ccc        write(*,*) 'iquad=*', iquad
ccc        call prinf('iquad=*',iquad,1)
ccc        pause
c
c
cc        call prinf('n=*',n,1)
cc        call prin2('xs=*',xs,n)
cc        call prin2('ws=*',ws,n)
cc        call prinf('iquadtype=*',iquadtype,1)
c
        do 1200 k=1,n
c
        if( iquadtype .eq. 1 ) then
        cxs(k)=sqrt(xs(k)**2-rk**2)
        cws(k)=ws(k)*xs(k)/sqrt(xs(k)**2-rk**2)
        endif
        
        if( iquadtype .eq. 2 ) then
        xmin=0
        xmax=100
        call linmap0(xs(k),u,uweight,xmin,xmax)
ccc           if( dble(rk) .gt. 1d0 ) then
ccc           xmid=dble(rk)
ccc           call linmap1(xs(k),u,uweight,xmin,xmax,xmid)
ccc           endif
c
        gamma=0
        v=u/(1+gamma*abs(dble(-ima*rk)))
        cxs(k)=u-ima*rk+abs(imag(-ima*rk))*ima*(v/(1+v))
        cws(k)=1+abs(imag(-ima*rk))*ima*((1+v)-v)/(1+v)**2
     $     /(1+gamma*abs(dble(-ima*rk)))
        cws(k)=cws(k)*ws(k)*uweight
        endif
c
        if( iquadtype .eq. 5 ) then
        xs1(k)=xs(k)-ima*atan(xs(k))
        ws1(k)=ws(k)*(1-ima/(1+xs(k)**2))
        cxs(k)=sqrt(xs1(k)**2-rk**2)
        cws(k)=ws1(k)*xs1(k)/sqrt(xs1(k)**2-rk**2)
        endif
c
 1200   continue
c
        return
        end

c
c
c
c
        subroutine linmap0(t,u,w,ax,bx)
        implicit real *8 (a-h,o-z)
c
c       this subroutine contructs a linear mapping from interval
c       [0,1] to [ax,bx] 
c        
        u=ax+t*(bx-ax)
        w=(bx-ax)
c
        return
        end
c
c
c
c
c
        subroutine linmap1(t,u,w,ax,bx,cx)
        implicit real *8 (a-h,o-z)
c
c       this subroutine contructs a piecewise linear mapping from interval
c       [0,1] to [ax,bx], such that [0,0.5]->[ax,cx], [0.5,1]->[cx,bx]
c        
        if( t .le. 0.5d0 ) then
        u=ax+2*t*(cx-ax)
        w=2*(cx-ax)
        endif
c
        if( t .gt. 0.5d0 ) then
        u=cx+2*(t-0.5d0)*(bx-cx)
        w=2*(bx-cx)
        endif
c
        return
        end
c
c
c
c
c

        subroutine hwts3dgetd(ier,rk,idomain)
        implicit real *8 (a-h,o-z)
        integer(8) ier, idomain
        complex *16 ima,rk
        dimension cx(100),cy(100)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
c       rk range [0, 16*pi] \times [0, 12*pi]
c
        nx=23
        ny=19
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
c
        rkrea=dble(rk)
        rkima=imag(rk)
c
        idomain = 0
        ier = 0
c
        ix=0
        do 1200 i=1,nx
        if( rkrea .ge. cx(i+1) ) ix=i+1
 1200   continue
c
 1300   continue
c
        iy=0
        do 2200 i=1,ny
        if( rkima .ge. cy(i+1) ) iy=i+1
 2200   continue
c
 2300   continue
c
ccc        write(*,*) rk
ccc        write(*,*) ix,iy
c
        if( ix .gt. nx ) ier = 4
        if( iy .gt. ny ) ier = 4
        if( ix .gt. nx ) return
        if( iy .gt. ny ) return
c       
        if( ix .eq. 0 ) ix = 1
        if( iy .eq. 0 ) iy = 1
c
        if( ix .gt. nx ) ix = nx
        if( iy .gt. ny ) iy = ny
c
        idomain = ix + (iy-1) * nx
c
c        write(*,*) ix,iy
cc        write(*,*) idomain
c
ccc        stop
c
        return
        end
c
c
c
c
        subroutine getrkinfo2(inum,ndomains,rkmin,rkmax)
        implicit real *8 (a-h,o-z)
        complex *16 rkmin,rkmax,ima
        complex *16 cx(100),cy(100)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
c       rk range [0, 16*pi] \times [0, 12*pi]
c
        nx=23
        ny=19
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
c
        ndomains=nx*ny
c
        ix=mod(inum-1,nx)+1
        iy=(inum-1)/nx+1
c
        rkmin=cx(ix)+ima*cy(iy)
        rkmax=cx(ix+1)+ima*cy(iy+1)
c
        return
        end
c
c
c
c
c
        subroutine hwts3p0(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        integer(8) idomain,n
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_3d_iprec0.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts3p1(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        integer(8) idomain,n
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_3d_iprec1.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts3p2(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        integer(8) idomain,n
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_3d_iprec2.txt'
c
        return
        end        
c        
c
c
c
c
        subroutine hwts3p3(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        integer(8) idomain,n
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_3d_iprec3.txt'
c
        return
        end
c        
c
c
c
c
        subroutine hwts3p4(idomain,xs,ws,n,iquadtype,err)
        implicit real *8 (a-h,o-z)
        integer(8) idomain,n
        dimension xs(*), ws(*)
c
        n = 0
        iquadtype = 0
c
        INCLUDE 'hwts_3d_iprec4.txt'
c
        return
        end
c
c
c
c
c
