c     testing code for FMM - tests charges and dipoles against
c     O(N^2) direct method
c
c
        implicit double precision (a-h,o-z)
        call test_part()
c
        stop
        end
c
c---------------------------------------------------------------------------
c
        subroutine test_part()
        implicit double precision (a-h,o-z)

        parameter(nsmax=1 000 000)
        parameter(ntmax=8 000 000)
        parameter(ntest = 1000)
        double precision xtmp(nsmax),ytmp(nsmax)

        dimension source(3,nsmax)
        double complex charge(nsmax)
        double complex dipstr(nsmax)
        dimension dipvec(3,nsmax)
        double complex pot(nsmax)
        double complex fld(3,nsmax)
        
c
        dimension target(3,ntmax)
        double complex pottarg(ntmax)
        double complex fldtarg(3,ntmax)

        double complex potex(ntest),fldex(3,ntest)
        double complex pottargex(ntest),fldtargex(3,ntest)

c
c
        double complex ima
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
c
cc         distributions: idist = 1, sources and targets in the 
c                                    volume of a cube
c                         idist = 2, sources and targets on the 
c                                    surface of a sphere
c                         idist = 3, sources on surface of a sphere
c                                    targets in the volume
c                                    (number of targets = nsources^1.5)
c                         idist = 4, sources on surface of a cube
c                                    with a dyadic mesh
c  


        call prin2('Enter n*',i,0)
        read *, n

        icase = 2
        idist= 2
        iprec = 2

        ifcharge = 1
        ifdipole = 0

        ifpot = 0
        iffld = 0

        ifpottarg = 1
        iffldtarg = 0


c
c
        if( idist .eq. 1 ) then
c
c       ... construct randomly located charge distribution on a unit 
c            cube
c
           nsource = icase*100000
           ntarget = nsource
           do i=1,nsource

              source(1,i)=hkrand(0)-0.5d0
              source(2,i)=hkrand(0)-0.5d0
              source(3,i)=hkrand(0)-0.5d0

              target(1,i) = hkrand(0)-0.5d0
              target(2,i) = hkrand(0)-0.5d0
              target(3,i) = hkrand(0)-0.5d0
           enddo
        endif
c
        if( idist .eq. 2 ) then
c
c       ... construct randomly located charge distribution on 
c             surface of sphere
c

           nsource = icase*100000

           ntarget = nsource
           r = 0.5d0
           do i=1,nsource
              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              source(1,i) = r*cos(thet)*cos(phi)
              source(2,i) = r*cos(thet)*sin(phi)
              source(3,i) = r*sin(thet)

              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              target(1,i) = r*cos(thet)*cos(phi)
              target(2,i) = r*cos(thet)*sin(phi)
              target(3,i) = r*sin(thet)

           enddo
        endif
c
        if( idist .eq. 3 ) then
c
c       ... construct randomly located charge distribution on 
c            a unit sphere and targets in the volume
c
          
           nsource = icase*10000 
           ntarget = nsource**(1.5d0)/5
           r = 0.5d0
           do i=1,nsource
              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              source(1,i) = r*cos(thet)*cos(phi)
              source(2,i) = r*cos(thet)*sin(phi)
              source(3,i) = r*sin(thet)

           enddo

           do i=1,ntarget
              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              rr=  hkrand(0)*r

              target(1,i) = rr*cos(thet)*cos(phi)
              target(2,i) = rr*cos(thet)*sin(phi)
              target(3,i) = rr*sin(thet)
           enddo
        endif

        if(idist.eq.4) then

c
c         ... consturct two face with dyadic refinement 
c         along edges and vertices
c         make sure face is of random length so that boxes 
c         and face grid doesn't align
           
           k = 16
           rpan = hkrand(0)

           nlev = icase*10


           call getface(nlev,rpan,k,xtmp,ytmp)
           nsource = (k*nlev)**2

           a1 = hkrand(0)
           b1 = hkrand(0)
           
           a2 = hkrand(0)
           b2 = hkrand(0)

           a3 = hkrand(0)
           b3 = hkrand(0)

           ntarget = nsource

           do i=1,nsource
              source(1,i) = a1*xtmp(i) + b1*ytmp(i) 
              source(2,i) = a2*xtmp(i) + b2*ytmp(i)
              source(3,i) = a3*xtmp(i) + b3*ytmp(i)

              target(1,i) = source(1,i)
              target(2,i) = source(2,i)
              target(3,i) = source(3,i)
           enddo

        endif


        do i=1,nsource
           charge(i) = 1
           dipstr(i) = 1
           dipvec(1,i) = hkrand(0)
           dipvec(2,i) = hkrand(0)
           dipvec(3,i) = hkrand(0)
        enddo

        do i=1,ntarget
           pottarg(i) = 0
           fldtarg(1,i) = 0
           fldtarg(2,i) = 0
           fldtarg(3,i) = 0
        enddo

        t1 = second()
c       Call memory management tree code

        if(idist.eq.1) then 
           if(iprec.eq.1) then
              idivflag = 3
              ndiv = 100
           endif

           if(iprec.eq.2) then
              idivflag = 3
              ndiv =150 
           endif
        endif

        if(idist.eq.2) then
           if(iprec.eq.2) then
              idivflag = 0
              ndiv = 200
           endif
        endif

        if(idist.eq.4) then
           if(iprec.eq.2) then
              idivflag = 0
              ndiv = 200
           endif
        endif

        do i=1,nsource
           pot(i) =0
           fld(1,i) = 0
           fld(2,i) = 0
           fld(3,i) = 0
        enddo

        do i=1,ntarget
           pottarg(i) =0
           fldtarg(1,i) = 0
           fldtarg(2,i) = 0
           fldtarg(3,i) = 0
        enddo


        call lfmm3dpartstost(ier,iprec,nsource,source,ifcharge,
     1          charge,ifdipole,dipstr,dipvec,ifpot,pot,
     2          iffld,fld,ntarget,target,ifpottarg,pottarg,
     3          iffldtarg,fldtarg)


        m=min(nsource,100)

        call geterr(nsource,source,ifcharge,charge,ifdipole,dipstr,
     1               dipvec,m,source,ifpot,pot,iffld,
     2               fld,rerr,rgerr)

        call prin2('source potential error=*',rerr,1)
        call prin2('source field error=*',rgerr,1)

        call geterr(nsource,source,ifcharge,charge,ifdipole,dipstr,
     1               dipvec,m,target,ifpottarg,pottarg,iffldtarg,
     2               fldtarg,rerr,rgerr)
        call prin2('target pot error=*',rerr,1)
        call prin2('target field error=*',rgerr,1)

        
        return
        end
     
c-----------------------------------------------------------------

        subroutine geterr(nsource,source,ifcharge,charge,ifdipole,
     1               dipstr,dipvec,m,target,ifpottarg,pottarg,
     2               iffldtarg,fldtarg,rerr,rgerr)

        implicit double precision (a-h,o-z)
        double precision source(3,*),target(3,*),dipvec(3,*)
        double complex charge(*),dipstr(*),pottarg(*),fldtarg(3,*)
        double complex pot2(m),fld2(3,m),ptemp,ftemp(3)

        rerr = 0
        rgerr = 0


        do i=1,m
        if (ifpottarg .eq. 1) pot2(i)=0
        if (iffldtarg .eq. 1) then
           fld2(1,i)=0
           fld2(2,i)=0
           fld2(3,i)=0
        endif
        enddo
c
        t1=second()
        do 8160 j=1,m
        do 8150 i=1,nsource


        rr = (source(1,i)-target(1,j))**2 + 
     2          (source(2,i)-target(2,j))**2 + 
     3          (source(3,i)-target(3,j))**2

        if(rr.gt.1.0d-28) then

        if( ifcharge .eq. 1 ) then
        call lpotfld3d(iffldtarg,source(1,i),charge(i),
     $     target(1,j),
     1     ptemp,ftemp)
        if (ifpottarg .eq. 1) pot2(j)=pot2(j)+ptemp
        if (iffldtarg .eq. 1) then
           fld2(1,j)=fld2(1,j)+ftemp(1)
           fld2(2,j)=fld2(2,j)+ftemp(2)
           fld2(3,j)=fld2(3,j)+ftemp(3)
        endif
        endif
        if (ifdipole .eq. 1) then
           call lpotfld3d_dp(iffldtarg,source(1,i),
     $     dipstr(i),dipvec(1,i),
     $     target(1,j),ptemp,ftemp)
           if (ifpottarg .eq. 1) pot2(j)=pot2(j)+ptemp
           if (iffldtarg .eq. 1) then
              fld2(1,j)=fld2(1,j)+ftemp(1)
              fld2(2,j)=fld2(2,j)+ftemp(2)
              fld2(3,j)=fld2(3,j)+ftemp(3)
           endif
        endif

        endif
c
 8150   continue
 8160   continue
c
        t2=second()

        
C
        call l3derror(pottarg,pot2,m,aerr,rerr)

        call prin2('relative L2 error in target potential=*',rerr,1)

     
        rgerr = 0
        call l3derror(fldtarg,fld2,3*m,aerr,rgerr)
        call prin2('relative L2 error in target gradient=*',rgerr,1)

        return
        end
c
c------------------------------------------------------
c
c
        subroutine l3derror(pot1,pot2,n,ae,re)
        implicit double precision (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        double complex pot1(n),pot2(n)
c
        d=0
        a=0
c

 1111 format(2(2x,d22.16))
        do i=1,n
        d=d+abs(pot1(i)-pot2(i))**2
        a=a+abs(pot1(i))**2
        enddo
c
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c
        ae=d
        re=d/a
c
        return
        end
c
c-------------------------------------------------
c

        subroutine getface(nlev,rpan,k,xtmp,ytmp)
        implicit double precision (a-h,o-z)
        double precision xtmp(*),ytmp(*),ts(k),uk(k*k),vk(k*k),wts(k)
        double precision, allocatable :: x(:)
        
        itype = 2
        call legeexps(itype,k,ts,uk,vk,wts)

        allocate(x(nlev*k))


        hcur = rpan

        ii = 1
        do i=1,nlev

           if(i.ne.nlev) tstart = hcur/2
           if(i.eq.nlev) tstart = 0
           tend = hcur
           do j=1,k
              x(ii) = tstart + (tend-tstart)*(ts(j)+1)/2
              ii = ii+1
           enddo
           hcur = hcur/2
        enddo

        n = nlev*k

        ii  =1
        do i=1,n
        do j=1,n

        xtmp(ii) = x(i)
        ytmp(ii) = x(j)
        ii = ii +1

        enddo
        enddo

        return
        end
c----------------------------------------------------------------------------------

