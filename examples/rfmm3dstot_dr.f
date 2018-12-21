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
        double precision charge(nsmax)
        double precision dipstr(nsmax)
        dimension dipvec(3,nsmax)
        
c
        dimension targ(3,ntmax)
        double precision pottarg(ntmax)
        double precision fldtarg(3,ntmax)

        double precision pottargex(ntest),fldtargex(3,ntest)

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
cc         distributions: idist = 1, sources and targs in the 
c                                    volume of a cube
c                         idist = 2, sources and targs on the 
c                                    surface of a sphere
c                         idist = 3, sources on surface of a sphere
c                                    targs in the volume
c                                    (number of targs = nsources^1.5)
c                         idist = 4, sources on surface of a cube
c                                    with a dyadic mesh
c  


        call prin2('Enter n*',i,0)
        read *, n

        icase = 1
        idist= 2
        iprec = 2

        ifcharge = 1
        ifdipole = 0


        ifpottarg = 1
        iffldtarg = 1


c
c
        if( idist .eq. 1 ) then
c
c       ... construct randomly located charge distribution on a unit 
c            cube
c
           nsource = icase*10000

           nsource = 200000
           ntarg = nsource
           do i=1,nsource

              source(1,i)=hkrand(0)-0.5d0
              source(2,i)=hkrand(0)-0.5d0
              source(3,i)=hkrand(0)-0.5d0

              targ(1,i) = hkrand(0)-0.5d0
              targ(2,i) = hkrand(0)-0.5d0
              targ(3,i) = hkrand(0)-0.5d0

              ii = i +nsource
              targ(1,ii) = hkrand(0)-0.5d0
              targ(2,ii) = hkrand(0)-0.5d0
              targ(3,ii) = hkrand(0)-0.5d0
           enddo
        endif
c
        if( idist .eq. 2 ) then
c
c       ... construct randomly located charge distribution on 
c             surface of sphere
c

           nsource = icase*100000

           ntarg = nsource
           r = 0.5d0
           do i=1,nsource
              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              source(1,i) = r*cos(thet)*cos(phi)
              source(2,i) = r*cos(thet)*sin(phi)
              source(3,i) = r*sin(thet)

              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              targ(1,i) = r*cos(thet)*cos(phi)
              targ(2,i) = r*cos(thet)*sin(phi)
              targ(3,i) = r*sin(thet)

              ii = i + nsource
              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              targ(1,ii) = r*cos(thet)*cos(phi)
              targ(2,ii) = r*cos(thet)*sin(phi)
              targ(3,ii) = r*sin(thet)


           enddo
        endif
c
        if( idist .eq. 3 ) then
c
c       ... construct randomly located charge distribution on 
c            a unit sphere and targs in the volume
c
          
           nsource = icase*10000 
           ntarg = nsource**(1.5d0)/5
           r = 0.5d0
           do i=1,nsource
              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              source(1,i) = r*cos(thet)*cos(phi)
              source(2,i) = r*cos(thet)*sin(phi)
              source(3,i) = r*sin(thet)

           enddo

           do i=1,ntarg
              thet = hkrand(0)*pi
              phi = hkrand(0)*2*pi

              rr=  hkrand(0)*r

              targ(1,i) = rr*cos(thet)*cos(phi)
              targ(2,i) = rr*cos(thet)*sin(phi)
              targ(3,i) = rr*sin(thet)
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

           ntarg = nsource

           do i=1,nsource
              source(1,i) = a1*xtmp(i) + b1*ytmp(i) 
              source(2,i) = a2*xtmp(i) + b2*ytmp(i)
              source(3,i) = a3*xtmp(i) + b3*ytmp(i)

              targ(1,i) = source(1,i)
              targ(2,i) = source(2,i)
              targ(3,i) = source(3,i)
           enddo

        endif


        do i=1,nsource
           charge(i) = 1
           dipstr(i) = 1
           dipvec(1,i) = hkrand(0)
           dipvec(2,i) = hkrand(0)
           dipvec(3,i) = hkrand(0)
        enddo

        t1 = second()


        do i=1,ntarg
           pottarg(i) =0
           fldtarg(1,i) = 0
           fldtarg(2,i) = 0
           fldtarg(3,i) = 0
        enddo


        call rfmm3dpartstot(ier,iprec,nsource,source,ifcharge,
     1          charge,ifdipole,dipstr,dipvec,
     2          ntarg,targ,ifpottarg,pottarg,
     3          iffldtarg,fldtarg)


        m=min(nsource,100)


        call geterr(nsource,source,ifcharge,charge,ifdipole,dipstr,
     1               dipvec,m,targ,ifpottarg,pottarg,iffldtarg,
     2               fldtarg,rerr,rgerr)
        if(ifpottarg.eq.1) call prin2('targ pot error=*',rerr,1)
        if(iffldtarg.eq.1) call prin2('targ field error=*',rgerr,1)

        
        return
        end
     
c-----------------------------------------------------------------

        subroutine geterr(nsource,source,ifcharge,charge,ifdipole,
     1               dipstr,dipvec,m,targ,ifpottarg,pottarg,
     2               iffldtarg,fldtarg,rerr,rgerr)

        implicit double precision (a-h,o-z)
        double precision source(3,*),targ(3,*),dipvec(3,*)
        double precision charge(*),dipstr(*),pottarg(*),fldtarg(3,*)
        double precision, allocatable :: pot2(:),fld2(:,:)
        double complex cpot2(m),cfld2(3,m),ptemp,ftemp(3)
        double complex chrg,dipl

        rerr = 0
        rgerr = 0

        allocate(pot2(m),fld2(3,m))


        do i=1,m
        if (ifpottarg .eq. 1) cpot2(i)=0
        if (iffldtarg .eq. 1) then
           cfld2(1,i)=0
           cfld2(2,i)=0
           cfld2(3,i)=0
        endif
        enddo
c
        t1=second()
        do 8160 j=1,m
        do 8150 i=1,nsource


        rr = (source(1,i)-targ(1,j))**2 + 
     2          (source(2,i)-targ(2,j))**2 + 
     3          (source(3,i)-targ(3,j))**2

        if(rr.gt.1.0d-28) then

        if( ifcharge .eq. 1 ) then

        chrg = charge(i)
        
        call lpotfld3d(iffldtarg,source(1,i),chrg,
     $     targ(1,j),
     1     ptemp,ftemp)
        if (ifpottarg .eq. 1) cpot2(j)=cpot2(j)+ptemp
        if (iffldtarg .eq. 1) then
           cfld2(1,j)=cfld2(1,j)+ftemp(1)
           cfld2(2,j)=cfld2(2,j)+ftemp(2)
           cfld2(3,j)=cfld2(3,j)+ftemp(3)
        endif
        endif
        if (ifdipole .eq. 1) then
           dipl = dipstr(i)
           call lpotfld3d_dp(iffldtarg,source(1,i),
     $     dipl,dipvec(1,i),
     $     targ(1,j),ptemp,ftemp)
           if (ifpottarg .eq. 1) cpot2(j)=cpot2(j)+ptemp
           if (iffldtarg .eq. 1) then
              cfld2(1,j)=cfld2(1,j)+ftemp(1)
              cfld2(2,j)=cfld2(2,j)+ftemp(2)
              cfld2(3,j)=cfld2(3,j)+ftemp(3)
           endif
        endif

        endif
c
 8150   continue
 8160   continue
c
        t2=second()

        if(ifpottarg.eq.1) then
           do i=1,m
              pot2(i) = real(cpot2(i))
           enddo
        endif

        if(iffldtarg.eq.1) then
            do i=1,m
               fld2(1,i) = real(cfld2(1,i))
               fld2(2,i) = real(cfld2(2,i))
               fld2(3,i) = real(cfld2(3,i))
            enddo
        endif

        
C
        rerr = 0
        if(ifpottarg.eq.1) call l3derror(pottarg,pot2,m,aerr,rerr)
     
        rgerr = 0
        if(iffldtarg.eq.1) call l3derror(fldtarg,fld2,3*m,aerr,rgerr)

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
        double precision pot1(n),pot2(n)
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

