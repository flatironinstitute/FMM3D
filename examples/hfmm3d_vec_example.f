      implicit none
      integer ns,nt,nd
      double precision, allocatable :: source(:,:),targ(:,:)
      double complex, allocatable :: charge(:,:),dipvec(:,:,:)
      double complex, allocatable :: pot(:,:)
      double complex, allocatable :: pottarg(:,:)

      double precision eps
      double complex eye,zk
      integer i,j,k,idim,ier
      double precision hkrand,pi,thet,phi
      

      data eye/(0.0d0,1.0d0)/

      pi = 4.0d0*atan(1.0d0)

c
cc      initialize printing routine
c
      call prini(6,13)
      write(*,*)
      write(*,*)
      write(*,*) "================================="
      call prin2("This code is an example fortran driver*",i,0)
      call prin2("On output, the code prints sample pot,pottarg*",i,0)
      write(*,*)
      write(*,*)

      zk = 2.2d0

      ns = 10000
      nt = 10000
      
      nd = 1

      allocate(source(3,ns))
      allocate(targ(3,nt))
      allocate(charge(nd,ns))
      allocate(dipvec(nd,3,ns))
      allocate(pot(nd,ns))
      allocate(pottarg(nd,nt))


      eps = 0.51d-3

      write(*,*) "=========================================="

c
c   
c       example demonstrating use of 
c        source to source+targ, charges+dipoles, pot
c



c
cc      generate sources uniformly on the sphere 
c
c
      do i=1,ns
        thet = hkrand(0)*pi
        phi = hkrand(0)*2*pi
        source(1,i) = sin(thet)*cos(phi) 
        source(2,i) = sin(thet)*sin(phi)
        source(3,i) = cos(thet)

        do idim=1,nd
          charge(idim,i) = hkrand(0) + eye*hkrand(0)
          dipvec(idim,1,i) = hkrand(0) + eye*hkrand(0)
          dipvec(idim,2,i) = hkrand(0) + eye*hkrand(0)
          dipvec(idim,3,i) = hkrand(0) + eye*hkrand(0)
          pot(idim,i) = 0
        enddo
      enddo

      do i=1,nt
        targ(1,i) = source(1,i) + 10.0d0
        targ(2,i) = source(2,i)
        targ(3,i) = source(3,i)

        do idim=1,nd
          pottarg(idim,i) = 0
        enddo
      enddo


       call hfmm3d_st_cd_p_vec(nd,eps,zk,ns,source,charge,
     1      dipvec,pot,nt,targ,pottarg,ier)

       call prin2("potential at sources=*",pot,12)
       call prin2("potential at targets=*",pottarg,12)


      stop
      end
c----------------------------------------------------------
c
cc
c
c
