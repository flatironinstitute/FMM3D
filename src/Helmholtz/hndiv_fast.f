      subroutine hndiv(dbsize,eps,ns,nt,ifcharge,ifdipole,ifpgh,
     1   ifpghtarg,ndiv,idivflag)
c
c
c       this subroutine estimates ndiv and idivflag 
c       based on geometry parameters
c       
c
      implicit none
      real *8 eps,dbsize,rfac
      integer *8 ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg,ndiv
      integer *8 idivflag

      idivflag = 0
      rfac = 1.0d0
      if (dbsize.gt.8) rfac = (dbsize/8)**1.5d0
      

       if(eps.ge.0.5d-0) then
         ndiv = ceiling(300*rfac)
       else if(eps.ge.0.5d-1) then
         ndiv = ceiling(300*rfac)
       else if(eps.ge.0.5d-2) then
         ndiv = ceiling(300*rfac)
       else if(eps.ge.0.5d-3) then
         ndiv = ceiling(300*rfac)
       else if(eps.ge.0.5d-6) then
         ndiv = ceiling(1000*rfac)
       else if(eps.ge.0.5d-9) then
         ndiv = ceiling(1000*rfac)
       else if(eps.ge.0.5d-12) then
         ndiv = ceiling(1000*rfac)
       else if(eps.ge.0.5d-15) then
         ndiv = ceiling(1000*rfac)
       else
         ndiv = ns+nt
       endif

       ndiv = min(ndiv, 10000)

      return
      end
