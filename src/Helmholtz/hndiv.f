      subroutine hndiv(eps,ns,nt,ifcharge,ifdipole,ifpgh,
     1   ifpghtarg,ndiv,idivflag)
c
c
c       this subroutine estimates ndiv and idivflag 
c       based on geometry parameters
c     
c       TODO: better way to choose ndiv based on boxsize
c
c
      implicit none
      real *8 eps
      integer(8) ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg,ndiv
      integer(8) idivflag

      idivflag = 0

       if(eps.ge.0.5d-0) then
         ndiv = 40
       else if(eps.ge.0.5d-1) then
         ndiv = 40
       else if(eps.ge.0.5d-2) then
         ndiv = 40
       else if(eps.ge.0.5d-3) then
         ndiv = 80
       else if(eps.ge.0.5d-6) then
         ndiv = 200
       else if(eps.ge.0.5d-9) then
         ndiv = 400
       else if(eps.ge.0.5d-12) then
         ndiv = 600
       else if(eps.ge.0.5d-15) then
         ndiv = 700
       else
         ndiv = ns+nt
       endif

      return
      end
