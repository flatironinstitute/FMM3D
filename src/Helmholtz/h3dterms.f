cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2010-07-10 23:01:41 -0400 (Sat, 10 Jul 2010) $
c    $Revision: 1082 $
c
c
c-----------------------------------------------------------------------------
c
c      h3dterms - determine number of terms in mpole expansions for box
c           of size "size" with Helmholtz parameter zk.
c
c      h3dterms_list2 - build the number of terms table for all boxes 
c           in list 2
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine h3dterms(size, zk, eps, nterms)
      implicit real *8 (a-h,o-z)
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Maximum number of terms is 1000, which 
c     works for boxes up to 160 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:2000), fhder(0:1)
c
c
      z1 = (zk*size)*1.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
c
      ntmax = 1000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*size) .lt. 1.0d0) rscale = cdabs(zk*size)
      call h3dall(ntmax,z1,rscale,hfun,ifder,fhder)
      z2 = (zk*size) * dsqrt(3d0)/2.0d0

      call besseljs3d(ntmax, z2, rscale, jfun, ifder, fjder)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
c
      enddo
c
c       ... computational box is too big, set nterms to 1000
c
        nterms=1000
c
      return
      end
c
c
c
c
c
      subroutine h3dterms_far(size, zk, eps, nterms)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk. 
c
c     The method is based on examining the decay of h_n * j_n.
c
c     This routine assumes slightly larger separation of boxes: the
c     first unit box is located at the origin and the second box is
c     located at (3,0,0).
c
c     Maximum number of terms is 1000, which 
c     works for boxes up to 160 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:2000), fhder(0:1)
c
      z1 = (zk*size)*2.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
ccc        if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 1000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*size) .lt. 1.0d0) rscale = cdabs(zk*size)
      call h3dall(ntmax,z1,rscale,hfun,ifder,fhder)
      z2 = (zk*size) * dsqrt(3d0)/2.d0
c
      ier1 = 0
c
      call besseljs3d(ntmax, z2, rscale, jfun, ifder, fjder)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
      enddo
c
c       ... computational box is too big, set nterms to 1000
c
        nterms=1000
c
      return
      end
c
c
c
c
