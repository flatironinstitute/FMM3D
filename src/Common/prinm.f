ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Multipole printing routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
      SUBROUTINE PRINM(MPOLE,NTERMS)
      implicit real *8 (a-h,o-z)
      real *8 MPOLE0(0:NTERMS,0:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,-nterms:NTERMS)
      INTEGER NTERMS
C
C     print out coefficients of multipole expansion
C
1000  FORMAT(6(2X,E11.5))
1001  FORMAT(/)
      DO 100 L = 0,NTERMS
         WRITE(6,1000)(MPOLE(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
100   CONTINUE
        return
C
C
C
C
      ENTRY PRINM0(MPOLE0,NTERMS)
      DO 200 L = 0,NTERMS
         WRITE(6,1000)(MPOLE0(L,M),M=0,L)
         WRITE(13,1000)(MPOLE0(L,M),M=0,L)
         WRITE(6,1001)
         WRITE(13,1001)
200   CONTINUE
c
      RETURN
c
c
c
      ENTRY PRINM1(MPOLE,NTERMS,NTERMS2)
      DO 300 L = 0,NTERMS2
         WRITE(6,1000)(MPOLE(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
300   CONTINUE
      RETURN
c
      ENTRY PRINMBORDER(MPOLE,NTERMS)
      DO 400 L = 0,NTERMS-1
         WRITE(6,1000)MPOLE(L,L), MPOLE(L+1,L)
         WRITE(13,1000)MPOLE(L,L), MPOLE(L+1,L)
         WRITE(6,1001)
         WRITE(13,1001)
400   CONTINUE
c
      RETURN
        end
c
c
