c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Printing routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C
        SUBROUTINE PRINI(IP1,IQ1)
        CHARACTER *1 MES(1), AA(1)
         save
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *8 A4(1)
ccc        INTEGER *4 IA(1)
        INTEGER IA(1)
        INTEGER *4 IA1(1)
        INTEGER *2 IA2(1)
        data IP/0/,IQ/0/
        IP=IP1
        IQ=IQ1
        RETURN

C
C
C
C
C
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C
C
C
C
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C
C
C
C
        ENTRY PRINQ(MES,A4,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C
C
C
C
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C
C
C
C
        ENTRY PRINF1(MES,IA1,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA1(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA1(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c
c
c
c
c
        SUBROUTINE MESSPR(MES,IP,IQ)
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
c
c
c
c
c
        subroutine msgmerge(a,b,c)
        character *1 a(*),b(*),c(*),ast
        data ast/'*'/
c
        do 1200 i=1,1000
c
        if(a(i) .eq. ast) goto 1400
        c(i)=a(i)       
        iadd=i
 1200 continue
c
 1400 continue
c
        do 1800 i=1,1000
c
        c(iadd+i)=b(i)
        if(b(i) .eq. ast) return
 1800 continue
        return
        end
c
c
c
c
c
c
C
C
      SUBROUTINE PRINM(MPOLE,NTERMS)
      implicit real *8 (a-h,o-z)
      real *8 MPOLE0(0:NTERMS,-nterms:NTERMS)
      real *8 MPOLE2(0:NTERMS,0:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,-nterms:NTERMS)
      INTEGER NTERMS
C
C     print out coefficients of multipole expansion
C
1000  FORMAT(6E12.5)
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
      ENTRY PRINM_TRUNC(MPOLE,NTERMS,NP)
      DO L = 0,NP
         WRITE(6,1000)(MPOLE(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
      ENDDO
      RETURN
C
C
C
C
      ENTRY PRINM0(MPOLE0,NTERMS)
      DO 200 L = 0,NTERMS
         WRITE(6,1000)(MPOLE0(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE0(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
200   CONTINUE
      RETURN
C
C
      ENTRY PRINM2(MPOLE2,NTERMS)
      DO L = 0,NTERMS
         WRITE(6,1000)(MPOLE2(L,M),M=0,L)
         WRITE(13,1000)(MPOLE2(L,M),M=0,L)
         WRITE(6,1001)
         WRITE(13,1001)
      ENDDO
c
c
      RETURN
      end

