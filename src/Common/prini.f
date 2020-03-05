c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning of the
c       printing subroutines proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c       This file contains 18 user-callable subroutines: prini, prin,
c       prin2, prin2_long, prinq, prinq_long, prinf, prinf_long, prinf2,
c       prina, prinl1, prini_flush, prini_off, prini_on, msgmerge,
c       fileflush, mach_zero, and length.  The following is a brief
c       description of these subroutines.
c
c   prini - Initialize the printing routines with two I/O units
c
c   prin - Print a real*4 array, 6 numbers per line
c
c   prin2 - Print a real*8 array, 6 numbers per line
c
c   prin2_long - Print a real*8 array, 2 numbers per line
c
c   prinq - Print a real*16 array, 6 numbers per line
c
c   prinq_long - Print a real*16 array, 2 numbers per line
c
c   prinf - Print an integer*4 array, 10 numbers per line
c
c   prinf_long - Print an integer*4 array, 6 numbers per line
c
c   prinf2 - Print an integer*2 array, 10 numbers per line
c
c   prina - Print a string
c
c   prinl1 - Print a logical*1 array, 20 logical values per line
c
c   prini_flush - Flush the output buffers for the I/O units provided to prini
c       
c   prini_off - Turn the printing off
c
c   prini_on - Turn the printing back on
c
c
c       This file also provides the following utility subroutines.
c
c   msgmerge - Merge two strings, omitting trailing blanks
c
c   fileflush - Flush the output buffer for the specified I/O unit
c
c   mach_zero - Compute the machine epsilon for real*8 calculations
c
c   length - Compute the length of a string minus trailing blanks
c 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        SUBROUTINE PRINI(IP1,IQ1)
        save
        CHARACTER MES*(*), AA(1)
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *8 A4(1)
ccc        REAL *16 A4(1)
        INTEGER *4 IA(1)
        INTEGER *2 IA2(1)
        logical *1 la(1)
        integer ison
        data ison/1/

        IP=IP1
        IQ=IQ1
  
        RETURN
  
C 
C 
C 
C 
C 
        ENTRY PRIN(MES,A,N)
        if (ison .le. 0) return
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
        if (ison .le. 0) return
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRIN2_long(MES,A2,N)
        if (ison .le. 0) return
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1450)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1450)(A2(J),J=1,N)
 1450 FORMAT(2(2X,E22.16))
        RETURN
C 
C 
C 
C 
        ENTRY PRINQ(MES,A4,N)
        if (ison .le. 0) return
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRINQ_long(MES,A4,N)
        if (ison .le. 0) return
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1550)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1550)(A4(J),J=1,N)
 1550 FORMAT(2(2X,E22.16))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF(MES,IA,N)
        if (ison .le. 0) return
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF_long(MES,IA,N)
        if (ison .le. 0) return
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1700)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1700)(IA(J),J=1,N)
 1700 FORMAT(6(2X,I11))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF2(MES,IA2,N)
        if (ison .le. 0) return
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C 
C 
C 
C 
        ENTRY PRINA(MES,AA,N)
        if (ison .le. 0) return
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
C 
C 
C 
C 
        entry prinl1(mes,la,n)
        if (ison .le. 0) return
        CALL MESSPR(MES,ip,iq)
 2200 format(20L3)
c
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2200)(lA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2200)(lA(J),J=1,N)
c
        return
c
c
c
c
        entry prini_flush()

        if (ip.ne.0 .and. ip.ne.6) call fileflush(ip)
        if (iq.ne.0 .and. iq.ne.6) call fileflush(iq)

        return
c
c
c
c
        entry prini_off()
        ison=0

        return
c
c
c
c
        entry prini_on()
        ison=1

        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE MESSPR(MES,IP,IQ)
        save
        CHARACTER mes*(*), ast*1
        data ast/'*'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        i1=length(mes)
        IF(MES(i1:i1).EQ. ast) i1=i1-1

         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) mes(1:i1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) mes(1:i1)
 1800 FORMAT(1X,A)
         RETURN
         END
C 
C 
C 
C 
C 
        subroutine msgmerge(a,b,c)
        save
        character a*(*),b*(*),c*(*),ast
        data ast/'*'/
c 
        la=length(a)
        if(a(la:la) .eq. ast) la=la-1

        c=a(1:la)

        lb=length(b)
        if(b(lb:lb) .eq. ast) lb=lb-1

        c(la+1:)=b(1:lb)

        return
        end
c 
c 
c 
c 
c 
  
        subroutine fileflush(iw)
        implicit real *8 (a-h,o-z)
c
ccc        call prinf('flushing unit',iw,1)
c 
        save
        close(iw)
        open(iw,status='old')
        do 1400 i=1,1000000
c 
        read(iw,1200,end=1600)
 1200 format(1a1)
 1400 continue
 1600 continue
c
        backspace iw
c 
        return
        end
  
c 
c 
c 
c 
c  
        subroutine mach_zero(zero_mach7)
        implicit real *8 (a-h,o-z)
        data isinit/-7/
        save

        zero_mach7=zero_mach
        if (isinit .eq. 1) return
c
c
ccc        call prinf('computing zero_mach...',0,0)

        zero_mach=100       
c
        d1=1.1d0
        d=1.11d0
        do 1200 i=1,10000
c

        d=d/2
        d2=d1+d
        d4=sin(sin(d1)-sin(d2))
c
        if(d4 .eq. 0) goto 1400
c
 1200 continue
 1400 continue
c
        zero_mach=d
        isinit=1
c
c
        zero_mach7=zero_mach

        return
        end
c
c
c
c
c
      INTEGER FUNCTION LENGTH(STRING)
*Returns length of string ignoring trailing blanks
      CHARACTER*(*) STRING
      DO 15, I = LEN(STRING), 1, -1
      IF(STRING(I:I) .NE. ' ') GO TO 20
15    CONTINUE
20    LENGTH = I
      END
        
