      SUBROUTINE QWEXTRAP(M,IDR)
************************************************************************
*   use a mixed polynomial/sqrt singular extrapolation of the form     *
*   QW(s)=A*S*S + B*S + C/SQRT(SL-S) (all three terms if nwc>3, the    *
*   last two terms if nwc=2, and only the last term if nwc=1           *
*   moved here from delr.f on 06-23-92 NF                              *
*   Date of last revision     Revision                                 *
*   ---------------------     --------                                 *
*   07-28-92 NF      rewrote to include changes to qtextrap           *
*   02-06-99 JY      modified subroutine to allow cavity to grow on    *
*                    both the back and face of the foil.               *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      N1=NC+NWC(M,IDR)
      NSC=NWC(M,IDR)

      RMAX=0.0

C      IF(NWC(M,IDR).GE.4)THEN
C.......3-point extrapolation
C        S3=SPZ(N1-2,M)-SPZ(N1-3,M)
C        S2=SPZ(N1-1,M)-SPZ(N1-3,M)
C        S1=SPZ(N1,M)  -SPZ(N1-3,M)
C        S0=SPZ(N1+1,M)-SPZ(N1-3,M)
C        SL=S0+HALF*FLS(M,IDR)*DZW(NSC+1,M)
C        IF(RZP(M).LE.RMAX)THEN
C          CALL QEXT1(S0,S1,S2,S3,SL,QW(M,1,IDR),QW(M,2,IDR),
C     *          QW(M,3,IDR),QW(M,4,IDR))
C        ELSE
C          CALL CUBEXT(S0,S1,S2,S3,QW(M,1,IDR),QW(M,2,IDR),
C     *          QW(M,3,IDR),QW(M,4,IDR))
C        ENDIF
C      ELSEIF(NWC(M,IDR).EQ.3)THEN
C.......3-point extrapolation
C        S2=SPZ(N1-1,M)-SPZ(N1-2,M)
C        S1=SPZ(N1,M)  -SPZ(N1-2,M)
C        S0=SPZ(N1+1,M)-SPZ(N1-2,M)
C        SL=S0+HALF*FLS(M,IDR)*DZW(NSC+1,M)
C        IF(RZP(M).LE.RMAX)THEN
C          CALL QEXT2(S0,S1,S2,SL,QW(M,1,IDR),QW(M,2,IDR),QW(M,3,IDR))
C        ELSE
C          CALL QUADEXT(S0,S1,S2,QW(M,1,IDR),QW(M,2,IDR),QW(M,3,IDR))
C        ENDIF
C        QW(M,4,IDR)=ZERO
C      ELSEIF(NWC(M,IDR).EQ.2)THEN
      IF(NWC(M,IDR).GE.2)THEN
C.......2-point extrapolation
        S0=SPZ(N1+1,M)-SPZ(N1-1,M)
        S1=SPZ(N1,M)  -SPZ(N1-1,M)
        SL=S0+HALF*FLS(M,IDR)*DZW(NSC+1,M)
        IF(HRZP(1,M).LE.RMAX)THEN
          CALL QEXT3(S0,S1,SL,QW(M,1,IDR),QW(M,2,IDR))
        ELSE
          QW(M,1,IDR)=S0/S1
          QW(M,2,IDR)=ONE-QW(M,1,IDR)
        ENDIF
        QW(M,3,IDR)=ZERO
        QW(M,4,IDR)=ZERO
      ELSE
C.......1-point extrapolation
        QW(M,1,IDR)=ONE
        QW(M,2,IDR)=ZERO
        QW(M,3,IDR)=ZERO
        QW(M,4,IDR)=ZERO
      ENDIF


      RETURN
      END
