      SUBROUTINE POSIT_SC2(CAVL2,JS)
************************************************************************
*                                                                      *
*   Subroutine POSIT determines the exact POSITion of the cavity       *
*   detachment and reattachment points and the relevant indices.       *
*                                                                      *
*   Author: Neal Fine  2-7-91                                          *
*   Date of last Revision                     Revision                 *
*   ---------------------                   ------------               *
*   JY110901                Copied from posit.f.  Modified for ISC=1.  *
*                           Please note that for ISC=1, DZL & DZR is   *
*                           only defined for split panels on the blade.*
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      DATA FTOLSC /0.01/

C-----------------------------------------------------------------------
C     define the chord length at the midspan of the strip.
C-----------------------------------------------------------------------
      CH=HALF*(CHORD(JS)+CHORD(JS+1))

      IDR=1
      XCTE=ONE
      DO 50 I=1,NTRA
         XCTE=XCTE+DZW(I,JS)*HALF/CH
         DSP=XCTE-CAVL2
         IF(I.EQ.NTRA)THEN
            NSPS(JS,IDR)=0
            NWC(JS,IDR)=I
            CAVLT(JS)=XCTE
C...........The next warning is added by JY020699....................
C            WRITE(*,2000) 
C2000       FORMAT('WARNING: Supercavity length exceeds max limit!!')
            GOTO 60
         ENDIF
         IF(DSP.LT.ZERO)THEN
            XCTEKM1=XCTE
            GOTO 50
         ELSEIF(DSP.EQ.ZERO)THEN
            NSPS(JS,IDR)=0
            NWC(JS,IDR)=I
            CAVLT(JS)=CAVL2
            GOTO 60
         ELSE
            NWC(JS,IDR)=I-1

            IF(I.EQ.1)THEN
               NSPS(JS,IDR)=0
               DXP=DZW(I,JS)
               TOL0=DSP/DXP
               IF(TOL0.GT.0.5) THEN
                  CAVLT(JS)=ONE
               ELSE
                  CAVLT(JS)=XCTE
                  NWC(JS,IDR)=1
               END IF
            ELSE
               DXP=XCTE-XCTEKM1
               TOL0=DSP/DXP
               TOL1=ONE-TOL0
               IF(TOL0.LT.FTOLSC)THEN
                  NWC(JS,IDR)=I
                  NSPS(JS,IDR)=0
                  CAVLT(JS)=XCTE
               ELSEIF(TOL1.LT.FTOLSC)THEN
                  NSPS(JS,IDR)=0
                  CAVLT(JS)=XCTEKM1
               ELSE
                  NSPS(JS,IDR)=1
                  CAVLT(JS)=CAVL2
                  FLS(JS,IDR)=TOL1
                  FRS(JS,IDR)=TOL0
               ENDIF
            ENDIF
            GOTO 60
         ENDIF
 50   CONTINUE
 60   CONTINUE

C....Both sides of the supercavity should equal.
      NWC(JS,2)=NWC(JS,1)
      NNWC(JS)=NWC(JS,1)
      NWDIR(JS)=1
      IF(NSPS(JS,1).GT.0) THEN
         NSPS(JS,2)=NSPS(JS,1)
         FRS(JS,2)=FRS(JS,1)
         FLS(JS,2)=FLS(JS,1)
      END IF

      RETURN
      END
