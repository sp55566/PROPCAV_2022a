      SUBROUTINE POSIT(CAVL1,JS,IDR)
************************************************************************
*                                                                      *
*   Subroutine POSIT determines the exact POSITion of the cavity       *
*   detachment and reattachment points and the relevant indices.       *
*                                                                      *
*   Author: Neal Fine  2-7-91                                          *
*   Date of last Revision                     Revision                 *
*   ---------------------                   ------------               *
*     05-07-91                  -included lamda and relevant indices   *
*     07-24-91                  -find the nearest panel to the right   *
*                                of the cavity length....call that     *
*                                panel # lcav                          *
*     03-19-91  NF       -re-written for psfcav                        *
*     07-20-92  NF       -changed FTOL to .0001 from .01               *
*     04-02-97  CM       -allow for a variable number of leading edge  *
*                         pre-detachment panels (NLEP(MR))             *
*     02-06-99  JY       -modified subroutine to allow cavity to grow  *
*                         on both the back and face of the foil.       *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      DATA FTOL /0.001/
      DATA FTOLSC /0.01/
     
C-----------------------------------------------------------------------
C     define the chord length at the midspan of the strip.
C-----------------------------------------------------------------------
      CH=HALF*(CHORD(JS)+CHORD(JS+1))

C-----------------------------------------------------------------------
C     determine the exact position of the cavity leading edge.
C-----------------------------------------------------------------------
      IF(IDR.EQ.2) THEN
         M0(JS,IDR)=NHP-NLEP(JS,IDXREV,IDR)
      ELSE IF(IDR.EQ.1) THEN
         M0(JS,IDR)=NHP+NLEP(JS,IDXREV,IDR)
      ENDIF

C-----------------------------------------------------------------------
C     determine the number of cavitating panels and the characteristics
C     of the split panel, if one exists
C-----------------------------------------------------------------------
      IF(CAVL1.LT.ONE)THEN
         DUMM=ARCLNG(NLEP(JS,IDXREV,IDR)+1,JS,IDR)/ARCLNG(NHP,JS,IDR)

C.......partial cavity (zero the supercavity indices)...................
         NWC(JS,IDR)=0
         NSPS(JS,IDR)=0
         SOP(JS,IDR)=ZERO
C.......negative cavity lengths are zeroed..............................
         IF(CAVL1.LE.DUMM.OR.NOCAV(JS,IDXREV,IDR).EQ.1)THEN
            MCAV=0
            GOTO 40
         ENDIF
C.......search for the nearest panel boundary to the right of the ......
C.......cavity t.e......................................................
         DO 30 I=1,NH
            J=NHP+I
C..........don't let the last panel on the blade be split!..............
            IF(I.EQ.NH)THEN
C..........If the cavity ends at the last blade panel, push it one......
C..........to the right because if not, we'll have stability problem,...
C..........and especially problem with applying IPK (JY082299)..........
C
C/s S.N.KIM... CAVL better be CAVL1. Let it be free to change, not be
C              pushed forward to 1.0+DWZ*HALF/CH. It would cause stupid
C              infinite loop whenever I.EQ.NH.
               CAVL(JS,IDR)=CAVL1 !ONE+DZW(1,JS)*HALF/CH
               NWC(JS,IDR)=1
               LCAV=NCP
               IF(IDR.EQ.2) THEN
                  MCAV=M0(JS,IDR)-1
               ELSE IF(IDR.EQ.1) THEN
                  MCAV=NCP-M0(JS,IDR)
               ENDIF
               NSPP(JS,IDR)=0
               SOP(JS,IDR)=ONE
               GOTO 40
            ENDIF
            XCTE=ARCLNG(I+1,JS,IDR)/ARCLNG(NHP,JS,IDR)
            DSP=XCTE-CAVL1
            IF(DSP.LT.ZERO)THEN
               XCTEKM1=XCTE
               GOTO 30
            ELSEIF(DSP.EQ.ZERO)THEN
               NSPP(JS,IDR)=0
               LCAV=J
               IF(IDR.EQ.2) THEN
                  MCAV=M0(JS,IDR)-1-NH+I
               ELSE IF(IDR.EQ.1) THEN
                  MCAV=LCAV-M0(JS,IDR)
               ENDIF
               CAVL(JS,IDR)=CAVL1
               GOTO 40
            ELSE
C.............Don't allow the first cavitating panel to split(JY060899).
               IF(I.EQ.NLEP(JS,IDXREV,IDR)+1)THEN
                  MCAV=0
                  LCAV=M0(JS,IDR)
                  NSPP(JS,IDR)=0
                  CAVL(JS,IDR)=ARCLNG(NLEP(JS,IDXREV,IDR)+1,JS,IDR)/
     *                 ARCLNG(NHP,JS,IDR)
                  GOTO 40
               ENDIF
               DXP=XCTE-XCTEKM1
               TOL0=DSP/DXP
               TOL1=ONE-TOL0
               IF(TOL0.LT.FTOL)THEN
                  NSPP(JS,IDR)=0
                  LCAV=J
                  IF(IDR.EQ.2) THEN
                     MCAV=M0(JS,IDR)-1-NH+I
                  ELSE IF(IDR.EQ.1) THEN
                     MCAV=LCAV-M0(JS,IDR)
                  ENDIF
                  CAVL(JS,IDR)=XCTE
                  GOTO 40
               ELSEIF(TOL1.LT.FTOL)THEN
                  NSPP(JS,IDR)=0
                  LCAV=J-1
                  IF(IDR.EQ.2) THEN
                     MCAV=M0(JS,IDR)-2-NH+I
                  ELSE IF(IDR.EQ.1) THEN
                     MCAV=LCAV-M0(JS,IDR)
                  ENDIF
                  CAVL(JS,IDR)=XCTEKM1
                  GOTO 40
               ELSE
                  LCAV=J
                  IF(IDR.EQ.2) THEN
                     MCAV=M0(JS,IDR)-2-NH+I
                  ELSE IF(IDR.EQ.1) THEN
                     MCAV=LCAV-M0(JS,IDR)-1
                  ENDIF
                  NSPP(JS,IDR)=1
                  CAVL(JS,IDR)=CAVL1
                  FLP(JS,IDR)=TOL1
                  FRP(JS,IDR)=TOL0
                  IF(IDR.EQ.2) THEN
                     LCAV1=M0(JS,IDR)-1-MCAV
                  ELSE IF(IDR.EQ.1) THEN
                     LCAV1=LCAV-1
                  ENDIF
                  DZL(JS,IDR)=FLP(JS,IDR)*DZ(LCAV1,JS)
                  DZR(JS,IDR)=FRP(JS,IDR)*DZ(LCAV1,JS)
                  GOTO 40
               ENDIF
            ENDIF
 30      CONTINUE
 40      CONTINUE

         IF(NOCAV(JS,IDXREV,IDR).EQ.1) THEN
            CAVL(JS,IDR)=ZERO
            JCV(JS,IDR)=0
            NSPP(JS,IDR)=0
            LCV(JS,IDR)=M0(JS,IDR)
         ELSE IF(CAVL1.LE.DUMM)THEN
            CAVL(JS,IDR)=DUMM
            JCV(JS,IDR)=0
            LCV(JS,IDR)=M0(JS,IDR)
            NSPP(JS,IDR)=0
         ELSE
            JCV(JS,IDR)=MCAV
            LCV(JS,IDR)=LCAV
         ENDIF

      ELSE IF(NOCAV(JS,IDXREV,IDR).EQ.0) THEN

C.......supercavity.....................................................

         IF(IDR.EQ.2)THEN
            JCV(JS,IDR)=M0(JS,IDR)-1
         ELSE
            JCV(JS,IDR)=NCP-M0(JS,IDR)
         ENDIF
         LCV(JS,IDR)=NCP
         SOP(JS,IDR)=ONE
         NSPP(JS,IDR)=0
         XCTE=ONE
         DO 50 I=1,NTRA
            XCTE=XCTE+DZW(I,JS)*HALF/CH
            DSP=XCTE-CAVL1
            IF(I.EQ.NTRA)THEN
               NSPS(JS,IDR)=0
               NWC(JS,IDR)=I
               CAVL(JS,IDR)=XCTE
C..............The next warning is added by JY020699....................
C               WRITE(*,2000) 
C2000          FORMAT('WARNING: Supercavity length exceeds max limit!!')
               GOTO 60
            ENDIF
            IF(DSP.LT.ZERO)THEN
               XCTEKM1=XCTE
               GOTO 50
            ELSEIF(DSP.EQ.ZERO)THEN
               NSPS(JS,IDR)=0
               NWC(JS,IDR)=I
               CAVL(JS,IDR)=CAVL1
               GOTO 60
            ELSE
               NWC(JS,IDR)=I-1
               IF(I.EQ.1)THEN
                  NSPS(JS,IDR)=0
C/s S.N.KIM... CAVL better be CAVL1. Let it be free to change, not be
C              pushed forward to XCTE. It would cause stupid infinite
C              loop whenever I.EQ.1.
                  CAVL(JS,IDR)=CAVL1 !XCTE
                  NWC(JS,IDR)=1
               ELSE
                  DXP=XCTE-XCTEKM1
                  TOL0=DSP/DXP
                  TOL1=ONE-TOL0
                  IF(TOL0.LT.FTOLSC)THEN
                     NWC(JS,IDR)=I
                     NSPS(JS,IDR)=0
                     CAVL(JS,IDR)=XCTE
                  ELSEIF(TOL1.LT.FTOLSC)THEN
                     NSPS(JS,IDR)=0
                     CAVL(JS,IDR)=XCTEKM1
                  ELSE
                     NSPS(JS,IDR)=1
                     CAVL(JS,IDR)=CAVL1
                     FLS(JS,IDR)=TOL1
                     FRS(JS,IDR)=TOL0
                     DZL(JS,IDR)=FLS(JS,IDR)*DZW(I,JS)
                     DZR(JS,IDR)=FRS(JS,IDR)*DZW(I,JS)
                  ENDIF
               ENDIF
               GOTO 60
            ENDIF
 50      CONTINUE
 60      CONTINUE
      ENDIF
 99   CONTINUE

C-----------------------------------------------------------------------
C     determine the number of cavitating panels within which the 
C     pressure model is applied **CURRENTLY PARTIAL CAVITATION ONLY**
C     
C     Note: The pressure model is not applied even for the partial 
C           cavitation case if LAMBDA=0.0 in *.geo.            JY020699
C-----------------------------------------------------------------------
      RLAMDA=RLAMDA1*CAVL(JS,IDR)
      IF(RLAMDA.EQ.ZERO)THEN
         RLAM(JS,IDR)=ZERO
      ELSE
         RLAM(JS,IDR)=RLAMDA1
      ENDIF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine posit>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END

