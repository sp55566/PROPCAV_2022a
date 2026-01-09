      SUBROUTINE POSIT_SC1(CAVL1,JS,IDR)
************************************************************************
*                                                                      *
*   Subroutine POSIT determines the exact POSITion of the cavity       *
*   detachment and reattachment points and the relevant indices.       *
*                                                                      *
*   Author: Neal Fine  2-7-91                                          *
*   Date of last Revision                     Revision                 *
*   ---------------------                   ------------               *
*   JY092601                Copied from posit.f.  Modified for ISC=1.  *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      DATA FTOL /0.001/
     
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

C....Cavity length at blade T.E.
      CAVLTE=ARCLNG(NHOLD+1,JS,IDR)/ARCLNG(NHP,JS,IDR)
      DUMM=ARCLNG(NLEP(JS,IDXREV,IDR)+1,JS,IDR)/ARCLNG(NHP,JS,IDR)

C-----------------------------------------------------------------------
C     determine the number of cavitating panels and the characteristics
C     of the split panel, if one exists
C-----------------------------------------------------------------------
      IF(NOCAV(JS,IDXREV,IDR).EQ.1) THEN
         CAVL(JS,IDR)=ZERO
         JCV(JS,IDR)=0
         LCV(JS,IDR)=M0(JS,IDR)
         NSPP(JS,IDR)=0
      ELSE IF(CAVL1.LT.CAVLTE)THEN
C.......negative cavity lengths are zeroed..............................
         IF(CAVL1.LE.DUMM)THEN
            MCAV=0
            GOTO 40
         ENDIF
C.......search for the nearest panel boundary to the right of the ......
C.......cavity t.e......................................................
         DO 30 I=1,NHOLD

            IF(IDR.EQ.1) THEN
               J=NHP+I
            ELSE
               J=NH-I
            END IF

C..........If the wetted space between cavity and SR is less than 4 
C..........panels, then let those cavitate too.
            IF(I.GE.NHOLD-3.AND.I.GE.NLEP(JS,IDXREV,IDR))THEN
               CAVL(JS,IDR)=CAVLTE

               IF(IDR.EQ.2) THEN
                  MCAV=M0(JS,IDR)-N0(IDR)
                  LCAV=N0(IDR)-1
               ELSE IF(IDR.EQ.1) THEN
                  MCAV=N0(IDR)-M0(JS,IDR)
                  LCAV=N0(IDR)
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
C                  CAVL(JS,IDR)=ARCLNG(NLEP(JS,IDXREV,IDR)+1,JS,IDR)/
C     *                 ARCLNG(NHP,JS,IDR)
                  CAVL(JS,IDR)=ZERO
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
                  IF(IDR.EQ.2) THEN
                     LCAV=J+1
                     MCAV=M0(JS,IDR)-2-NH+I
                  ELSE IF(IDR.EQ.1) THEN
                     LCAV=J-1
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

         IF(CAVL1.LE.DUMM.OR.MCAV.EQ.0)THEN
C            NOCAV(JS,IDXREV,IDR)=1
            CAVL(JS,IDR)=ZERO
            JCV(JS,IDR)=0
            LCV(JS,IDR)=M0(JS,IDR)
            NSPP(JS,IDR)=0
         ELSE
            JCV(JS,IDR)=MCAV
            LCV(JS,IDR)=LCAV
         ENDIF

      ELSE 

         IF(IDR.EQ.2)THEN
            JCV(JS,IDR)=M0(JS,IDR)-N0(IDR)
            LCV(JS,IDR)=N0(IDR)-1
         ELSE
            JCV(JS,IDR)=N0(IDR)-M0(JS,IDR)
            LCV(JS,IDR)=N0(IDR)
         ENDIF
         IF(NLEP(JS,IDXREV,IDR).LT.NHOLD) SOP(JS,IDR)=ONE
         NSPP(JS,IDR)=0
         CAVL(JS,IDR)=CAVLTE

      END IF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine posit>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END

