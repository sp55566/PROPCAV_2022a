      SUBROUTINE SURFACE
************************************************************************
*     Subroutine created to keep track of panels that are cut by the   *
*     free-surface.                                                    *
*                                                                      *
*     Date              Comments                                       *
*     ----------------  --------------------------------------------   *
*     JY061999          subroutine created.                            *
************************************************************************
      
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'
       DIMENSION XN(2),YN(2),ZN(2),MS(2)
       DIMENSION YB1(NBPZ,MBPZ),ZB1(NBPZ,MBPZ),YH1(NHPZ,MHPZ),
     *      ZH1(NHPZ,MHPZ),NEL(NSTEP)

C-----------------------------------------------------------------------
C     Loop over all time steps
C-----------------------------------------------------------------------
       DO 10 KK=1,NTPREV

          DT1=-DELTAT*(KK-1)

C-----------------------------------------------------------------------
C     Check for submerged panels on the key blade
C-----------------------------------------------------------------------
          DO 20 M=1,MR+1
             DO 30 N=1,NC+1
                YB1(N,M)=YB(N,M)*COS(DT1)-ZB(N,M)*SIN(DT1)
                ZB1(N,M)=YB(N,M)*SIN(DT1)+ZB(N,M)*COS(DT1)
 30          CONTINUE
 20       CONTINUE
          
          NEL(KK)=0
          
          DO 40 IDR=1,2
 
             DO 50 M=1,MR
                ICB(M,IDR,KK)=0
                ICNT=0
                ITAG=0

                IF(KK.LT.NTPREV/2) THEN
                   NI1=1
                   NI2=NH
                   NJMP=1
                   IF(IDR.EQ.1) NLAST=NC
                   IF(IDR.EQ.2) NLAST=1
                ELSE
                   NI1=NH
                   NI2=1
                   NJMP=-1
                   IF(IDR.EQ.1) NLAST=NHP
                   IF(IDR.EQ.2) NLAST=NH
                END IF

                DO 60 N=NI1,NI2,NJMP

                   IF(IDR.EQ.1) THEN
                      N1=NH+N
                      N2=N1+1
                      NIDX=N1
                   ELSE IF(IDR.EQ.2) THEN
                      N1=NH-N+2
                      N2=N1-1
                      NIDX=N2
                   END IF

C.................Check for split panels chordwise
                   XN(1)=(XB(N1,M)+XB(N1,M+1))/2.
                   YN(1)=(YB1(N1,M)+YB1(N1,M+1))/2.
                   ZN(1)=(ZB1(N1,M)+ZB1(N1,M+1))/2.
                   XN(2)=(XB(N2,M)+XB(N2,M+1))/2.
                   YN(2)=(YB1(N2,M)+YB1(N2,M+1))/2.
                   ZN(2)=(ZB1(N2,M)+ZB1(N2,M+1))/2.
                   
                   CALL CKSPLIT(YFS,XN,YN,ZN,SFAC,I1)

C.................If IMG=1, only perform calculations on fully
C.................submerged panels. (JY012300)
                   IF(IMG.EQ.1) THEN
                      CK1=YB1(N1,M)-YFS
                      CK2=YB1(N1,M+1)-YFS
                      CK3=YB1(N2,M)-YFS
                      CK4=YB1(N2,M+1)-YFS
                      CKTOL=1.E-5
                      IF(CK1.GT.CKTOL.OR.CK2.GT.CKTOL.OR.CK3.GT.CKTOL.
     *                     OR.CK4.GT.CKTOL) THEN
                         SFAC=0.0
                      END IF
                   END IF

                   IF(N.EQ.NI1) THEN

                      NIDX1=NIDX
                      IF(SFAC.GE.0.5) THEN
                         IFE=1
                         ICNT=ICNT+1
                         MS(ICNT)=NIDX1
                         ICB(M,IDR,KK)=1
                      ELSE
                         IFE=0
                      END IF

                   ELSE

                      IF(SFAC.GT.0.0.AND.SFAC.LT.1.0) THEN
                         IF(ICNT.EQ.0) ICB(M,IDR,KK)=1
                         ICNT=ICNT+1
                         IF(ICNT.EQ.1) ITAG=1
                         IF(SFAC.GE.0.5) THEN
                            MS(ICNT)=NIDX
                         ELSE
                            IF(I1.EQ.1) THEN
                               IF(IDR.EQ.1) MS(ICNT)=NIDX-1
                               IF(IDR.EQ.2) MS(ICNT)=NIDX+1
                            ELSE IF(I1.EQ.2) THEN
                               IF(IDR.EQ.1) MS(ICNT)=NIDX+1
                               IF(IDR.EQ.2) MS(ICNT)=NIDX-1
                            END IF
                         END IF
                         IF(ICNT.EQ.2) GO TO 55
                      ELSE
                         IF(SFAC.NE.SFACP.AND.ITAG.NE.1) THEN
                            IF(ICNT.EQ.0) ICB(M,IDR,KK)=1
                            ICNT=ICNT+1
                            IF(SFACP.EQ.1.0) MS(ICNT)=NIDXP
                            IF(SFACP.EQ.0.0) MS(ICNT)=NIDX
                            IF(ICNT.EQ.2) GO TO 55
                         ELSE IF(SFAC.EQ.SFACP.AND.N.EQ.NI2) THEN
                            IF(SFAC.EQ.1.0) THEN
                               IF(ICNT.EQ.0) ICB(M,IDR,KK)=1
                               ICNT=ICNT+1
                               MS(ICNT)=NLAST
                            END IF
                         END IF
                      END IF

                   END IF

                   NIDXP=NIDX
                   SFACP=SFAC
                   
 60             CONTINUE

 55             CONTINUE

                IF(MS(1).EQ.MS(2)) ICB(M,IDR,KK)=0

                IF(ICB(M,IDR,KK).EQ.1) THEN
                   IF(IDR.EQ.1) THEN
                      IC(1,M,KK)=MIN(MS(1),MS(2))
                      IC(2,M,KK)=MAX(MS(1),MS(2))
                   ELSE
                      IW(1,M,KK)=MIN(MS(1),MS(2))
                      IW(2,M,KK)=MAX(MS(1),MS(2))
                   END IF

                   IF(ISC.EQ.1) THEN
                      IF(IDR.EQ.2.AND.IW(2,M,KK).LT.N0(2)) 
     *                     ICB(M,2,KK)=0
                      IF(IDR.EQ.1.AND.IC(1,M,KK).GE.N0(1))
     *                     ICB(M,1,KK)=0
                   END IF

                END IF
                
 50          CONTINUE

 40       CONTINUE

          DO M=1,MR

C...........Flag the submerged panels
             IF(ICB(M,1,KK).NE.ICB(M,2,KK)) THEN
                DO IDR=1,2
                   ICB(M,IDR,KK)=0
                END DO
             END IF

             DO N=1,NC
                ISUBM(INDEXB(N,M),KK)=0
             END DO

             DO IDR=1,2
                ISUB1(M,IDR,KK)=0
             END DO

             IF(ICB(M,2,KK).EQ.1) THEN
                DO N=IW(1,M,KK),IW(2,M,KK)
                   ISUBM(INDEXB(N,M),KK)=1
                   ISUB1(M,2,KK)=ISUB1(M,2,KK)+1
                END DO
             END IF
             IF(ICB(M,1,KK).EQ.1) THEN
                DO N=IC(1,M,KK),IC(2,M,KK)
                   ISUBM(INDEXB(N,M),KK)=1  
                   ISUB1(M,1,KK)=ISUB1(M,1,KK)+1    
                END DO
             END IF

C..........Sum the total number of submerged elements on the blade
             DO IDR=1,2
                IF(ICB(M,IDR,KK).EQ.1) THEN
                   IF(IDR.EQ.1) THEN
                      N1=IC(1,M,KK)
                      N2=IC(2,M,KK)
                      KI=1
                   ELSE
                      N1=IW(2,M,KK)
                      N2=IW(1,M,KK)
                      KI=-1
                   END IF

                   DO N=N1,N2,KI
                      NEL(KK)=NEL(KK)+1
                   END DO
                END IF
             END DO

          END DO

C-----------------------------------------------------------------------
C     Check for submerged panels on the key wake
C-----------------------------------------------------------------------
          DO 80 M=1,MR+1
             DO 90 N=1,NWMIN+1
                YW1(N,M)=YW(N,M)*COS(DT1)-ZW(N,M)*SIN(DT1)
                ZW1(N,M)=YW(N,M)*SIN(DT1)+ZW(N,M)*COS(DT1)
 90          CONTINUE
 80       CONTINUE

          DO M=1,MR
             ICW(M,KK)=0
             MSW(M,KK)=0

             IF(ICB(M,1,KK).EQ.0.OR.ICB(M,2,KK).EQ.0) GO TO 105
             IF(IC(2,M,KK).NE.NC.OR.IW(1,M,KK).NE.1) GO TO 105

             NI1=1
             NI2=NWMIN
             NJMP=1
             
             DO 110 N=NI1,NI2,NJMP
                
                NIDX=N

                XN(1)=(XW(N,M)+XW(N,M+1))/2.
                YN(1)=(YW1(N,M)+YW1(N,M+1))/2.
                ZN(1)=(ZW1(N,M)+ZW1(N,M+1))/2.
                XN(2)=(XW(N+1,M)+XW(N+1,M+1))/2.
                YN(2)=(YW1(N+1,M)+YW1(N+1,M+1))/2.
                ZN(2)=(ZW1(N+1,M)+ZW1(N+1,M+1))/2.

                CALL CKSPLIT(YFS,XN,YN,ZN,SFAC,I1)

                IF(N.EQ.NI1) THEN

                   IF(SFAC.GE.0.5) THEN
                      ICW(M,KK)=1
                      MSW(M,KK)=1
                   ELSE
                      GO TO 105
                   END IF

                ELSE

                   IF(SFAC.GT.0.0.AND.SFAC.LT.1.0) THEN
                      ICW(M,KK)=1
                      IF(SFAC.GE.0.5) THEN
                         MSW(M,KK)=N
                      ELSE
                         IF(I1.EQ.1) THEN
                            MSW(M,KK)=N-1
                         ELSE IF(I1.EQ.2) THEN
                            MSW(M,KK)=N+1
                         END IF
                      END IF
                      GO TO 105
                   ELSE
                      IF(SFAC.NE.SFACP) THEN
                         ICW(M,KK)=1
                         IF(SFACP.EQ.1.0) MSW(M,KK)=NIDXP
                         IF(SFACP.EQ.0.0) MSW(M,KK)=N
                         GO TO 105
                      ELSE IF(SFAC.EQ.SFACP.AND.N.EQ.NI2) THEN
                         IF(SFAC.EQ.1.0) THEN
                            ICW(M,KK)=1
                            MSW(M,KK)=NWMIN
                         END IF
                      END IF
                   END IF

                END IF

                NIDXP=NIDX
                SFACP=SFAC

 110         CONTINUE
          
 105         CONTINUE

C...........If IMG=1, only perform calculations on fully
C...........submerged panels. (JY012300)
             IF(IMG.EQ.1.AND.ICW(M,KK).EQ.1) THEN
 106            N1=MSW(M,KK)
                N2=MSW(M,KK)+1

                CK1=YW1(N1,M)-YFS
                CK2=YW1(N1,M+1)-YFS
                CK3=YW1(N2,M)-YFS
                CK4=YW1(N2,M+1)-YFS
                CKTOL=1.E-5
                IF(CK1.GT.CKTOL.OR.CK2.GT.CKTOL.OR.CK3.GT.CKTOL.OR.
     *               CK4.GT.CKTOL) THEN
                   MSW(M,KK)=MSW(M,KK)-1
                   IF(MSW(M,KK).EQ.0) THEN
                      ICW(M,KK)=0
                      GO TO 108
                   END IF
                   GO TO 106
                ELSE
                   GO TO 108
                END IF
             END IF

 108         CONTINUE

C...........Sum the number of submerged elements on the wake
             IF(ICW(M,KK).EQ.1) THEN
                N1=1
                N2=MSW(M,KK)           
                DO N=N1,N2
                   NEL(KK)=NEL(KK)+1
                END DO
             END IF
          END DO

C-----------------------------------------------------------------------
C     Check for submerged panels on the key hub strip
C-----------------------------------------------------------------------
          DO M=1,MHBT+1
             DO N=1,NHBX+1
                YH1(N,M)=YH(N,M)*COS(DT1)-ZH(N,M)*SIN(DT1)
                ZH1(N,M)=YH(N,M)*SIN(DT1)+ZH(N,M)*COS(DT1)
             END DO
          END DO

          DO N=1,NHBX
             ISUBH(N,KK)=0
             DO M=1,MHBT
                I=INDEXH(N,M)

C..............Check for split panels chordwise
                XN(1)=(XH(N,M)+XH(N,M+1))/2.
                YN(1)=(YH1(N,M)+YH1(N,M+1))/2.
                ZN(1)=(ZH1(N,M)+ZH1(N,M+1))/2.
                XN(2)=(XH(N+1,M)+XH(N+1,M+1))/2.
                YN(2)=(YH1(N+1,M)+YH1(N+1,M+1))/2.
                ZN(2)=(ZH1(N+1,M)+ZH1(N+1,M+1))/2.
                
                CALL CKSPLIT(YFS,XN,YN,ZN,SFAC,I1)

C..............If IMG=1, only perform calculations on fully
C..............submerged panels. (JY012300)
                IF(IMG.EQ.1) THEN
                   CK1=YH1(N,M)-YFS
                   CK2=YH1(N,M+1)-YFS
                   CK3=YH1(N+1,M)-YFS
                   CK4=YH1(N+1,M+1)-YFS
                   CKTOL=1.E-5
                   IF(CK1.GT.CKTOL.OR.CK2.GT.CKTOL.OR.CK3.GT.CKTOL.
     *                  OR.CK4.GT.CKTOL) THEN
                      SFAC=0.0
                   END IF
                END IF

C..............Sum the number of submerged panels on the hub.  Also
C..............tag submerged panels
                IF(SFAC.EQ.1.) THEN
                   ISUBM(I,KK)=1
                   ISUBH(N,KK)=ISUBH(N,KK)+1
                   NEL(KK)=NEL(KK)+1
                ELSE
                   ISUBM(I,KK)=0
                END IF
             END DO
          END DO

 10    CONTINUE

C.....Plot submerged panels
       CALL SUBMPLT(NEL)

       RETURN
       END


       SUBROUTINE CKSPLIT(YFS,XN,YN,ZN,SFAC,I1)

       DIMENSION XN(2),YN(2),ZN(2)

C.....Let YN(I2)>YN(I1)
       IF(YN(2).GE.YN(1)) THEN
          I2=2
          I1=1
       ELSE
          I2=1
          I1=2
       END IF

C.....CUT is the portion of YN(I2)-YN(I1) that's cut by
C.....the free surface
       IF(YN(2).NE.YN(1)) THEN
          CUT=(YFS-YN(I1))/(YN(I2)-YN(I1))
          IF(CUT.LT.0.0) THEN
             SFAC=0.
          ELSE IF(CUT.GT.1.0) THEN
             SFAC=1.
          ELSE
             XSUB=XN(I1)+CUT*(XN(I2)-XN(I1))
             YSUB=YFS
             ZSUB=ZN(I1)+CUT*(ZN(I2)-ZN(I1))
             XLCUT=SQRT((XSUB-XN(I1))**2.+(YSUB-YN(I1))**2.
     *            +(ZSUB-ZN(I1))**2.)
             XL=SQRT((XN(1)-XN(2))**2.+(YN(1)-YN(2))**2.+
     *            (ZN(1)-ZN(2))**2.)
             SFAC=XLCUT/XL
          END IF
       ELSE
          IF(YN(2).LE.YFS) THEN
             SFAC=1.
          ELSE IF(YN(2).GT.YFS) THEN
             SFAC=0.
          END IF
       END IF

       IF(SFAC.NE.0.0.AND.ABS(SFAC).LE.1.E-5) THEN
          SFAC=0.0
       END IF
       IF(SFAC.NE.1.0.AND.ABS(1.0-SFAC).LE.1.E-5) THEN
          SFAC=1.0
       END IF

       RETURN
       END


