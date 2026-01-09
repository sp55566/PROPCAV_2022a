        SUBROUTINE DELR
************************************************************************
*                                                                      *
*  Subroutine DELR computes delta(r) for a given guess of l(r)         *
*                                                                      *
*  Author: Neal Fine  July 15, 1991                                    *
*                                                                      *
*  Date of last Revision                     Revision                  *
*  ---------------------                 ----------------              *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C-----------------------------------------------------------------------
C     set up the matrix equation
C-----------------------------------------------------------------------
      ISR=1
      IF(IFACE.EQ.2) ISR=2
      DTN=DELTAT*ADVCO/PI

      DO 10 M=1,MR         

         CH=HALF*(CHORD(M)+CHORD(M+1))

         NNWC(M)=0

C.......initialize variables.........................................
         DO IDR=1,2
            DO 30 K=1,4
               CT(M,K,IDR)=ZERO
               DT(M,K,IDR)=ZERO
               QW(M,K,IDR)=ZERO
               QT(M,K,IDR)=ZERO
 30         CONTINUE
            DLISP(M,IDR)=ZERO
            QSPR(M,IDR)=ZERO
            FRP(M,IDR)=ZERO
            FLP(M,IDR)=ZERO
            DZL(M,IDR)=ZERO
            DZR(M,IDR)=ZERO
            NSPP(M,IDR)=0
            FRS(M,IDR)=ZERO
            FLS(M,IDR)=ZERO
            NSPS(M,IDR)=0
            SOP(M,IDR)=ZERO
            JCV(M,IDR)=0
            NWC(M,IDR)=0
         END DO

         DO 20 I=1,ISR
            IF((IFACE.EQ.0).OR.(I.EQ.1.AND.IFACE.EQ.2)) THEN
               IDR=1
               KI=1
               ISF=0
            ELSE
               IDR=2
               KI=-1
               ISF=1
            END IF
            AK1=FLOAT(KI)
         
            IF(ISP.EQ.1) THEN
              IF(ISUB1(M,IDR,IDXREV).EQ.0) GO TO 20
            ENDIF

C..........Special treatment for surface piering propellers (JY112199)
            IF(ISP.EQ.1) THEN
               NWDIR(M)=1
               IF(ICB(M,1,IDXREV).EQ.1) THEN
                  M0(M,1)=NHP+NLEP(M,IDXREV,1)
                  JCV(M,1)=MAX(0,IC(2,M,IDXREV)-M0(M,1)+1)
                  LCV(M,1)=M0(M,1)+JCV(M,1)
                  IF(JCV(M,1).EQ.0) THEN
                     NOCAV(M,IDXREV,1)=1
                  ELSE
                     NOCAV(M,IDXREV,1)=0
                     IF(IC(2,M,IDXREV).EQ.NC.AND.
     *                    LCV(M,1).EQ.NCP) SOP(M,1)=ONE
                  END IF 
               END IF
               IF(ICW(M,IDXREV).EQ.1) THEN
                  IF(JCV(M,1).GT.0) SOP(M,1)=ONE
                  NWC(M,1)=N1SUB+(MSW(M,IDXREV)-1)*NWSUB1
                  NNWC(M)=NWC(M,1)
               END IF  

            ELSE  !         ISP!= 1

C.............find the true cavity length and define relavent indices
               CAVL1=CAVL(M,IDR)
               CALL POSIT(CAVL1,M,IDR)

C.............print warning and stop program if ICON=8, but it's not
C.............supercavity on both sides. (JY110100)
               IF(ICON.EQ.8) THEN
 811              FORMAT('Program stop: Not supercavitating!')
                  IF(SOP(M,IDR).EQ.ZERO) THEN
                     WRITE(*,*) 'Problem at strip #',M,' ,IDR=',IDR
                     WRITE(*,811)
                     STOP
                  END IF
               END IF

C.....If statement added for face & back supercavitation (JY053000)
               IF(IFACE.EQ.2.AND.IDR.EQ.2) THEN
                  IF(NWC(M,1).GT.0.AND.NWC(M,2).GT.0) THEN
                     IF(NWC(M,1).GE.NWC(M,2)) THEN
                        NWC(M,2)=NWC(M,1)
                        CAVL(M,2)=CAVL(M,1)
                        SOP(M,2)=SOP(M,1)
                        NNWC(M)=NWC(M,1)
                        NWDIR(M)=1
                        IF(NSPS(M,1).GT.0) THEN
                           NSPS(M,2)=NSPS(M,1)
                           FRS(M,2)=FRS(M,1)
                           FLS(M,2)=FLS(M,1)
                           DZL(M,2)=DZL(M,1)
                           DZR(M,2)=DZR(M,1)
                        END IF
                     ELSE
                        NWC(M,1)=NWC(M,2)
                        CAVL(M,1)=CAVL(M,2)
                        SOP(M,1)=SOP(M,2)
                        NNWC(M)=NWC(M,2)
                        NWDIR(M)=2
                        IF(NSPS(M,2).GT.0) THEN
                           NSPS(M,1)=NSPS(M,2)
                           FRS(M,1)=FRS(M,2)
                           FLS(M,1)=FLS(M,2)
                           DZL(M,1)=DZL(M,2)
                           DZR(M,1)=DZR(M,2)
                        END IF
                     END IF
                  END IF
               END IF

            END IF

            M0M=M0(M,IDR)
            M0M1=M0M-1
            JCAV=JCV(M,IDR)
            LCAV=M0M+KI*(JCAV+NSPP(M,IDR))
            IF(NSPP(M,IDR).EQ.1) ISPLIT=LCAV-KI*1-ISF

C-----------------------------------------------------------------------
C     Redefine the source strength, to account for the source/sink left
C     by a collapsing cavity or moving cavity detachment.       JY092699
C-----------------------------------------------------------------------
            ISW=0

            IF(IAN .EQ. 2 .AND. ISCAV .EQ. 2
     %           .AND. NREV .EQ. NNREV) ISW = 1

C..........Temporarily disregard DH/DT term for ISP=1 & IFACE=2
            IF(ISP.EQ.1) ISW=1

            IF(ISW.EQ.1) GO TO 40

            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.2) GO TO 40

            IF(JCAV.EQ.0) THEN
               DO N=1,NH
                  DHDT=HALF*(HTP(N,M,IDR)+HTP(N+1,M,IDR))/DTN
                  J=INDEXB(NH+ISF+KI*N,M)
                  DPDNC(J)=DPDNC(J)-DHDT
               END DO
            ELSE
               DO N=1,NLEP(M,IDXREV,IDR)
                  DHDT=HALF*(HTP(N,M,IDR)+HTP(N+1,M,IDR))/DTN
                  J=INDEXB(NH+ISF+KI*N,M)
                  DPDNC(J)=DPDNC(J)-DHDT
               END DO
               N1=M0M1+ISF+KI*(JCAV+NSPP(M,IDR)+1)
               IF(IDR.EQ.1) THEN
                  N2=NC
                  NJ=1
                  N3=N1-NH
               ELSE IF(IDR.EQ.2) THEN
                  N2=1
                  NJ=-1
                  N3=NHP-N1
               END IF
               DO N=N1,N2,NJ
                  DHDT=HALF*(HTP(N3,M,IDR)+HTP(N3+1,M,IDR))/DTN
                  N3=N3+1
                  J=INDEXB(N,M)
                  DPDNC(J)=DPDNC(J)-DHDT
               END DO
            END IF

 40         CONTINUE

            IF(JCV(M,IDR).EQ.0) GOTO 60

            NN1=M0M1+ISF

C..........determine the extrapolation constants for phi0...............
            IWLE=NH

            IF(IFACE.EQ.2.AND.JCV(M,1).GT.0.AND.JCV(M,2).GT.0) 
     *           IWLE=NLEP(M,IDXREV,1)+NLEP(M,IDXREV,2)
            IF(ISP.EQ.1) THEN
               IF(IC(1,M,IDXREV)-1.EQ.IW(2,M,IDXREV)) THEN
                  IWLE=NLEP(M,IDXREV,IDR)+
     *                 IW(2,M,IDXREV)-IW(1,M,IDXREV)+1
               ELSE
                  IWLE=NLEP(M,IDXREV,IDR)-(IC(1,M,IDXREV)-NHP)
               END IF
            END IF

            IF(IWLE.GE.3) THEN
               S2=HALF*(DS(NN1-KI*2,M)+DS(NN1-KI*1,M))
               S1=S2+HALF*(DS(NN1-KI*1,M)+DS(NN1,M))
               S0=S1+HALF*DS(NN1,M)
               CALL EXTRAP(S0,S1,S2,C1,C2,C3,C4)
               CT(M,1,IDR)=C1
               CT(M,2,IDR)=C2
               CT(M,3,IDR)=C3
               CT(M,4,IDR)=C4
            ELSE IF(IWLE.EQ.2) THEN
               S1=HALF*(DS(NN1-KI*1,M)+DS(NN1,M))
               S0=S1+HALF*DS(NN1,M)
               CALL EXTRAP2(S0,S1,C1,C2,C4)
               CT(M,1,IDR)=C1
               CT(M,2,IDR)=C2
               CT(M,3,IDR)=0.
               CT(M,4,IDR)=C4               
            ELSE IF(IWLE.EQ.1) THEN
               S0=HALF*DS(NN1,M)
               CT(M,1,IDR)=1.
               CT(M,2,IDR)=0.
               CT(M,3,IDR)=0.
               CT(M,4,IDR)=S0
            END IF

C..........define pressure-law quantities...............................
            SLPR=AK1*(SZ(M0M1+ISF+KI*(JCAV+1-ISF),M)-SZ(NHP,M))
            STPR=SLPR*(ONE-RLAM(M,IDR))

C..........redefine last SPZ to be at center of the left split panel....
            IF(NSPP(M,IDR).EQ.1) THEN
               J1=M0M1+ISF+KI*JCAV
               SPZ(J1+KI*1,M)=SPZ(J1,M)+AK1*HALF*(DS(J1,M)+
     *              DZL(M,IDR))
               SLPR=AK1*(SPZ(J1+KI*1,M)+HALF*DZL(M,IDR)-
     *              SZ(NHP,M))  
               STPR=SLPR*(ONE-RLAM(M,IDR))
            ENDIF
            IF(NSPS(M,IDR).EQ.1)THEN
               J1=NC+NWC(M,IDR)
               SPZ(J1+1,M)=SPZ(J1,M)+HALF*(DZW(NWC(M,IDR),M)+DZL(M,IDR))
            ENDIF

C..........determine PHI1 for the dynamic boundary condition............
            CALL COMPPHI1(M,IDR)

C..........define the extrapolation constants for split-panel wake .....
C..........source.......................................................
            IF(NSPS(M,IDR).EQ.1) CALL QWEXTRAP(M,IDR)

            IF(NSPP(M,IDR).EQ.0) GOTO 60

C..........source strength on the right side of the split panel due to..
C..........collapsing cavity............................................
            J=INDEXB(ISPLIT,M)
            QSPR(M,IDR)=DPDN(J)

            ISW=0

            IF(IAN .EQ. 2 .AND. ISCAV .EQ. 2
     %           .AND. NREV .EQ. NNREV) ISW = 1

C..........Temporarily disregard DH/DT term for ISP=1 & IFACE=2
            IF(ISP.EQ.1) ISW=1

            IF(ISW.EQ.1) GO TO 70

            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.2) GO TO 70

            N=NLEP(M,IDXREV,IDR)+JCAV+1
            HTENM1=HTP(N,M,IDR)+(HTP(N+1,M,IDR)-HTP(N,M,IDR))
     *           *DZL(M,IDR)/DZ(ISPLIT,M)
            DHDT=HALF*(HTENM1+HTP(N+1,M,IDR))/DTN
            QSPR(M,IDR)=QSPR(M,IDR)-DHDT

 70         CONTINUE

C..........dipole strength (- phi0) on the left side of the split panel.
            DLISP(M,IDR)=PHI1(JCAV+1,M,IDR)

C..........extrapolation constants for the right dipole.................
            CALL DTEXTRAP(M,IDR)

C..........extrapolation constant for split panel left source on the ...
C..........foil.........................................................
            CALL QTEXTRAP(M,IDR)

C..........extrapolate dpdv to left-split panel.........................
            CALL DPDVEXTRAP(M,IDR)

 60         CONTINUE

 20      CONTINUE

 10   CONTINUE

      NBW=0
      DO 80 M=1,MR

C........determine the side of the blade that is supercaviting..........
         IF((NWC(M,2).EQ.0.AND.NWC(M,1).GT.0).OR.IFACE.EQ.0) THEN
            NNWC(M)=NWC(M,1)
            NWDIR(M)=1
         ELSE IF((NWC(M,1).EQ.0.AND.NWC(M,2).GT.0).OR.IFACE.EQ.1) THEN
            NNWC(M)=NWC(M,2)
            NWDIR(M)=2
         END IF

C-----------------------------------------------------------------------
C     Add the effect of the trailing sinks (due to collapsing of the
C     supercavity) on the wake sheet after where the cavity ended.  
C     Introduced new variable QSSR(M) to represent the strength of 
C     the trailing sink on the right side of the split panel.   JY060800
C-----------------------------------------------------------------------
         QSSR(M)=ZERO
         IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.2) GO TO 90

c         ISW=1
         ISW=0  

C.......Temporarily disregard DH/DT term for surface piercing
C.......propellers. (JY061800)
         IF(ISP.EQ.1.OR.ISW.EQ.1) GO TO 90
         
         NSC=NNWC(M)
         IF(NSC.GT.0) THEN
            IDR=NWDIR(M)
            NSP=NSPS(M,IDR)
         ELSE
            NSP=0
         END IF
         
         IF(NSC.GT.0.AND.NSP.EQ.1) THEN
            N=NSC+1
            HTENM1=HTWP(N,M)+(HTWP(N+1,M)-HTWP(N,M))
     *           *DZL(M,IDR)/DZW(NSC+1,M)
            DHDT=HALF*(HTENM1+HTWP(N+1,M))/DTN
            QSSR(M)=QSSR(M)-DHDT
         END IF
         DO N=NSC+NSP+1,NTRA
            J=(MR-M)*NTRA+N
            DHDT=HALF*(HTWP(N,M)+HTWP(N+1,M))/DTN
            SORW(J)=-DHDT
         END DO

 90      CONTINUE

C-----------------------------------------------------------------------
C     determine the number of applications of Greens formula 
C-----------------------------------------------------------------------

C.......Only solve the system of equations for panels that are 
C.......submerged for surface piercing propellers. (JY112199)
         IF(ISP.EQ.1) THEN
            NPERM(M)=ISUB1(M,1,IDXREV)+ISUB1(M,2,IDXREV)+NNWC(M)
            NBW=NBW+NPERM(M)
         ELSE
            NBW=NBW+NC-NSPP(M,1)-NSPP(M,2)+NNWC(M)
         END IF

 80   CONTINUE

C....Only solve for the submerged panels in the case of ISP=1.(JY112499)
      IF(ISP.EQ.1) THEN
         IF(NBW.EQ.0) THEN
            CALL CAVSOL0_SP
            GO TO 100            
         ELSE
            CALL CAVSET_SP
            CALL CAVSOL_SP
            CALL CAVHT_SP
            CALL PRSDIF_SP
            GO TO 100
         END IF
      END IF

C-----------------------------------------------------------------------
C     set up the linear system
C-----------------------------------------------------------------------
      CALL CAVSET

C-----------------------------------------------------------------------
C     solve the linear system
C-----------------------------------------------------------------------
      CALL CAVSOL

C-----------------------------------------------------------------------
C     compute the cavity pressure                               JY080500
C-----------------------------------------------------------------------
      CALL PRSDIF(1)

C-----------------------------------------------------------------------
C     iterate to determine the correct DPDVB.                   JY080600
C-----------------------------------------------------------------------
c      IF(IFACE.NE.2) CALL ITER2

C-----------------------------------------------------------------------
C     compute delta(r)
C-----------------------------------------------------------------------
      CALL CAVHT

 100  CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine DELR>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END
