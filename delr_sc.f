        SUBROUTINE DELR_SC
************************************************************************
*                                                                      *
*  Subroutine DELR computes delta(r) for a given guess of l(r)         *
*                                                                      *
*  Author: Neal Fine  July 15, 1991                                    *
*                                                                      *
*  Date of last Revision                     Revision                  *
*  ---------------------                 ----------------              *
*  JY092601                Copied from delr.f.  Modified for ISC=1.    *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C-----------------------------------------------------------------------
C     set up the matrix equation
C-----------------------------------------------------------------------
      DTN=DELTAT*ADVCO/PI

      DO 10 M=1,MR         

         CH=HALF*(CHORD(M)+CHORD(M+1))
         NNWC(M)=0

C.......initialize variables.........................................
         DO IDR=1,2

            DO 30 K=1,4
               CT(M,K,IDR)=ZERO
               CT2(M,K,IDR)=ZERO
               DT(M,K,IDR)=ZERO
               QW(M,K,IDR)=ZERO
               QT(M,K,IDR)=ZERO
 30         CONTINUE
            DLISP(M,IDR)=ZERO
            QSPR(M,IDR)=ZERO
            FRP(M,IDR)=ZERO
            FLP(M,IDR)=ZERO
            NSPP(M,IDR)=0
            FRS(M,IDR)=ZERO
            FLS(M,IDR)=ZERO
            DZL(M,IDR)=ZERO
            DZR(M,IDR)=ZERO
            NSPS(M,IDR)=0
            SOP(M,IDR)=ZERO
            JCV(M,IDR)=0
            NWC(M,IDR)=0

            IF(ISP.EQ.1) THEN
              IF (ISUB1(M,IDR,IDXREV).EQ.0) GO TO 60
            END IF

            IF(IDR.EQ.1) THEN
               KI=1
               ISF=0
            ELSE
               KI=-1
               ISF=1
            END IF
            AK1=FLOAT(KI)
         
C..........find the true cavity length and define relavent indices
            IF(ISP.EQ.0) THEN
               CAVL1=CAVL(M,IDR)
               CAVL2=CAVLT(M)
               CALL POSIT_SC1(CAVL1,M,IDR)
               CALL POSIT_SC2(CAVL2,M)
            ELSE
               NWDIR(M)=1
               IF(ICB(M,1,IDXREV).EQ.1) THEN
                  M0(M,1)=NHP+NLEP(M,IDXREV,1)
                  NMAX=MIN(IC(2,M,IDXREV),N0(1)-1)
                  JCV(M,1)=MAX(0,NMAX-M0(M,1)+1)
                  LCV(M,1)=M0(M,1)+JCV(M,1)
                  IF(JCV(M,1).EQ.0) THEN
                     NOCAV(M,IDXREV,1)=1
                  ELSE
                     NOCAV(M,IDXREV,1)=0
                     IF(ISUBM(INDEXB(N0(1),M),IDXREV).EQ.1.
     *                    AND.LCV(M,1).EQ.N0(1)) SOP(M,1)=ONE
                  END IF                     
               END IF
               IF(ICW(M,IDXREV).EQ.1) THEN
                  NWC(M,1)=N1SUB+(MSW(M,IDXREV)-1)*NWSUB1                 
                  NNWC(M)=NWC(M,1)
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
            ISW=1

            IF(IAN .EQ. 2 .AND. ISCAV .EQ. 2
     %           .AND. NREV .EQ. NNREV) ISW = 1

C..........Temporarily disregard DH/DT term for surface piercing
C..........propellers. (JY112199)
            IF(ISP.EQ.1) ISW=1

            IF(ISW.EQ.1) GO TO 40

            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.4) GO TO 40

            IF(JCAV.EQ.0) THEN
               DO N=1,NHOLD
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
                  N2=N0(1)-1
                  NJ=1
               ELSE IF(IDR.EQ.2) THEN
                  N2=N0(2)
                  NJ=-1
               END IF
               N3=NLEP(M,IDXREV,IDR)+JCAV+NSPP(M,IDR)+1
               DO N=N1,N2,NJ
                  DHDT=HALF*(HTP(N3,M,IDR)+HTP(N3+1,M,IDR))/DTN
                  N3=N3+1
                  J=INDEXB(N,M)
                  DPDNC(J)=DPDNC(J)-DHDT
               END DO
            END IF

 40         CONTINUE

            IF(ISP.EQ.1) THEN
               IFL1=0
               DO N=1,NSR2
                  IF(IDR.EQ.2) THEN
                     N1=N0(2)-N
                  ELSE
                     N1=N0(1)-1+N
                  END IF
                  IF(ISUBM(INDEXB(N1,M),IDXREV).EQ.1) THEN
                     IFL1=1
                     GO TO 112
                  END IF
               END DO
               GO TO 111
 112           CONTINUE
            END IF

C..........Determine extrapolation constant for PHI0T
            IF(SOP(M,IDR).EQ.ZERO) THEN
               IF(ISP.EQ.0) THEN
                  IWLE=NHOLD
                  IF(JCV(M,IDR).GT.0) IWLE=IWLE-(NLEP(M,IDXREV,IDR)
     *                 +JCV(M,IDR))
               ELSE
                  IWLE=0
                  IF(ICB(M,IDR,IDXREV).EQ.1) THEN
                     IF(IDR.EQ.2) THEN
                        IF(IW(2,M,IDXREV).GE.N0(2).AND.
     *                       IW(1,M,IDXREV).LT.N0(2)) IWLE=
     *                       IW(2,M,IDXREV)-N0(2)+1
                     ELSE
                        IF(IC(1,M,IDXREV).LT.N0(1).AND.
     *                       IC(2,M,IDXREV).GE.N0(1)) IWLE=
     *                       N0(1)-IC(1,M,IDXREV)
                     END IF
                  END IF
               END IF

               N0M=N0(IDR)
               N0M1=N0M-1
               NN1=N0M1+ISF
               IF(IWLE.GE.3) THEN
                  S2=HALF*(DS(NN1-KI*2,M)+DS(NN1-KI*1,M))
                  S1=S2+HALF*(DS(NN1-KI*1,M)+DS(NN1,M))
                  S0=S1+HALF*DS(NN1,M)
                  CALL EXTRAP(S0,S1,S2,C1,C2,C3,C4)
                  CT2(M,1,IDR)=C1
                  CT2(M,2,IDR)=C2
                  CT2(M,3,IDR)=C3
                  CT2(M,4,IDR)=C4
               ELSE IF(IWLE.EQ.2) THEN
                  S1=HALF*(DS(NN1-KI*1,M)+DS(NN1,M))
                  S0=S1+HALF*DS(NN1,M)
                  CALL EXTRAP2(S0,S1,C1,C2,C4)
                  CT2(M,1,IDR)=C1
                  CT2(M,2,IDR)=C2
                  CT2(M,3,IDR)=0.
                  CT2(M,4,IDR)=C4               
               ELSE IF(IWLE.EQ.1) THEN
                  S0=HALF*DS(NN1,M)
                  CT2(M,1,IDR)=1.
                  CT2(M,2,IDR)=0.
                  CT2(M,3,IDR)=0.
                  CT2(M,4,IDR)=S0
               END IF
            END IF

            IF(NSPS(M,IDR).EQ.1)THEN
               J1=NC+NWC(M,IDR)
               SPZ(J1+1,M)=SPZ(J1,M)+HALF*DZW(NWC(M,IDR),M)*
     *              (ONE+FLS(M,IDR))
            ENDIF

            CALL COMPPHI1_SC2(M,IDR)

C..........define the extrapolation constants for split-panel wake .....
C..........source.......................................................
            IF(NSPS(M,IDR).EQ.1) CALL QWEXTRAP(M,IDR)

 111        CONTINUE

            IF(JCV(M,IDR).EQ.0) GOTO 60

            NN1=M0M1+ISF

C..........determine the extrapolation constants for phi0...............
            IWLE=NHOLD
            IF(IFACE.EQ.2.AND.JCV(M,1).GT.0.AND.JCV(M,2).GT.0) 
     *           IWLE=NLEP(M,IDXREV,1)+NLEP(M,IDXREV,2)
            IF(ISP.EQ.1) THEN
               IF(IDR.EQ.2) THEN
                  IWLE=0
               ELSE
                  IF(IC(1,M,IDXREV)-1.EQ.IW(2,M,IDXREV)) THEN
                     IWLE=NLEP(M,IDXREV,IDR)+
     *                    IW(2,M,IDXREV)-IW(1,M,IDXREV)+1
                  ELSE
                     IWLE=NLEP(M,IDXREV,IDR)-(IC(1,M,IDXREV)-NHP)
                 END IF
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

C..........determine PHI1 for the dynamic boundary condition............
            CALL COMPPHI1_SC1(M,IDR)

            IF(NSPP(M,IDR).EQ.0) GOTO 60

C..........source strength on the right side of the split panel due to..
C..........collapsing cavity............................................
            J=INDEXB(ISPLIT,M)
            QSPR(M,IDR)=DPDN(J)

            ISW=1

            IF(IAN .EQ. 2 .AND. ISCAV .EQ. 2
     %           .AND. NREV .EQ. NNREV) ISW = 1

C..........Temporarily disregard DH/DT term for surface piercing
C..........propellers. (JY061800)
            IF(ISP.EQ.1) ISW=1

            IF(ISW.EQ.1) GO TO 70

            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.4) GO TO 70

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

         END DO

 10   CONTINUE

      NBW=0

      DO 80 M=1,MR

C-----------------------------------------------------------------------
C     Add the effect of the trailing sinks (due to collapsing of the
C     supercavity) on the wake sheet after where the cavity ended.  
C     Introduced new variable QSSR(M) to represent the strength of 
C     the trailing sink on the right side of the split panel.   JY060800
C-----------------------------------------------------------------------
         QSSR(M)=ZERO
         IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.3) GO TO 90

         ISW=1

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
            HTENM1=HTWP(N,M)+(HTWP(N+1,M)-HTWP(N,M))*FLS(M,IDR)
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
C-----------------------------------------------------------------------
C     determine the number of applications of Greens formula for ISP=1
C-----------------------------------------------------------------------
         IF(ISP.EQ.1) THEN
            NPERM(M)=ISUB1(M,1,IDXREV)+ISUB1(M,2,IDXREV)+NNWC(M)
            NBW=NBW+NPERM(M)
         ELSE
            NBW=NBW+NC-NSPP(M,1)-NSPP(M,2)+NNWC(M)
         END IF

 80   CONTINUE

C....Only solve for the submerged panels in the case of ISP=1.(JY112499)
      IF(ISP.EQ.1) THEN
         NTSOL=NBW
         IF(IHUB.NE.0) THEN
            DO I1=1,NPANH
               I=I1+NPANB
               IF(ISUBM(I,IDXREV).EQ.1) THEN
                  NTSOL=NTSOL+1
               END IF
            END DO
         END IF

         IF(NTSOL.EQ.0) THEN
            CALL CAVSOL0_SP
            GO TO 100            
         ELSE
            CALL CAVSET_SP
            CALL CAVSOL_SP 
            IF(NBW.GT.0) THEN
               CALL CAVHT_SP
               CALL PRSDIF_SP
            END IF
            GO TO 100
         END IF
      END IF

C-----------------------------------------------------------------------
C     set up the linear system
C-----------------------------------------------------------------------
      CALL CAVSET_SC

C-----------------------------------------------------------------------
C     solve the linear system
C-----------------------------------------------------------------------
      CALL CAVSOL_SC

C-----------------------------------------------------------------------
C     compute the cavity pressure                               JY080500
C-----------------------------------------------------------------------
      CALL PRSDIF_SC(1)

C-----------------------------------------------------------------------
C     compute delta(r)
C-----------------------------------------------------------------------
      CALL CAVHT_SC

 100  CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine DELR>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END

