C-----------------------------------------------------------------------------C
C                       CREATED BY YIRAN SU 10/14/2015                        C
C                 FOR COUPLING PROPCAV WITH AXISYMETRIC FLUENT                C
C                                  PC2NS                                      C
C-----------------------------------------------------------------------------C
C                       MODIFIED BY YIRAN SU 07/08/2016                       C
C                   TO HANDLE UNSTEADY/NON-AXISYMETRIC CASES                  C
C-----------------------------------------------------------------------------C

      MODULE M_PC2NS
      IMPLICIT NONE
      INTEGER PC_SPR ! time step per revolution
      INTEGER PC_NNN
      INTEGER PC_MAXIT
      REAL,ALLOCATABLE:: PC_BFX(:,:,:),PC_BFR(:,:,:),PC_BFT(:,:,:)
      REAL,ALLOCATABLE:: PC_SRCS(:,:,:)
      REAL,ALLOCATABLE:: PC_IVX(:,:,:,:),PC_IVR(:,:,:,:),PC_IVT(:,:,:,:)
      END MODULE M_PC2NS
      SUBROUTINE PC2NS
      USE M_BLADEM
      USE M_PC2NS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      REAL,DIMENSION(NHP,MRP)::BXX,BRR,BYY,BZZ
      REAL,DIMENSION(NH,MR)::BFX,BFR,BFT,SRCS
      REAL,DIMENSION(NC,MR)::INDX,INDR,INDT
      REAL,ALLOCATABLE:: DTT(:,:,:)
C ----- allocate and initiate M_PC2NS module
      IF (ISTEADY.NE.0) THEN
        PC_SPR=NTPREV
      ELSE
        PC_SPR=1
      END IF
      PC_NNN=15
      ALLOCATE(PC_BFX(NH,MR,PC_SPR),PC_BFR(NH,MR,PC_SPR))
      ALLOCATE(PC_BFT(NH,MR,PC_SPR),PC_SRCS(NH,MR,PC_SPR))
      ALLOCATE(PC_IVX(NH,MR,PC_NNN,PC_SPR),PC_IVR(NH,MR,PC_NNN,PC_SPR))
      ALLOCATE(PC_IVT(NH,MR,PC_NNN,PC_SPR))
      ALLOCATE(DTT(NH,MR,PC_NNN))
C ----- generate mean camber surface geometry
      DO J=1,MRP
        DO I=1,NHP
          BXX(I,J)=XBM(I,J)
          BYY(I,J)=YBM(I,J)
          BZZ(I,J)=ZBM(I,J)
          BRR(I,J)=SQRT(YBM(I,J)**2+ZBM(I,J)**2)
        END DO
      END DO

C ---- find the relative location between two sides of the blade
      K1=INDEXB(NH+5,1)
      K2=INDEXB(NH-5+1,1)
      IT0=0
      IF ((ATAN2(XCT(K1,3),XCT(K1,2))).GE.(ATAN2(XCT(K2,3),XCT(K2,2))))
     &      IT0=1
      WRITE(*,*) 'In PC2NS, IT0=',IT0
C ----- calculate ind vel (outer), body force
      DO K=1,PC_SPR
        CALL STRENGTHREN(K) ! RENEW POT AND DPDN VALUE
        TSTEP=DELTAT*REAL(K-1)
        IDXREV=K
        CALL INFLOW
        CALL PRSDIF(1)

        CALL INDV_BLD(INDX,INDR,INDT)
        DO I=1,NH
          DO J=1,MR
            IF (IT0.EQ.1) THEN
              K1=K+PC_SPR-NINT(DELK/2.0/PI*REAL(PC_SPR))
              IF(K1.GT.PC_SPR) K1=K1-PC_SPR
              K2=K
            ELSE
              K1=K
              K2=K+PC_SPR-NINT(DELK/2.0/PI*REAL(PC_SPR))
              IF(K2.GT.PC_SPR) K2=K2-PC_SPR
            END IF
            PC_IVX(I,J,1,K1)=INDX(NH+I,J)
            PC_IVR(I,J,1,K1)=INDR(NH+I,J)
            PC_IVT(I,J,1,K1)=INDT(NH+I,J)
            PC_IVX(I,J,PC_NNN,K2)=INDX(NH-I+1,J)
            PC_IVR(I,J,PC_NNN,K2)=INDR(NH-I+1,J)
            PC_IVT(I,J,PC_NNN,K2)=INDT(NH-I+1,J)
          END DO
        END DO ! not correct

        CALL BFCALC(BFX,BFR,BFT,SRCS)
        DO I=1,NH
          DO J=1,MR
            PC_BFX(I,J,K)=BFX(I,J)
            PC_BFR(I,J,K)=BFR(I,J)
            PC_BFT(I,J,K)=BFT(I,J)
            PC_SRCS(I,J,K)=SRCS(I,J)
          END DO
        END DO
      END DO

C ----- calculate induced velocity
      CALL INDVCALC(NC,MR,NPANZ,XCT,DELK,DTT)

C ----- generating plots
      IF (PC_SPR.EQ.1) THEN
        CALL PLT_INV_OP5
        CALL PLT_BF_OP5(BXX,BRR)
        CALL PLT_CAM_OP5(NHP,MRP,BXX,BRR)
      ELSE
        CALL PLT_INV_OP6(DTT) ! induced velocity
        CALL PLT_BF_OP6(BXX,BYY,BZZ) !body force plot
        CALL PLT_CAM_OP6(NHP,MRP,BXX,BYY,BZZ,BRR)   !camber geometry plot
      END IF

      END SUBROUTINE


      SUBROUTINE INDVCALC(NC,MR,NPANZ,XCT,DELK,DTT)
      USE M_PC2NS
      IMPLICIT NONE

      REAL DELK,XCT(NPANZ,3)

      REAL,DIMENSION(NC/2,MR,PC_NNN):: DXX,DRR,DTT
      REAL,DIMENSION((NC/2)*MR*(PC_NNN-2),3)::CRDITF
      REAL,DIMENSION((NC/2)*MR*(PC_NNN-2),3,PC_SPR)::IDVITF

      REAL XX,RR,TT,YY,ZZ,VVX,VVR,VVT,VVY,VVZ,XX2,RR2,TT2

      INTEGER I,J,K,K1,K2,NC,MR,M,N,L,NPANZ,INDEXB,NH,NHP,NCP,MRP,IT

      NH=NC/2
      NHP=NH+1
      NCP=NC/2+1
      MRP=MR+1

C ---- Calculate the coordinates
      K1=INDEXB(NH+5,1)
      K2=INDEXB(NH-5+1,1)
      TT=ATAN2(XCT(K1,3),XCT(K1,2))
      TT2=ATAN2(XCT(K2,3),XCT(K2,2))
      IF ((ABS(TT-TT2)).GT.(3.0)) THEN
        WRITE(*,*) 'WARNING:(indvcalc) wrong tangential coordinates!'
      END IF
      IT=0
      IF (TT.GE.TT2) IT=1


      DO I=1,NH
        DO J=1,MR
          K1=INDEXB(NH+I,J)
          K2=INDEXB(NH-I+1,J)
          XX=XCT(K1,1)
          RR=SQRT(XCT(K1,2)**2+XCT(K1,3)**2)
          TT=ATAN2(XCT(K1,3),XCT(K1,2))
          XX2=XCT(K2,1)
          RR2=SQRT(XCT(K2,2)**2+XCT(K2,3)**2)
          TT2=ATAN2(XCT(K2,3),XCT(K2,2))
          IF ((ABS(TT-TT2)).GT.(3.0)) THEN
            WRITE(*,*) 'WARNING:(indvcalc) wrong tangential coordinates!'
          END IF

          IF (IT.EQ.1) THEN
            TT=TT-DELK
          ELSE
            TT2=TT2-DELK
          END IF

          DO K=1,PC_NNN
            DXX(I,J,K)=(XX*REAL(PC_NNN-K)+XX2*REAL(K-1))/REAL(PC_NNN-1)
            DRR(I,J,K)=(RR*REAL(PC_NNN-K)+RR2*REAL(K-1))/REAL(PC_NNN-1)
            DTT(I,J,K)=(TT*REAL(PC_NNN-K)+TT2*REAL(K-1))/REAL(PC_NNN-1)
          END DO

        END DO
      END DO

C --- CALCULATE INNER INDUCED VELOCITY
      L=0
      DO I=1,NH
        DO J=1,MR
          DO K=2,(PC_NNN-1)
            L=L+1
            CRDITF(L,1)=DXX(I,J,K)
            CRDITF(L,2)=DRR(I,J,K)
            CRDITF(L,3)=DTT(I,J,K)
          END DO
        END DO
      END DO
      CALL INDV_INNER(CRDITF,IDVITF,(PC_NNN-2)*NH*MR,PC_SPR)
      L=0
      DO I=1,NH
        DO J=1,MR
          DO K=2,(PC_NNN-1)
            L=L+1
            DO IT=1,PC_SPR
              PC_IVX(I,J,K,IT)=IDVITF(L,1,IT)
              PC_IVR(I,J,K,IT)=IDVITF(L,2,IT)
              PC_IVT(I,J,K,IT)=IDVITF(L,3,IT)
            END DO
          END DO
        END DO
      END DO
C ---- Plot 3D induced velocity
      OPEN(725,FILE='Test_Indv_3d.plt', STATUS='UNKNOWN')
      WRITE(725,*) 'VARIABLES="X","Y","Z","VX","VY","VZ","VR","VT"'
      DO IT=1,PC_SPR
        DO I=1,NH
          WRITE(725,*) 'ZONE T="SLICE',I,'",SOLUTIONTIME=',REAL(IT),',I='
     *                  ,MR,',J=',PC_NNN,',K=',1
          DO K=1,PC_NNN
            DO J=1,MR
              XX=DXX(I,J,K)
              RR=DRR(I,J,K)
              TT=DTT(I,J,K)
              YY=RR*COS(TT)
              ZZ=RR*SIN(TT)
              VVX=PC_IVX(I,J,K,IT)
              VVR=PC_IVR(I,J,K,IT)
              VVT=PC_IVT(I,J,K,IT)
              VVY=VVR*COS(TT)-VVT*SIN(TT)
              VVZ=VVR*SIN(TT)+VVT*COS(TT)
              WRITE(725,*) XX,YY,ZZ,VVX,VVY,VVZ,VVR,VVT
            END DO
          END DO
        END DO
      END DO
      CLOSE(725)
      RETURN
      END SUBROUTINE
      SUBROUTINE INDV_INNER(CRDITF,IDVITF,RYTOT,PC_SPR)
C ===========================================================
C This subroutine calculates the induced velocities.
C Currently working on!!
C ===========================================================
      USE OMP_LIB
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INTEGER::RYTOT,PC_SPR,PCNTPOS
      REAL,DIMENSION(RYTOT,3)::CRDITF
      REAL,DIMENSION(RYTOT,3,PC_SPR)::IDVITF,IDVITF1
      REAL::SSX,SSY,SSZ,SSR,SST
      REAL::PV_XV(4),PV_YV(4),PV_SIDE(4),PV_S(15)
      INTEGER::NNNN,I,J,K,L,II,JJ,IT
      REAL::TB_POT(NBLADE,NPANEL,PC_SPR),TB_DPDN(NBLADE,NPANEL,PC_SPR)
      REAL::TB_TEMP5(NBLADE,NWMIN*MR,PC_SPR)
      REAL::TB_TEMPD4(NBLADE,NDWK*MDUCT,PC_SPR)
      INTEGER,ALLOCATABLE::I_A(:,:),I2_A(:,:)
      INTEGER::om_ct1,om_ct2,om_ctr
      REAL::om_dt

C - Induced velocity set to zero
      DO J=1,RYTOT
        DO L=1,3
          DO IT=1,PC_SPR
            IDVITF(J,L,IT)=0.0E0
            IDVITF1(J,L,IT)=0.0E0
          END DO
        END DO
      END DO

C - Read source/dipole strength
      NREAD = NWMIN*MR
      NREADD = NDWK*MDUCT
      DO IT=1,PC_SPR
        DO K=1,NBLADE
          IREC = PCNTPOS(K,IT)
          CALL READ2(45,IREC,POT,NPANEL)
          CALL READ2(47,IREC,DPDN,NPANEL)
          DO I=1,NPANEL
            TB_POT(K,I,IT)=POT(I)
            TB_DPDN(K,I,IT)=DPDN(I)
          END DO

          CALL READ2(46,IREC,TEMP5,NREAD)
          DO I=1,NREAD
            TB_TEMP5(K,I,IT)=TEMP5(I)
          END DO

          IF (IDUCT.NE.0) THEN
            CALL READ2(49,IREC,TEMPD4,NREADD)
            DO I=1,NREADD
              TB_TEMPD4(K,I,IT)=TEMPD4(I)
            END DO
          END IF
        END DO
      END DO

!$    om_ctr = omp_get_max_threads()
!$    write(*,*) 'Number of threads is: ',om_ctr
      call system_clock(om_ct1,om_ctr)

C =============================================================================
C                    Begin wall panel induced velocity!
C =============================================================================

      WRITE(*,*) 'Begin calculating blade/hub/casing induced velocity! '

cS.KIM -- The parallel part below is NOT done correctly, so I turned it off.
c!$OMP PARALLEL PRIVATE(I,J,II,PV_XV,PV_YV,PV_SIDE,PV_S,SSX,SSY,SSZ,SSR,SST,
c!$OMP&                 K,XLOC,YLOC,ZLOC,FS,FD,FSX,FSY,FDX,FDY,FDZ,VDX,VDY,
c!$OMP&                 VDZ,VSX,VSY,VSZ,VX0,VY0,VZ0,VR0,VT0,IT)
c!$OMP& SHARED (NPANEL,RYTOT,XVP,YVP,SID,SS,NBLADE,CRDITF,DELK,DIR,XCT,
c!$OMP&         CHRLEPS,TB_POT,TB_DPDN,PI,PC_SPR)
c!$OMP& REDUCTION (+:IDVITF1)
c!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC,5000)
      DO I=1,NPANEL
        DO J=1,RYTOT

C -- panel geometric parameter
          DO II=1,4
            PV_XV(II) = XVP(I,II)
            PV_YV(II) = YVP(I,II)
            PV_SIDE(II) = SID(I,II)
          END DO
          DO II=1,15
            PV_S(II) = SS(I,II)
          END DO

          SSX=CRDITF(J,1)
          SSR=CRDITF(J,2)

          DO K=1,NBLADE
C - just in case, ISTEADY=0 for current runs

            SST=CRDITF(J,3)+DELK*REAL(K-1)
            SSY=SSR*COS(SST)
            SSZ=SSR*SIN(SST)

            XLOC=(SSX-XCT(I,1))*DIR(I,1,1)+(SSY-XCT(I,2))*DIR(I,1,2)
     *            +(SSZ-XCT(I,3))*DIR(I,1,3)
            YLOC=(SSX-XCT(I,1))*DIR(I,2,1)+(SSY-XCT(I,2))*DIR(I,2,2)
     *            +(SSZ-XCT(I,3))*DIR(I,2,3)
            ZLOC=(SSX-XCT(I,1))*DIR(I,3,1)+(SSY-XCT(I,2))*DIR(I,3,2)
     *            +(SSZ-XCT(I,3))*DIR(I,3,3)

            CALL RPAN_PV(XLOC,YLOC,ZLOC,CHRLEPS(I),FS,FD,FSX,FSY,
     *                 FDX,FDY,FDZ,1,1,PV_XV,PV_YV,PV_S,PV_SIDE)

            IF(ABS(FD).GT.6.28) THEN
              FD=0.0
            END IF

            VDX = FDX*DIR(I,1,1)+FDY*DIR(I,2,1)+FDZ*DIR(I,3,1)
            VDY = FDX*DIR(I,1,2)+FDY*DIR(I,2,2)+FDZ*DIR(I,3,2)
            VDZ = FDX*DIR(I,1,3)+FDY*DIR(I,2,3)+FDZ*DIR(I,3,3)
            VSX = FSX*DIR(I,1,1)+FSY*DIR(I,2,1)-FD*DIR(I,3,1)
            VSY = FSX*DIR(I,1,2)+FSY*DIR(I,2,2)-FD*DIR(I,3,2)
            VSZ = FSX*DIR(I,1,3)+FSY*DIR(I,2,3)-FD*DIR(I,3,3)

            DO IT=1,PC_SPR
              VX0 = -1.0*(TB_POT(K,I,IT)*VDX-TB_DPDN(K,I,IT)*VSX)/4.0/PI
              VY0 = -1.0*(TB_POT(K,I,IT)*VDY-TB_DPDN(K,I,IT)*VSY)/4.0/PI
              VZ0 = -1.0*(TB_POT(K,I,IT)*VDZ-TB_DPDN(K,I,IT)*VSZ)/4.0/PI
              VR0 = VZ0*SIN(SST)+VY0*COS(SST)
              VT0 = VZ0*COS(SST)-VY0*SIN(SST)

              IDVITF1(J,1,IT) = IDVITF1(J,1,IT)+VX0
              IDVITF1(J,2,IT) = IDVITF1(J,2,IT)+VR0
              IDVITF1(J,3,IT) = IDVITF1(J,3,IT)+VT0
            END DO
          END DO
        END DO
      END DO
c!$OMP END DO
c!$OMP END PARALLEL
      DO J=1,RYTOT
        DO L=1,3
          DO IT=1,PC_SPR
            IDVITF(J,L,IT)=IDVITF(J,L,IT)+IDVITF1(J,L,IT)
            IDVITF1(J,L,IT)=0.0E0
          END DO
        END DO
      END DO



C =============================================================================
C                    Begin wake induced velocity!
C =============================================================================

      WRITE(*,*) 'Begin calculating wake induced velocity! '

      ALLOCATE (I_A(MR,1000),I2_A(MR,1000))
      DO M=1,MR
        DO N=1,NWPAN(M)
          I_A(M,N)=IDXWAK(N,M)
          I2_A(M,N)=INDEXW2(N,M)
        END DO
      END DO
c!$OMP PARALLEL PRIVATE(J,K,M,N,I,I2,II,SSX,SSR,SST,SSY,SSZ,XLOC,YLOC,ZLOC,FS,FD,
c!$OMP&                 FSX,FSY,FDX,FDY,FDZ,PV_XV,PV_YV,PV_S,PV_SIDE,
c!$OMP&                 VDX,VDY,VDZ,VX0,VY0,VZ0,VR0,VT0,IT)
c!$OMP& SHARED(NBLADE,RYTOT,CRDITF,DELK,DIRW,XCTW,CHRLEPS,TB_TEMP5,PI,I_A,
c!$OMP&        I2_A,XVPW,YVPW,SIDW,SSW,MR,NWPAN,NREAD,TB_T5,PC_SPR)
c!$OMP& REDUCTION (+:IDVITF1)
c!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC,5000)
      DO J=1,RYTOT
        DO K=1,NBLADE
          DO M=1,MR
            DO N=1,NWPAN(M)
              I= I_A(M,N)
              I2 = I2_A(M,N)
              DO II=1,4
                PV_XV(II) = XVPW(I,II)
                PV_YV(II) = YVPW(I,II)
                PV_SIDE(II) = SIDW(I,II)
              END DO
              DO II=1,15
                PV_S(II) = SSW(I,II)
              END DO

              SSX=CRDITF(J,1)
              SSR=CRDITF(J,2)
              SST=CRDITF(J,3)+DELK*REAL(K-1)
              SSY=SSR*COS(SST)
              SSZ=SSR*SIN(SST)

              XLOC=(SSX-XCTW(I,1))*DIRW(I,1,1)+(SSY-XCTW(I,2))
     *              *DIRW(I,1,2)+(SSZ-XCTW(I,3))*DIRW(I,1,3)
              YLOC=(SSX-XCTW(I,1))*DIRW(I,2,1)+(SSY-XCTW(I,2))
     *              *DIRW(I,2,2)+(SSZ-XCTW(I,3))*DIRW(I,2,3)
              ZLOC=(SSX-XCTW(I,1))*DIRW(I,3,1)+(SSY-XCTW(I,2))
     *              *DIRW(I,3,2)+(SSZ-XCTW(I,3))*DIRW(I,3,3)

              CALL RPAN_PV(XLOC,YLOC,ZLOC,CHRLEPS(I),FS,FD,FSX,FSY,
     *                 FDX,FDY,FDZ,1,1,PV_XV,PV_YV,PV_S,PV_SIDE)

              IF(ABS(FD).GT.6.28) THEN
                FD=0.0
              END IF

              VDX = FDX*DIRW(I,1,1)+FDY*DIRW(I,2,1)+FDZ*DIRW(I,3,1)
              VDY = FDX*DIRW(I,1,2)+FDY*DIRW(I,2,2)+FDZ*DIRW(I,3,2)
              VDZ = FDX*DIRW(I,1,3)+FDY*DIRW(I,2,3)+FDZ*DIRW(I,3,3)

              DO IT=1,PC_SPR
                TB_T5=TB_TEMP5(K,I2,IT)
                IF (I2.GT.NREAD) TB_T5=0.0

                VX0 = -1.0*(TB_T5*VDX)/4.0/PI
                VY0 = -1.0*(TB_T5*VDY)/4.0/PI
                VZ0 = -1.0*(TB_T5*VDZ)/4.0/PI
                VR0 = VZ0*SIN(SST)+VY0*COS(SST)
                VT0 = VZ0*COS(SST)-VY0*SIN(SST)
                IDVITF1(J,1,IT) = IDVITF1(J,1,IT)+VX0
                IDVITF1(J,2,IT) = IDVITF1(J,2,IT)+VR0
                IDVITF1(J,3,IT) = IDVITF1(J,3,IT)+VT0
              END DO
            END DO
          END DO
        END DO
      END DO
c!$OMP END DO
c!$OMP END PARALLEL
      DO J=1,RYTOT
        DO L=1,3
          DO IT=1,PC_SPR
            IDVITF(J,L,IT)=IDVITF(J,L,IT)+IDVITF1(J,L,IT)
            IDVITF1(J,L,IT)=0.0E0
          END DO
        END DO
      END DO

C =============================================================================
C                    Begin duct wake induced velocity!
C =============================================================================
      IF (IDUCT.EQ.1) THEN
        WRITE(*,*) 'Begin calculating DUCT wake induced velocity! '

        DEALLOCATE(I_A,I2_A)
        ALLOCATE (I_A(MDUCT,1000))
        DO M=1,MDUCT
          DO N=1,NDWK
            I_A(M,N)=INDEXWD(N,M)
          END DO
        END DO
c!$OMP PARALLEL PRIVATE(J,K,M,N,I,I2,SSX,SSR,SST,SSY,SSZ,XLOC,YLOC,ZLOC,FS,FD,
c!$OMP&                 FSX,FSY,FDX,FDY,FDZ,PV_XV,PV_YV,PV_S,PV_SIDE,
c!$OMP&                 VDX,VDY,VDZ,VX0,VY0,VZ0,VR0,VT0,IT)
c!$OMP& SHARED(NBLADE,RYTOT,MDUCT,NDWK,CRDITF,DELK,DIRWD,XCTWD,CHRLEPS,TB_TEMPD4,PI,I_A,
c!$OMP&        I2_A,XVPWD,YVPWD,SIDWD,SSWD,PC_SPR)
c!$OMP& REDUCTION (+:IDVITF1)
c!$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC,5000)
        DO K=1,NBLADE
          DO M=1,MDUCT
            DO N=1,NDWK
              DO J=1,RYTOT
                I=I_A(M,N)
                I2=I

                DO II=1,4
                  PV_XV(II) = XVPWD(I,II)
                  PV_YV(II) = YVPWD(I,II)
                  PV_SIDE(II) = SIDWD(I,II)
                END DO
                DO II=1,15
                  PV_S(II) = SSWD(I,II)
                END DO

                SSX=CRDITF(J,1)
                SSR=CRDITF(J,2)
                SST=CRDITF(J,3)+DELK*REAL(K-1)
                SSY=SSR*COS(SST)
                SSZ=SSR*SIN(SST)

                XLOC=(SSX-XCTWD(I,1))*DIRWD(I,1,1)+(SSY-XCTWD(I,2))
     *              *DIRWD(I,1,2)+(SSZ-XCTWD(I,3))*DIRWD(I,1,3)
                YLOC=(SSX-XCTWD(I,1))*DIRWD(I,2,1)+(SSY-XCTWD(I,2))
     *              *DIRWD(I,2,2)+(SSZ-XCTWD(I,3))*DIRWD(I,2,3)
                ZLOC=(SSX-XCTWD(I,1))*DIRWD(I,3,1)+(SSY-XCTWD(I,2))
     *              *DIRWD(I,3,2)+(SSZ-XCTWD(I,3))*DIRWD(I,3,3)

                CALL RPAN_PV(XLOC,YLOC,ZLOC,CHRLEPS(I),FS,FD,FSX,FSY,
     *                        FDX,FDY,FDZ,1,1,PV_XV,PV_YV,PV_S,PV_SIDE)

                IF(ABS(FD).GT.6.28) THEN
                  FD=0.0
                END IF

                VDX = FDX*DIRWD(I,1,1)+FDY*DIRWD(I,2,1)+FDZ*DIRWD(I,3,1)
                VDY = FDX*DIRWD(I,1,2)+FDY*DIRWD(I,2,2)+FDZ*DIRWD(I,3,2)
                VDZ = FDX*DIRWD(I,1,3)+FDY*DIRWD(I,2,3)+FDZ*DIRWD(I,3,3)
                DO IT=1,PC_SPR
                  VX0 = -1.0*(TB_TEMPD4(K,I2,IT)*VDX)/4.0/PI
                  VY0 = -1.0*(TB_TEMPD4(K,I2,IT)*VDY)/4.0/PI
                  VZ0 = -1.0*(TB_TEMPD4(K,I2,IT)*VDZ)/4.0/PI
                  VR0 = VZ0*SIN(SST)+VY0*COS(SST)
                  VT0 = VZ0*COS(SST)-VY0*SIN(SST)
                  IDVITF1(J,1,IT) = IDVITF1(J,1,IT)+VX0
                  IDVITF1(J,2,IT) = IDVITF1(J,2,IT)+VR0
                  IDVITF1(J,3,IT) = IDVITF1(J,3,IT)+VT0
                END DO
              END DO
            END DO
          END DO
        END DO
c!$OMP END DO
c!$OMP END PARALLEL
        DO J=1,RYTOT
          DO L=1,3
            DO IT=1,PC_SPR
              IDVITF(J,L,IT)=IDVITF(J,L,IT)+IDVITF1(J,L,IT)
              IDVITF1(J,L,IT)=0.0E0
            END DO
          END DO
        END DO
      END IF

C =========== finished =====================
      WRITE(*,*) 'Finished calculating wake induced velocity'
C plot time cost
      call system_clock(om_ct2)
      om_dt=(real(om_ct2)-real(om_ct1))/real(om_ctr)
      write(*,*) 'Time cost for calculating induced velocity:',om_dt

      RETURN
      END


      SUBROUTINE INDV_BLD(INDX,INDR,INDT)
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      REAL,DIMENSION(NC,MR)::INDX,INDR,INDT

      DO M=1,MR
        DO N=1,NC
          L=INDEXB(N,M)
          UXI=DPDUB(N,M)
          UETA=(DPDVB(N,M)-DPDUB(N,M)*SINPHI(L)) / COSPHI(L)
          UNORM=DPDN(L)
          UTXX=UXI*DIR(L,1,1)+UETA*DIR(L,2,1)+UNORM*DIR(L,3,1)
          UTYY=UXI*DIR(L,1,2)+UETA*DIR(L,2,2)+UNORM*DIR(L,3,2)
          UTZZ=UXI*DIR(L,1,3)+UETA*DIR(L,2,3)+UNORM*DIR(L,3,3)
          TPTT=ATAN2(XCT(L,3),XCT(L,2))
          INDX(N,M)=UTXX
          INDR(N,M)=UTZZ*SIN(TPTT)+UTYY*COS(TPTT)
          INDT(N,M)=UTZZ*COS(TPTT)-UTYY*SIN(TPTT)
        END DO
      END DO

      RETURN
      END
      SUBROUTINE BFCALC(BFX,BFR,BFT,SRCS)
      USE M_BLADEM
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      REAL,DIMENSION(NH,MR)::BFX,BFR,BFT,SRCS,BMR,TT_AREA
      REAL,DIMENSION(NHP,MRP)::BXX,BRR
      DOUBLE PRECISION,DIMENSION(3):: TT_V1,TT_V2,TT_V3

      DO J=1,MRP
        DO I=1,NHP
          BXX(I,J)=XBM(I,J)
          BRR(I,J)=SQRT(YBM(I,J)**2+ZBM(I,J)**2)
        END DO
      END DO
C --- Calculate body force on blade surface
      DO J=1,MR
        DO I = 1,NH
          I2=I+1
          J2=J+1
          TC_XX=(XBM(I,J)+XBM(I,J2)+XBM(I2,J)+XBM(I2,J2))/4.0E0
          TC_YY=(YBM(I,J)+YBM(I,J2)+YBM(I2,J)+YBM(I2,J2))/4.0E0
          TC_ZZ=(ZBM(I,J)+ZBM(I,J2)+ZBM(I2,J)+ZBM(I2,J2))/4.0E0
          TC_RR=SQRT(TC_YY**2+TC_ZZ**2)
          TC_TT=ATAN2(TC_ZZ,TC_YY)

          TT_V1(1)=DBLE((BXX(I2,J)+BXX(I2,J2)-BXX(I,J)-BXX(I,J2))/2.0E0)
          TT_V1(2)=DBLE((BRR(I2,J)+BRR(I2,J2)-BRR(I,J)-BRR(I,J2))/2.0E0)
          TT_V1(3)=0.0D0
          TT_V2(1)=DBLE((BXX(I,J2)+BXX(I2,J2)-BXX(I,J)-BXX(I2,J))/2.0E0)
          TT_V2(2)=DBLE((BRR(I,J2)+BRR(I2,J2)-BRR(I,J)-BRR(I2,J))/2.0E0)
          TT_V2(3)=0.0D0
          CALL EXPROD(TT_V1,TT_V2,TT_V3)
          TC_AREA=REAL(DSQRT(TT_V3(1)**2+TT_V3(2)**2+TT_V3(3)**2))
          TT_AREA(I,J)=TC_AREA

          L1=INDEXB(NH+I,J)
          L2=INDEXB(NHP-I,J)

          TC_PX=(CPB(NH+I,J)*VEL(L1,1)*SS(L1,1)+CPB(NHP-I,J)*VEL(L2,1)
     *            *SS(L2,1))/TC_AREA
          TC_PY=(CPB(NH+I,J)*VEL(L1,2)*SS(L1,1)+CPB(NHP-I,J)*VEL(L2,2)
     *            *SS(L2,1))/TC_AREA
          TC_PZ=(CPB(NH+I,J)*VEL(L1,3)*SS(L1,1)+CPB(NHP-I,J)*VEL(L2,3)
     *            *SS(L2,1))/TC_AREA
          TC_SC=(DPDN(L1)*SS(L1,1)+DPDN(L2)*SS(L2,1))/TC_AREA
c method 3 for calculating the thickness source term
*          TT_1=SQRT(VEL(L1,2)**2+VEL(L1,3)**2)
*          TT_2=SQRT(VEL(L2,2)**2+VEL(L2,3)**2)
*          T_SC=(DPDN(L1)*SS(L1,1)*TT_1+DPDN(L2)*SS(L2,1)*TT_2)/TC_AREA

          TC_PX=TC_PX/2.0E0/2.0/PI/TC_RR*NBLADE
          TC_PY=TC_PY/2.0E0/2.0/PI/TC_RR*NBLADE
          TC_PZ=TC_PZ/2.0E0/2.0/PI/TC_RR*NBLADE
          SRCS(I,J)=-1.0*TC_SC/2.0/PI/TC_RR*NBLADE   !1/s or kg/m^3/s since density=1
          BFX(I,J)=-1.0*TC_PX  !N/m^3
          BFR(I,J)=-1.0*TC_PY*COS(TC_TT)-TC_PZ*SIN(TC_TT)  !N/m^3
          BFT(I,J)=-1.0*TC_PZ*COS(TC_TT)+TC_PY*SIN(TC_TT)  !N/m^3

C--calculate R-moment on the blade
          TC_MR=0.0E0
          TC_TT=TC_TT*TC_RR
          L=INDEXB(NH+I,J)
          TB_XX=XCT(L,1)
          TB_TT=ATAN2(XCT(L,3),XCT(L,2))
          TB_RR=SQRT(XCT(L,3)**2+XCT(L,2)**2)
          TB_NX=VEL(L,1)
          TB_NT=VEL(L,3)*COS(T_TT)-VEL(L,2)*SIN(T_TT)
          TB_NN=SQRT(TB_NT**2+TB_NX**2)
          TB_NT=TB_NT/TB_NN
          TB_NX=TB_NX/TB_NN
          TB_TT=TB_TT*TB_RR
          TB_DD=(TB_TT-TC_TT)*TB_NX-(TB_XX-TC_XX)*TB_NT  ! left - right +
          TC_MR=TC_MR+CPB(NH+I,J)*SS(L,1)*TB_DD  ! counter-clockwise +

          L=INDEXB(NHP-I,J)
          TB_XX=XCT(L,1)
          TB_TT=ATAN2(XCT(L,3),XCT(L,2))
          TB_RR=SQRT(XCT(L,3)**2+XCT(L,2)**2)
          TB_NX=VEL(L,1)
          TB_NT=VEL(L,3)*COS(T_TT)-VEL(L,2)*SIN(T_TT)
          TB_NN=SQRT(TB_NT**2+TB_NX**2)
          TB_NT=TB_NT/TB_NN
          TB_NX=TB_NX/TB_NN
          TB_TT=TB_TT*TB_RR
          TB_DD=(TB_TT-TC_TT)*TB_NX-(TB_XX-TC_XX)*TB_NT  ! left - right +
          TC_MR=TC_MR+CPB(NHP-I,J)*SS(L,1)*TB_DD  ! counter-clockwise +

          BMR(I,J)=-1.0*TC_MR/TC_AREA*2.0E0/2.0/PI/TC_RR*NBLADE  !(N.m)/(m^3)
C-- finish R-moment calculation
        END DO
      END DO

C -- Average of the first 3 panels on the LE
      DO J=1,MR
        TC_AREA=TT_AREA(1,J)+TT_AREA(2,J)+TT_AREA(3,J)
        TC_PX=(BFX(1,J)*TT_AREA(1,J)+BFX(2,J)*TT_AREA(2,J)
     *        +BFX(3,J)*TT_AREA(3,J))/TC_AREA
        TC_PR=(BFR(1,J)*TT_AREA(1,J)+BFR(2,J)*TT_AREA(2,J)
     *        +BFR(3,J)*TT_AREA(3,J))/TC_AREA
        TC_PT=(BFT(1,J)*TT_AREA(1,J)+BFT(2,J)*TT_AREA(2,J)
     *        +BFT(3,J)*TT_AREA(3,J))/TC_AREA
        DO I = 1,3
          BFX(I,J)=TC_PX
          BFR(I,J)=TC_PR
          BFT(I,J)=TC_PT
        END DO
      END DO

      RETURN
      END SUBROUTINE


      SUBROUTINE PLT_BF_OP5(BXX,BRR)
      USE M_PC2NS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      REAL,DIMENSION(NHP,MRP)::BXX,BRR

      REAL RVT
      LOGICAL IRVT
      RVT=1.0
      INQUIRE(FILE='rvt.inp',EXIST=IRVT)
      IF (IRVT.EQ.(.TRUE.)) THEN
        RVT=-1.0
        WRITE(*,*) 'Left handed propeller in calculating body force!'
      END IF

C ... Test plot curvilinear coordinates and body force
      OPEN(785,FILE='test_realcamber_bf.plt',STATUS='UNKNOWN')
      WRITE(785,*) 'VARIABLES="X","R","BFX","BFR","BFT","SOURCE"'
      DO IT=1,PC_SPR
        WRITE(785,*) 'ZONE T="REAL_cam_bf",SOLUTIONTIME=',REAL(IT),',I=',NHP,',J=',MRP,',K=1,'
        WRITE(785,*) 'DATAPACKING=BLOCK, VARLOCATION=(3=CELLCENTERED,4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED)'
        WRITE(785,*) ((BXX(I,J),I=1,NHP),J=1,MRP)
        WRITE(785,*) ((BRR(I,J),I=1,NHP),J=1,MRP)
        WRITE(785,*) ((PC_BFX(I,J,IT),I=1,NH),J=1,MR)
        WRITE(785,*) ((PC_BFR(I,J,IT),I=1,NH),J=1,MR)
        WRITE(785,*) ((RVT*PC_BFT(I,J,IT),I=1,NH),J=1,MR)
        WRITE(785,*) ((PC_SRCS(I,J,IT),I=1,NH),J=1,MR)
      END DO
      CLOSE(785)

      OPEN(714,FILE='prop_bf.ip', STATUS='UNKNOWN')
      WRITE(714,*) "2"
      WRITE(714,*) "2"
      WRITE(714,*) MR*NH
      WRITE(714,*) "4"
      WRITE(714,*) "uds-0"
      WRITE(714,*) "uds-1"
      WRITE(714,*) "uds-2"
      WRITE(714,*) "uds-3"
C --- Plot file prop_bf.ip and AREA.dat
103   FORMAT(F16.6)
      DO J=1,MR
        DO I=1,NH
        XX=(BXX(I,J)+BXX(I+1,J)+BXX(I+1,J+1)+BXX(I,J+1))
        XX=XX/4.0
        WRITE(714,103) XX
        END DO
      END DO
      DO J=1,MR
        DO I=1,NH
        RR=(BRR(I,J)+BRR(I+1,J)+BRR(I+1,J+1)+BRR(I,J+1))
        RR=RR/4.0
        WRITE(714,103) RR
        END DO
      END DO
      WRITE(714,103) ((PC_BFX(N,M,1),N=1,NH),M=1,MR)
      WRITE(714,103) ((PC_BFR(N,M,1),N=1,NH),M=1,MR)
      WRITE(714,103) ((RVT*PC_BFT(N,M,1),N=1,NH),M=1,MR)
      WRITE(714,103) ((PC_SRCS(N,M,1),N=1,NH),M=1,MR)
      CLOSE(714)
      END SUBROUTINE
      SUBROUTINE PLT_BF_OP6(BXX,BYY,BZZ)
      USE M_PC2NS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      REAL,DIMENSION(NHP,MRP)::BXX,BYY,BZZ
      REAL,DIMENSION(NH,MR,PC_SPR)::PCX,PCY,PCZ,PFX,PFY,PFZ,PSC

      REAL RVT
      LOGICAL IRVT
      RVT=1.0
      INQUIRE(FILE='rvt.inp',EXIST=IRVT)
      IF (IRVT.EQ.(.TRUE.)) THEN
        RVT=-1.0
        WRITE(*,*) 'Left handed propeller in calculating body force!'
      END IF

C ... Test plot curvilinear coordinates and body force
      OPEN(785,FILE='Test_BF_blade.plt',STATUS='UNKNOWN')
      WRITE(785,*) 'VARIABLES="X","Y","Z","BFX","BFR","BFT","SOURCE"'
      DO IT=1,PC_SPR
        WRITE(785,*) 'ZONE T="REAL_cam_bf",SOLUTIONTIME=',REAL(IT),',I=',NHP,',J=',MRP,',K=1,'
        WRITE(785,*) 'DATAPACKING=BLOCK, VARLOCATION=(4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED,7=CELLCENTERED)'
        WRITE(785,*) ((BXX(I,J),I=1,NHP),J=1,MRP)
        WRITE(785,*) ((BYY(I,J),I=1,NHP),J=1,MRP)
        WRITE(785,*) ((BZZ(I,J),I=1,NHP),J=1,MRP)
        WRITE(785,*) ((PC_BFX(I,J,IT),I=1,NH),J=1,MR)
        WRITE(785,*) ((PC_BFR(I,J,IT),I=1,NH),J=1,MR)
        WRITE(785,*) ((RVT*PC_BFT(I,J,IT),I=1,NH),J=1,MR)
        WRITE(785,*) ((PC_SRCS(I,J,IT),I=1,NH),J=1,MR)
      END DO
      CLOSE(785)

      OPEN(777,FILE='Test_BF_disk.plt',STATUS='UNKNOWN')
      WRITE(777,*) 'VARIABLES="X","Y","Z","BFX","BFR","BFT","SRCS"'
      WRITE(777,*) 'ZONE T = "BLD",I=',NH,',J=',MR,',K=',NTPREV+1
      DO K1=1,NTPREV+1
        K=K1
        IF (K.EQ.(NTPREV+1)) K=1
        DO J=1,MR
          DO I=1,NH
            L=INDEXB(NH+I,J)
            TT_RR=SQRT(XCT(L,2)**2+XCT(L,3)**2)
            TT_TT=ATAN2(XCT(L,3),XCT(L,2))-REAL(K-1)*2.0*PI/REAL(NTPREV)
            XX=XCT(L,1)
            YY=TT_RR*COS(TT_TT)
            ZZ=TT_RR*SIN(TT_TT)

            VVX=PC_BFX(I,J,K)
            VVR=PC_BFR(I,J,K)
            VVT=RVT*PC_BFT(I,J,K)
            VVY=PC_SRCS(I,J,K)

            WRITE(777,*) XX,YY,ZZ,VVX,VVR,VVT,VVY
          END DO
        END DO
      END DO
      CLOSE(777)

      OPEN(714,FILE='prop_bf.ip', STATUS='UNKNOWN')
      WRITE(714,*) "2"
      WRITE(714,*) "3"
      WRITE(714,*) MR*NH*PC_SPR
      WRITE(714,*) "4"
      WRITE(714,*) "uds-0"
      WRITE(714,*) "uds-1"
      WRITE(714,*) "uds-2"
      WRITE(714,*) "uds-3"
C --- Plot file prop_bf.ip and AREA.dat
103   FORMAT(F16.6)

      DO I=1,NH
        DO J=1,MR
          DO K=1,PC_SPR
            PCX(I,J,K)=(BXX(I,J)+BXX(I+1,J)+BXX(I+1,J+1)+BXX(I,J+1))/4.0
            YY=(BYY(I,J)+BYY(I+1,J)+BYY(I+1,J+1)+BYY(I,J+1))/4.0
            ZZ=(BZZ(I,J)+BZZ(I+1,J)+BZZ(I+1,J+1)+BZZ(I,J+1))/4.0
            TT_RR=SQRT(YY**2+ZZ**2)
            TT_TT=ATAN2(ZZ,YY)+REAL(K-1)*2.0*PI/REAL(NTPREV)
            PCY(I,J,K)=TT_RR*COS(TT_TT)
            PCZ(I,J,K)=TT_RR*SIN(TT_TT)*RVT
            L=NTPREV+2-K
            IF (K.EQ.1) L=1
            PSC(I,J,K)=PC_SRCS(I,J,L)
            PFX(I,J,K)=PC_BFX(I,J,L)
            PFY(I,J,K)=PC_BFR(I,J,L)*COS(TT_TT)
     &                    -PC_BFT(I,J,L)*SIN(TT_TT)
            PFZ(I,J,K)=RVT*PC_BFR(I,J,L)*SIN(TT_TT)
     &                    +RVT*PC_BFT(I,J,L)*COS(TT_TT)
          END DO
        END DO
      END DO

      WRITE(714,103) (((PCX(I,J,K),I=1,NH),J=1,MR),K=1,PC_SPR)
      WRITE(714,103) (((PCY(I,J,K),I=1,NH),J=1,MR),K=1,PC_SPR)
      WRITE(714,103) (((PCZ(I,J,K),I=1,NH),J=1,MR),K=1,PC_SPR)
      WRITE(714,103) (((PFX(I,J,K),I=1,NH),J=1,MR),K=1,PC_SPR)
      WRITE(714,103) (((PFY(I,J,K),I=1,NH),J=1,MR),K=1,PC_SPR)
      WRITE(714,103) (((PFZ(I,J,K),I=1,NH),J=1,MR),K=1,PC_SPR)
      WRITE(714,103) (((PSC(I,J,K),I=1,NH),J=1,MR),K=1,PC_SPR)

      CLOSE(714)
      END SUBROUTINE
      SUBROUTINE PLT_CAM_OP5(NCAM,MCAM,BXX,BRR)
      IMPLICIT NONE
      INTEGER NCAM,MCAM,N,M,I,J
      REAL BXX(NCAM,MCAM),BRR(NCAM,MCAM)
      REAL AREA(NCAM-1,MCAM-1),X(5),Y(5),A
      DO N = 1,(NCAM-1)
        DO M = 1,(MCAM-1)
          X(1)=BXX(N,M)
          Y(1)=BRR(N,M)
          X(2)=BXX(N+1,M)
          Y(2)=BRR(N+1,M)
          X(3)=BXX(N+1,M+1)
          Y(3)=BRR(N+1,M+1)
          X(4)=BXX(N,M+1)
          Y(4)=BRR(N,M+1)
          X(5)=BXX(N,M)
          Y(5)=BRR(N,M)
          CALL POLYGON(5,X,Y,A)
          AREA(N,M)=A
        END DO
      END DO

C AREA.dat file
      OPEN(710,FILE='AREA.dat',STATUS='UNKNOWN')
      DO J=1,(MCAM-1)
        DO I=1,(NCAM-1)
          WRITE(710,*) J,I,AREA(I,J),1.0
        END DO
      END DO
      CLOSE(710)

c ... test plot camber panel area
      OPEN(715,FILE='test_area.plt',STATUS='UNKNOWN')
      WRITE(715,*) 'VARIABLES="X","R","Area"'
      WRITE(715,*) 'ZONE T="camArea",I=',NCAM,',J=',MCAM,',K=1,'
      WRITE(715,*) 'DATAPACKING=BLOCK, VARLOCATION=(3=CELLCENTERED)'
      WRITE(715,*) ((BXX(I,J),I=1,NCAM),J=1,MCAM)
      WRITE(715,*) ((BRR(I,J),I=1,NCAM),J=1,MCAM)
      WRITE(715,*) ((AREA(I,J),I=1,NCAM-1),J=1,MCAM-1)
      CLOSE(715)

C ---- Generate camber surface and plot it
      OPEN(711,FILE='geopanel.dat', STATUS='UNKNOWN')
      WRITE(711,*) 'VARIABLES="X","Y","Z"'
      WRITE(711,*) 'ZONE T="panel",I=',NCAM,',J=',MCAM,',F=POINT'
      WRITE(711,'(1A,5I,5I)') '#',NCAM,MCAM
      DO J=1,MCAM
        DO I=1,NCAM
          WRITE(711,*) BXX(I,J),BRR(I,J),0.0
        END DO
      END DO
      CLOSE(711)

      RETURN
      END SUBROUTINE
      SUBROUTINE PLT_CAM_OP6(NCAM,MCAM,BXX,BYY,BZZ,BRR)
      IMPLICIT NONE
      INTEGER NCAM,MCAM,N,M,I,J
      REAL,DIMENSION(NCAM,MCAM):: BXX,BYY,BZZ,BRR
      REAL AREA(NCAM-1,MCAM-1),X(5),Y(5),A

      REAL RVT
      LOGICAL IRVT
      RVT=1.0
      INQUIRE(FILE='rvt.inp',EXIST=IRVT)
      IF (IRVT.EQ.(.TRUE.)) THEN
        RVT=-1.0
      END IF

      DO N = 1,(NCAM-1)
        DO M = 1,(MCAM-1)
          X(1)=BXX(N,M)
          Y(1)=BRR(N,M)
          X(2)=BXX(N+1,M)
          Y(2)=BRR(N+1,M)
          X(3)=BXX(N+1,M+1)
          Y(3)=BRR(N+1,M+1)
          X(4)=BXX(N,M+1)
          Y(4)=BRR(N,M+1)
          X(5)=BXX(N,M)
          Y(5)=BRR(N,M)
          CALL POLYGON(5,X,Y,A)
          AREA(N,M)=ABS(A)
        END DO
      END DO

C AREA.dat file
      OPEN(710,FILE='AREA.dat',STATUS='UNKNOWN')
      DO J=1,(MCAM-1)
        DO I=1,(NCAM-1)
          WRITE(710,*) J,I,AREA(I,J),1.0
        END DO
      END DO
      CLOSE(710)

c ... test plot camber panel area
      OPEN(715,FILE='test_area.plt',STATUS='UNKNOWN')
      WRITE(715,*) 'VARIABLES="X","Y","Z","Area"'
      WRITE(715,*) 'ZONE T="camArea",I=',NCAM,',J=',MCAM,',K=1,'
      WRITE(715,*) 'DATAPACKING=BLOCK, VARLOCATION=(4=CELLCENTERED)'
      WRITE(715,*) ((BXX(I,J),I=1,NCAM),J=1,MCAM)
      WRITE(715,*) ((BYY(I,J),I=1,NCAM),J=1,MCAM)
      WRITE(715,*) ((BZZ(I,J),I=1,NCAM),J=1,MCAM)
      WRITE(715,*) ((AREA(I,J),I=1,NCAM-1),J=1,MCAM-1)
      CLOSE(715)

C ---- Generate camber surface and plot it
      OPEN(711,FILE='geopanel.dat', STATUS='UNKNOWN')
      WRITE(711,*) 'VARIABLES="X","Y","Z"'
      WRITE(711,*) 'ZONE T="panel",I=',NCAM,',J=',MCAM,',F=POINT'
      WRITE(711,'(1A,5I,5I)') '#',NCAM,MCAM
      DO J=1,MCAM
        DO I=1,NCAM
          WRITE(711,*) BXX(I,J),BYY(I,J),BZZ(I,J)*RVT
        END DO
      END DO
      CLOSE(711)

      RETURN
      END SUBROUTINE
      SUBROUTINE PLT_INV_OP5
      USE M_PC2NS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      REAL,DIMENSION(NH,MR):: AVX,AVR,AVT

      REAL RVT
      LOGICAL IRVT
      RVT=1.0
      INQUIRE(FILE='rvt.inp',EXIST=IRVT)
      IF (IRVT.EQ.(.TRUE.)) THEN
        RVT=-1.0
        WRITE(*,*) 'Left handed propeller induced velocity!'
      END IF

C --- GET CIRCUMFENTIALLY AVERAGED VEL
      DO I=1,NH
        DO J=1,MR
          VVX=0.0E0
          VVR=0.0E0
          VVT=0.0E0
          DO K=1,PC_NNN
            VVX=VVX+PC_IVX(I,J,K,1)
            VVR=VVR+PC_IVR(I,J,K,1)
            VVT=VVT+PC_IVT(I,J,K,1)
          END DO
          VVX=VVX*2.0E0-PC_IVX(I,J,1,1)-PC_IVX(I,J,PC_NNN,1)
          VVR=VVR*2.0E0-PC_IVR(I,J,1,1)-PC_IVR(I,J,PC_NNN,1)
          VVT=VVT*2.0E0-PC_IVT(I,J,1,1)-PC_IVT(I,J,PC_NNN,1)
          AVX(I,J)=VVX/2.0E0/REAL(PC_NNN-1)!*TT-1.0*(1.0-TT)
          AVR(I,J)=VVR/2.0E0/REAL(PC_NNN-1)!*TT
          AVT(I,J)=VVT/2.0E0/REAL(PC_NNN-1)!*TT-3.1416/0.889*RR*(1.0-TT)
        END DO
      END DO

C --- Interpolate the perturbation velocity back to blade control point
      OPEN(729,FILE='prop_inv.plt',STATUS='UNKNOWN')
      WRITE(729,*) 'VARIABLES = "X","R","T","VX","VR","VT"'
      WRITE(729,*) 'ZONE T = "BLD",I=',NC,',J=',MR,',K=',1
      OPEN(728,FILE='Test_indv.plt',STATUS='UNKNOWN')
      WRITE(728,*) 'VARIABLES="X","Y","Z","VX","VR","VT"'
      WRITE(728,*) 'ZONE T = "BLD",I=',NC,',J=',MR,',K=',1
      WRITE(729,*) NC,MR,1
      DO J=1,MR
        DO I=1,NC
          L=INDEXB(I,J)
          XX=XCT(L,1)
          YY=XCT(L,2)
          ZZ=XCT(L,3)
          RR=SQRT(YY**2+ZZ**2)
          IF (I.GT.NH) THEN
            K=I-NH
          ELSE
            K=NH-I+1
          END IF
          VVX=AVX(K,J)
          VVR=AVR(K,J)
          VVT=RVT*AVT(K,J)
          WRITE(729,*) XX,RR,0.0,VVX,VVR,VVT
          WRITE(728,*) XX,YY,ZZ,VVX,VVR,VVT
        END DO
      END DO
      CLOSE(729)
      CLOSE(728)

      OPEN(713,FILE='prop_EFF_XYZ.dat', STATUS='UNKNOWN')


      WRITE(713,*) NC*MR,1
      DO M=1,MR
        DO N=1,NC
          L=INDEXB(N,M)
          TT_XX=XCT(L,1)
          TT_RR=SQRT(XCT(L,2)**2+XCT(L,3)**2)
C Yiran Su 20161018 evaluate total velocity on the camber panel control point
          L=INDEXB(NC+1-N,M)
          TT_XX=(TT_XX+XCT(L,1))*0.5
          TT_RR=(TT_RR+SQRT(XCT(L,2)**2+XCT(L,3)**2))*0.5
C Finished
          WRITE(713,*) TT_XX,TT_RR,0.0D0
        END DO
      END DO
      CLOSE(713)

      END SUBROUTINE
      SUBROUTINE PLT_INV_OP6(DTT)
      USE M_PC2NS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      REAL,DIMENSION(NC,MR,NTPREV):: AVX,AVR,AVT
      REAL,DIMENSION(NH,MR,PC_NNN):: DTT

      REAL RVT
      LOGICAL IRVT
      RVT=1.0
      INQUIRE(FILE='rvt.inp',EXIST=IRVT)
      IF (IRVT.EQ.(.TRUE.)) THEN
        RVT=-1.0
        WRITE(*,*) 'Left handed propeller induced velocity!'
      END IF

*      WRITE(*,*) NTPREV,PC_NNN,DTT(1,1,1),PI
C --- GET TIME AVERAGED VEL
      DO I=1,NC
        IF (I.GT.NH) I2=I-NH
        IF (I.LE.NH) I2=NH-I+1
        DO J=1,MR
          L=INDEXB(I,J)
          DO K=1,NTPREV ! time-wise, k can be seen as the time step, reverse direction of theta
            VVX=0.0E0
            VVR=0.0E0
            VVT=0.0E0
            PC_T0=ATAN2(XCT(L,3),XCT(L,2))/PI*180.0-REAL(K-1)*360.0/REAL(NTPREV)
            DO L=1,PC_NNN
              IF ((L.EQ.1).OR.(L.EQ.PC_NNN)) THEN
                PC_WT=1.0/(REAL(PC_NNN-1))/2.0
              ELSE
                PC_WT=1.0/(REAL(PC_NNN-1))
              END IF
              PC_T=DTT(I2,J,L)/PI*180.0
              PC_I=PC_T-PC_T0

              PC_I=PC_I*REAL(NTPREV)/360.0e0+1.0e0
              I0=FLOOR(PC_I)
              IF (I0.GT.PC_I) WRITE(*,*) '??',PC_I,I0
              WW0=PC_I-REAL(I0)
              IF (I0.GT.NTPREV) I0=I0-NTPREV
              IF (I0.GT.NTPREV) I0=I0-NTPREV
              IF (I0.GT.NTPREV) I0=I0-NTPREV
              IF (I0.LT.1) I0=I0+NTPREV
              IF (I0.LT.1) I0=I0+NTPREV
              IF (I0.LT.1) I0=I0+NTPREV
              WW1=1.0-WW0
              I1=I0+1
              IF (I1.GT.NTPREV) I1=1
              IF (I0.GE.(NTPREV+1)) WRITE(*,*) '?????',pc_i
*        IF ((I.EQ.18).AND.(J.EQ.6)) THEN
*          IF ((K.LT.4).OR.(K.GT.55)) THEN
*            WRITE(*,*) I0,WW0,PC_IVX(I,J,L,I0),PC_IVX(I,J,L,I1),PC_WT
*          END IF
*        END IF
             VVX=VVX+(PC_IVX(I2,J,L,I0)*WW1+PC_IVX(I2,J,L,I1)*WW0)*PC_WT
             VVR=VVR+(PC_IVR(I2,J,L,I0)*WW1+PC_IVR(I2,J,L,I1)*WW0)*PC_WT
             VVT=VVT+(PC_IVT(I2,J,L,I0)*WW1+PC_IVT(I2,J,L,I1)*WW0)*PC_WT
            END DO
            AVX(I,J,K)=VVX
            AVR(I,J,K)=VVR
            AVT(I,J,K)=VVT
          END DO
        END DO
      END DO

C --- correction #1=(#2+#ntprev)/2
      DO I=1,NC
        DO J=1,MR
          AVX(I,J,1)=(AVX(I,J,2)+AVX(I,J,NTPREV))/2.0
          AVR(I,J,1)=(AVR(I,J,2)+AVR(I,J,NTPREV))/2.0
          AVT(I,J,1)=(AVT(I,J,2)+AVT(I,J,NTPREV))/2.0
        END DO
      END DO
C --- correction hub base section extrapolation
      DO I=1,NC
        DO K=1,NTPREV
          T_TT =      (XCT(INDEXB(I,1),2)-XCT(INDEXB(I,2),2))
     &            /   (XCT(INDEXB(I,2),2)-XCT(INDEXB(I,3),2))
          AVX(I,1,K)=AVX(I,2,K)+(AVX(I,2,K)-AVX(I,3,K))*T_TT
          AVR(I,1,K)=AVR(I,2,K)+(AVR(I,2,K)-AVR(I,3,K))*T_TT
          AVT(I,1,K)=AVT(I,2,K)+(AVT(I,2,K)-AVT(I,3,K))*T_TT
        END DO
      END DO

C --- PLOT
      OPEN(729,FILE='prop_inv.plt',STATUS='UNKNOWN')
      WRITE(729,*) 'VARIABLES = "X","Y","Z","VX","VR","VT"'
      WRITE(729,*) 'ZONE T = "BLD",I=',NC,',J=',MR,',K=',NTPREV
      WRITE(729,*) NC,MR,NTPREV
      DO K=1,NTPREV
        DO J=1,MR
          DO I=1,NC
            L=INDEXB(I,J)
            TT_RR=SQRT(XCT(L,2)**2+XCT(L,3)**2)
            TT_TT=ATAN2(XCT(L,3),XCT(L,2))-REAL(K-1)*2.0*PI/REAL(NTPREV)
            XX=XCT(L,1)
            YY=TT_RR*COS(TT_TT)
            ZZ=TT_RR*SIN(TT_TT)*RVT
            WRITE(729,*) XX,YY,ZZ,AVX(I,J,K),AVR(I,J,K),RVT*AVT(I,J,K)
          END DO
        END DO
      END DO
      CLOSE(729)

      OPEN(728,FILE='Test_indv_disk.plt',STATUS='UNKNOWN')
      WRITE(728,*) 'VARIABLES="X","Y","Z","VX","VY","VZ","VR","VT"'
      WRITE(728,*) 'ZONE T = "BLD",I=',NH,',J=',MR,',K=',NTPREV+1
      DO K1=1,NTPREV+1
        K=K1
        IF (K.EQ.(NTPREV+1)) K=1
        DO J=1,MR
          DO I=1,NH
            L=INDEXB(NH+I,J)
            TT_RR=SQRT(XCT(L,2)**2+XCT(L,3)**2)
            TT_TT=ATAN2(XCT(L,3),XCT(L,2))-REAL(K-1)*2.0*PI/REAL(NTPREV)
            XX=XCT(L,1)
            YY=TT_RR*COS(TT_TT)
            ZZ=TT_RR*SIN(TT_TT)

            VVX=AVX(I,J,K)
            VVR=AVR(I,J,K)
            VVT=AVT(I,J,K)
            VVY=VVR*COS(TT_TT)-VVT*SIN(TT_TT)
            VVZ=VVR*SIN(TT_TT)+VVT*COS(TT_TT)

            WRITE(728,*) XX,YY,ZZ,VVX,VVY,VVZ,VVR,VVT
          END DO
        END DO
      END DO
      CLOSE(728)

      OPEN(727,FILE='Test_indv_blade.plt',STATUS='UNKNOWN')
      WRITE(727,*) 'VARIABLES="X","Y","Z","VX","VR","VT"'
      DO K=1,NTPREV
        WRITE(727,*) 'ZONE T ="TIME',K,'",I=',NC,',J=',MR,',K=',1
        DO J=1,MR
          DO I=1,NC
            L=INDEXB(I,J)
            TT_RR=SQRT(XCT(L,2)**2+XCT(L,3)**2)
            TT_TT=ATAN2(XCT(L,3),XCT(L,2))-REAL(K-1)*2.0*PI/REAL(NTPREV)
            XX=XCT(L,1)
            YY=TT_RR*COS(TT_TT)
            ZZ=TT_RR*SIN(TT_TT)
            VVX=AVX(I,J,K)
            VVR=AVR(I,J,K)
            VVT=AVT(I,J,K)
            WRITE(727,*) XX,YY,ZZ,VVX,VVR,VVT
          END DO
        END DO
      END DO
      CLOSE(727)

      OPEN(713,FILE='prop_EFF_XYZ.dat', STATUS='UNKNOWN')
      WRITE(713,*) NC*MR*NTPREV,1
      DO K=1,NTPREV
        DO M=1,MR
          DO N=1,NC
            L=INDEXB(N,M)
            TT_XX=XCT(L,1)
            TT_RR=SQRT(XCT(L,2)**2+XCT(L,3)**2)
            TT_TT=ATAN2(XCT(L,3),XCT(L,2))-REAL(K-1)*2.0*PI/REAL(NTPREV)

C Yiran Su 20161018 evaluate total velocity on the camber panel control point
            L=INDEXB(NC+1-N,M)
            TT_XX=(TT_XX+XCT(L,1))*0.5
            TT_RR=(TT_RR+SQRT(XCT(L,2)**2+XCT(L,3)**2))*0.5
            TT_TT=(TT_TT + ATAN2(XCT(L,3),XCT(L,2))
     &                   - REAL(K-1)*2.0*PI/REAL(NTPREV) )*0.5
C Finished

            WRITE(713,*) TT_XX,TT_RR*COS(TT_TT),TT_RR*SIN(TT_TT)*RVT
          END DO
        END DO
      END DO
      CLOSE(713)

      END SUBROUTINE


      SUBROUTINE STRENGTHREN(IT)
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INTEGER PCNTPOS
      NREAD = NWMIN*MR
      NREADD = NDWK*MDUCT
      IREC = PCNTPOS(1,IT)
      CALL READ2(45,IREC,POT,NPANEL)
      CALL READ2(47,IREC,DPDN,NPANEL)
      CALL READ2(46,IREC,TEMP5,NREAD)
      IF (IDUCT.NE.0) THEN
        CALL READ2(49,IREC,TEMPD4,NREADD)
      END IF
      END SUBROUTINE
      FUNCTION PCNTPOS(KB,KEYPOS)
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INTEGER PCNTPOS
      NTOT = 360 / NDLTAT
      NADV = 360 / NBLADE / NDLTAT
      PCNTPOS=KEYPOS+(KB-1)*NADV
      IF(PCNTPOS.GT.NTOT) THEN
         PCNTPOS=PCNTPOS-NTOT
      END IF
      RETURN
      END
      SUBROUTINE POLYGON(NP,X,Y,A)
C ===== calculating polygon area ======
C
C      p4  ----------- p3
C          |         |
C          |         |     ^ y (or 'r')
C          |         |     |
C          |         |     |
C   p5/p1  ----------- p2  --->x
C
C ====================================
      implicit none
      real*8, parameter :: eps = 1.0d-6
      integer np
      real x(np),y(np)
      real a
      real*8 chrlen, s
      real xa(np),ya(np)
      integer i

      chrlen = 0.0d0
      do i = 1, np - 1
        s =sqrt((x(i+1)-x(i))**2+(y(i+1)-y(i))**2)
        chrlen = max(chrlen, s)
      end do

      if (chrlen .lt. eps) then
        chrlen = eps
      end if

      xa = x/chrlen
      ya = y/chrlen

      a=0.0
      do i = 1, np - 1
        a=a+xa(i)*ya(i+1)-xa(i+1)*ya(i)
      end do
      a = a * 0.5

      a = a*(chrlen*chrlen)
      END SUBROUTINE
      LOGICAL FUNCTION IISD(P1,P2,P3,P4,PP)
C----------------------------------------
C               T3
C           P4------ P3
C           |  .     |
C        T4 |   PP   | T2
C           |        |
C           P1------ P2
C               T1
C----------------------------------------
      IMPLICIT NONE
      REAL,DIMENSION(2)::P1,P2,P3,P4,PP,PS,PE
      REAL::T1,T2,T3,T4,A,B,C

      IISD=(.FALSE.)
C - EQUATION
C (Y-YS)(XE-XS)=(X-XS)(YE-YS)
C  AX+BY+C=0  A=YE-YS  B=-1(XE-XS) C=YS(XE-XS)-XS(YE-YS)
      PS=P1
      PE=P2
      A=PE(2)-PS(2)
      B=PS(1)-PE(1)
      C=PS(2)*(PE(1)-PS(1))-PS(1)*(PE(2)-PS(2))
      T1=A*PP(1)+B*PP(2)+C

      PS=P2
      PE=P3
      A=PE(2)-PS(2)
      B=PS(1)-PE(1)
      C=PS(2)*(PE(1)-PS(1))-PS(1)*(PE(2)-PS(2))
      T2=A*PP(1)+B*PP(2)+C

      PS=P3
      PE=P4
      A=PE(2)-PS(2)
      B=PS(1)-PE(1)
      C=PS(2)*(PE(1)-PS(1))-PS(1)*(PE(2)-PS(2))
      T3=A*PP(1)+B*PP(2)+C

      PS=P4
      PE=P1
      A=PE(2)-PS(2)
      B=PS(1)-PE(1)
      C=PS(2)*(PE(1)-PS(1))-PS(1)*(PE(2)-PS(2))
      T4=A*PP(1)+B*PP(2)+C

      IF (((T1*T3).GE.0).AND.((T2*T4).GE.0)) IISD=(.TRUE.)

      RETURN
      END FUNCTION
      SUBROUTINE RPAN_PV(Xi,Yi,Zi,CHRLENSI,FSS,FDD,FSXI,FSYI,FDXI,FDYI,
     %                 FDZI,ID,IFLAG,XV1,YV1,SQ1,SIDE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL XI,YI,ZI,CHRLENSI,FSS,FDD,FSXI,FSYI,FDXI,FDYI,FDZI
      REAL XV1(4),YV1(4),SQ1(15),SIDE1(4)
      DIMENSION XV(4),YV(4),S(15),SIDE(4)
      DIMENSION R(4),RR(4),RI(4),XRI(4),YRI(4),FE(4),
     *          B(5),XMXV(4),YMYV(4),N1(4)

C----------------------------------------------------------------------
C
C  DATA ARRAY B FOR 6D RATIONAL APPROXIMATION OF ARCTANGENT
C  SEE HART, ET AL, COMPUTER APPROXIMATIONS, WILEY, 1968, TABLE 5090
C
C----------------------------------------------------------------------
      DATA B/ 2.4091197D-01, 3.7851122D+00, 5.6770721D+00,
     *        5.6772854D+00, 5.6770747D+00/
C      DATA PI/ 3.1415927D+00 /, PI2/ 1.5707963D+00 /
C      DATA TWOPI/ 6.2831853D+00 /
      DATA ONE/ 1.D+00/, A3/ 3.D+00/, A5/ 5.D+00/
      DATA A7/ 7.D+00/, A9/ 9.D+00/, A11/ 11.D+00/
      DATA A14/ 14.D+00/, A35/ 35.D+00/
      DATA A49/ 49.D+00/, A63/ 63.D+00/, A99/ 99.D+00/
C      DATA ONE10/.1D+00/, ONE6/.1666667D+00/, ONE3/ .3333333D+00/
C      DATA FIVE3/ 1.666667D+00/, SEVEN3/2.333333D+00/, ZERO/ 0.0D0/
      DATA ZERO/ 0.0D0/
      DATA FT3/ 4.666667D+00/, TOL/ 1D-08 /
      DATA N1/ 2, 3, 4, 1 /

      pi = dacos(-1.d0)
      pi2 = 0.5d0*pi
      twopi = 2.d0*pi
      one10 = 0.1d0
      one6 = 1.d0/6.d0
      one3 = 1.d0/3.d0
      five3 = 5.d0/3.d0
      seven3 = 7.d0/3.d0

C --- Change Input variable to double precision

      x = dble(xi)
      y = dble(yi)
      z = dble(zi)
      chrlens = dble(chrlensi)
      do i = 1 , 4
      xv(i) = dble(xv1(i))
      yv(i) = dble(yv1(i))
      side(i) = dble(side1(i))
      enddo
      do i = 1 , 15
        s(i) = dble(sq1(i))
      enddo

      XMXC=X-S(6)
      YMYC=Y-S(2)
      XX=XMXC*XMXC
      YY=YMYC*YMYC
      ZZ=Z*Z
      RRC=XX+YY+ZZ
      IF (RRC.LT.100.d0*CHRLENS) THEN
         IF(IFLAG.EQ.1) THEN
            GO TO 11
         ELSE IF (IFLAG.EQ.0) THEN
            IFLAG=2
            RETURN
         END IF
      END IF
C----------------------------------------------------------------------
C
C   TWO-TERM MULTIPOLE EXPANSION INCLUDING SECOND MOMENTS
C
C----------------------------------------------------------------------
      R2=ONE/RRC
      R1=DSQRT(R2)
      R3=R1*R2
      R5=R3*R2
      ZR2=Z*R2
      XY=XMXC*YMYC
      SS1=S(1)*R1
      SS3=-(S(3)+S(10))*R3
      SS5=(XX*S(10)+XY*S(7)+YY*S(3))*R5
      FS=SS1+ONE3*SS3+SS5
      FDSUM=SS1+SS3+A5*SS5
      FD=ZR2*FDSUM

C----------------------------------------------------------------------
C    DELETE NEXT  14  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 8
      RSS3=R2*SS1
      SSX3=-XMXC*RSS3
      SSY3=-YMYC*RSS3
      SSX5=(XMXC*(S(3)+A3*S(10))+YMYC*S(7))*R5
      SSY5=(YMYC*(S(10)+A3*S(3))+XMXC*S(7))*R5
      A5R2=A5*R2
      RSS7=-A5R2*SS5
      SSX7=XMXC*RSS7
      SSY7=YMYC*RSS7
      FSX=SSX3+SSX5+SSX7
      FSY=SSY3+SSY5+SSY7
      FDX=ZR2*(A3*SSX3+A5*SSX5+A7*SSX7)
      FDY=ZR2*(A3*SSY3+A5*SSY5+A7*SSY7)
      ZZR4=ZR2*ZR2
      FDZ=R2*FDSUM-ZZR4*(A3*SS1+A5*SS3+A35*SS5)
  8   IF (RRC.GT.15.d0*CHRLENS) GO TO 99
C----------------------------------------------------------------------
C
C    THIRD AND FOURTH MOMENTS ADDED FOR RRC/AREA BETWEEN 40 AND 150
C
C----------------------------------------------------------------------
      S914=S(9)+S(14)
      S813=S(8)+S(13)
      S411=S(4)+S(11)
      S512=S(5)+S(12)
      S1215=S(12)+S(15)
      R7=R5*R2
      R9=R7*R2
      SS5=(-XMXC*S813-YMYC*S411+ONE10*(S512+S1215))*R5
      SS7=(FIVE3*((XMXC*XX*S(13)+YMYC*YY*S(4))+A3*XY*(XMXC*S(11)
     *+YMYC*S(8)))-XX*S1215-YY*S512-XY*S914)*R7
      SS9=(A7*(ONE6*(XX*XX*S(15)+YY*YY*S(5))+XX*YY*S(12))
     *+SEVEN3*XY*(XX*S(14)+YY*S(9)))*R9
      FS=FS+SS5+SS7+SS9
      FDSUM=A5*SS5+A7*SS7+A9*SS9
      FD=FD+ZR2*FDSUM

C----------------------------------------------------------------------
C    DELETE NEXT  20  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 99
      TXY=XY+XY
      SSX5=-S813*R5
      SSY5=-S411*R5
      RSS7=A5R2*SS5
      SSX7=(A5*(XX*S(13)+TXY*S(11)+YY*S(8))-S1215*(XMXC+XMXC)
     *-YMYC*S914)*R7-XMXC*RSS7
      SSY7=(A5*(YY*S(4)+XX*S(11)+TXY*S(8))-S512*(YMYC+YMYC)
     *-XMXC*S914)*R7-YMYC*RSS7
      RSS9=A7*SS7*R2
      SSX9=(FT3*XMXC*XX*S(15)+A14*XMXC*YY*S(12)+A49*YMYC*(XX*S(14)
     *+ONE3*YY*S(9)))*R9-XMXC*RSS9
      SSY9=(FT3*YMYC*YY*S(5)+A14*YMYC*XX*S(12)+A49*XMXC*(YY*S(9)
     *+ONE3*XX*S(14)))*R9-YMYC*RSS9
      RSS11=A9*SS9*R2
      SSX11=-XMXC*RSS11
      SSY11=-YMYC*RSS11
      FSX=FSX+SSX5+SSX7+SSX9+SSX11
      FSY=FSY+SSY5+SSY7+SSY9+SSY11
      FDX=FDX+ZR2*(A5*SSX5+A7*SSX7+A9*SSX9+A11*SSX11)
      FDY=FDY+ZR2*(A5*SSY5+A7*SSY7+A9*SSY9+A11*SSY11)
      FDZ=FDZ+R2*FDSUM-ZZR4*(A35*SS5+A63*SS7+A99*SS9)
      GO TO 99
C----------------------------------------------------------------------
C
C    NEAR-FIELD SECTION USES EXACT FORMULATION
C      SET Z=TOL IF Z.LT.TOL TO AVOID INDETERMINACY ON PANEL
C      ZVTX IS USED TO DETERMINE PROXIMITY TO VERTEX NORMALS
C      MFLAG=1 IF NEAR VERTEX NORMALS
C
C----------------------------------------------------------------------
  11  FD=0.0d0
      FS=0.0d0
      ABZ=DABS(Z)
      IF(ABZ.GT.TOL) GO TO 12
        Z=TOL
        ZZ=Z*Z
        ABZ=TOL
  12  ZVTX=1.005D0*ABZ
      MFLAG=0
C----------------------------------------------------------------------
C    DELETE NEXT 5 LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      FSX=0.0d0
      FSY=0.0d0
      FDX=0.0d0
      FDY=0.0d0
      FDZ=0.0d0
C----------------------------------------------------------------------
C
C    LOOP FOR CORNER FUNCTIONS
C
C----------------------------------------------------------------------
      DO N=1,4
        XMXV(N)=X-XV(N)
        YMYV(N)=Y-YV(N)
        XX=XMXV(N)*XMXV(N)
        YY=YMYV(N)*YMYV(N)
        FE(N)=ZZ+XX
        RR(N)=FE(N)+YY
        R(N)=DSQRT(RR(N))
        IF (R(N).LT.ZVTX) MFLAG=1
        RI(N)=ONE/R(N)
C----------------------------------------------------------------------
C    DELETE NEXT   3   LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 13
        XRI(N)=RI(N)*XMXV(N)
        YRI(N)=RI(N)*YMYV(N)
   13   CONTINUE
      END DO
C----------------------------------------------------------------------
C
C    LOOP FOR SIDE FUNCTIONS AND SUMS OVER FOUR SIDES
C
C----------------------------------------------------------------------
      DO N=1,4
        IF (SIDE(N).LT.TOL) GO TO 33
        SIDI=ONE/SIDE(N)
        CT=(XV(N1(N))-XV(N))*SIDI
        ST=(YV(N1(N))-YV(N))*SIDI
        V=XMXV(N)*ST-YMYV(N)*CT
        VV=V*V
        RADS=VV+ZZ
        U1=XMXV(N)*CT+YMYV(N)*ST
        U2=XMXV(N1(N))*CT+YMYV(N1(N))*ST
        RSUM=R(N)+R(N1(N))
        FLAG=RI(N)*RI(N1(N))*U1*U2
C----------------------------------------------------------------------
C
C       FLAG=1 ON EXTENSIONS, -1 ON SIDES
C         IN FOLLOWING SUBSECTIONS FS,FSX,FSY,FDZ ARE EVALUATED FROM
C           LAST FORM OF (3.9) IN NORMAL CASE
C         ELSE
C           FIRST FORM OF (3.9) IF NEAR SIDE OF PANEL
C
C----------------------------------------------------------------------
        IF (FLAG.GT.-.99D0) THEN
          RSP=RSUM+SIDE(N)
          RSM=RSUM-SIDE(N)
          FLN=DLOG(RSP/RSM)
        ELSE
          RU1=R(N)+U1
          RU2=R(N1(N))-U2
          RADI=ONE/RADS
          FLN=DLOG(RU1*RU2*RADI)
        ENDIF

        FS=FS+V*FLN

C----------------------------------------------------------------------
C    DELETE FOLLOWING LINES TO NEXT ENDIF TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 14
        IF (FLAG.GT.-.99D0) THEN
          FAC=V*(SIDE(N)+SIDE(N))/(RSP*RSM)
          FSX=FSX+FLN*ST-FAC*(XRI(N)+XRI(N1(N)))
          FSY=FSY-FLN*CT-FAC*(YRI(N)+YRI(N1(N)))
          FDZ=FDZ-FAC*(RI(N)+RI(N1(N)))
        ELSE
          RU1I=ONE/RU1
          RU2I=ONE/RU2
          FA=RU1I-RU2I
          FB=-(V+V)*RADI
          FSX=FSX+FLN*ST+V*(FA*CT+FB*ST+RU1I*XRI(N)+RU2I*XRI(N1(N)))
          FSY=FSY-FLN*CT+V*(FA*ST-FB*CT+RU1I*YRI(N)+RU2I*YRI(N1(N)))
          FDZ=FDZ+FB+V*(RU1I*RI(N)+RU2I*RI(N1(N)))
        ENDIF
C----------------------------------------------------------------------
C
C        IN FOLLOWING SUBSECTIONS FACTORS IN (2.15) ARE EVALUATED FROM
C          (2.7) IN NORMAL CASE
C        ELSE
C          (2.14) IF NEAR NORMAL TO A VERTEX
C
C----------------------------------------------------------------------
   14   IF (MFLAG.EQ.0) THEN
          S1=V*R(N)
          C1=ABZ*U1
          S2=V*R(N1(N))
          C2=ABZ*U2
        ELSE
          FH1=XMXV(N)*YMYV(N)
          FH2=XMXV(N1(N))*YMYV(N1(N))
          S1=FE(N)*ST-FH1*CT
          C1=ABZ*R(N)*CT
          S2=FE(N1(N))*ST-FH2*CT
          C2=ABZ*R(N1(N))*CT
        ENDIF
        S12=S1*C2-S2*C1
        C12=C1*C2+S1*S2
C----------------------------------------------------------------------
C
C    EVALUATE THIRD ARCTANGENT IN (2.15)
C         ANGLE (MODULO PI) BETWEEN -PI/4 AND PI/4
C       ELSE
C         USE INVERSE COTANGENT AND ADD/SUBTRACT PI/2
C
C----------------------------------------------------------------------
        IF (DABS(S12).LE.DABS(C12)) THEN
          U=S12/C12
          IF (C12.LT.ZERO) FD=FD+DSIGN(PI,S12)
        ELSE
          U=-C12/S12
          FD=FD+DSIGN(PI2,S12)
        ENDIF
        UU=U*U
        FD=FD+U*((B(1)*UU+B(2))*UU+B(3))/((UU+B(4))*UU+B(5))
C----------------------------------------------------------------------
C
C    FOLLOWING THREE SECTIONS EVALUATE FDX,FDY FOR
C       FIELD POINT NEAR NORMAL TO VERTEX (MFLAG=1)
C       NORMAL CASE (FLAG.LT.0.99)
C       NEAR EXTENSIONS OF SIDES (FLAG.GE.0.99)
C
C    DELETE ALL FOLLOWING LINES ABOVE LABEL 33 TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 33
        IF (MFLAG.EQ.0) GO TO 20
          FAC=C1/((C1*C1+S1*S1)*RR(N))
        IF(Z.LT.ZERO) FAC=-FAC   ! Fix sign of derivatives, Newman
          FDX=FDX+(RR(N)*V+FH1*U1)*FAC
          FDY=FDY-FE(N)*U1*FAC
          FAC=C2/((C2*C2+S2*S2)*RR(N1(N)))
        IF(Z.LT.ZERO) FAC=-FAC   ! Fix sign of derivatives, Newman
          FDX=FDX-(RR(N1(N))*V+FH2*U2)*FAC
          FDY=FDY+FE(N1(N))*U2*FAC
          GO TO 33
  20    IF (FLAG.LT.0.99D0) THEN
          U1V=U1*V
          FAC=Z/(C1*C1+S1*S1)
          FDX=FDX+(U1V*XRI(N)+R(N)*YMYV(N))*FAC
          FDY=FDY+(U1V*YRI(N)-R(N)*XMXV(N))*FAC
          U2V=U2*V
          FAC=Z/(C2*C2+S2*S2)
          FDX=FDX-(U2V*XRI(N1(N))+R(N1(N))*YMYV(N1(N)))*FAC
          FDY=FDY-(U2V*YRI(N1(N))-R(N1(N))*XMXV(N1(N)))*FAC
        ELSE
          ZS=Z*SIDE(N)
          USUM=U1+U2
          VRADS=V*RADS
          SFAC=VRADS*USUM
          SFS=-SFAC*ZS
          SFA=SFAC*C12
          CFAC=U2*R(N)+U1*R(N1(N))
          SFB=SFAC*CFAC
          CCF=C12*CFAC
          PA=(CCF+CCF)*VRADS-SFA*RSUM-SFB*ZZ*USUM
          PB=CCF*USUM*(VV+VV+RADS)-SFB*(S1+S1)*R(N1(N))
          PC=-SFA*U2-SFB*VV*R(N1(N))
          PD=-SFA*U1-SFB*VV*R(N)
          FAC=ZS/(CCF*CCF+SFS*SFS)
          FDX=FDX-(PA*CT+PB*ST+PC*XRI(N)+PD*XRI(N1(N)))*FAC
          FDY=FDY-(PA*ST-PB*CT+PC*YRI(N)+PD*YRI(N1(N)))*FAC
        ENDIF
  33    CONTINUE
      END DO
      IF (FD.LT.ZERO) FD=FD+TWOPI
      IF (Z.LT.ZERO) FD=-FD
      FS=FS-Z*FD

C----------------------------------------------------------------------
C    DELETE NEXT  3  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 99
      FSX=FSX-Z*FDX
      FSY=FSY-Z*FDY

99    CONTINUE

      fdd = sngl(fd)
      fss = sngl(fs)

      if(id .ne. 0) then
        fsxi = sngl(fsx)
        fsyi = sngl(fsy)
        fdxi = sngl(fdx)
        fdyi = sngl(fdy)
        fdzi = sngl(fdz)
      endif

      RETURN
      END
