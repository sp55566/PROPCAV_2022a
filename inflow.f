      SUBROUTINE INFLOW
************************************************************************
*     INFLOW: INFLOW velocities                                        *
*      --- Compute the source strength of each blade at current time   *
*            step                                                      *
*      --- Compute the total inflow velocities at the key blade        *
*                                                                      *
* Date             Revision                                            *
* ----             --------                                            *
* Feb, 2008    Hong changed the inflow conditions on duct and tunnel.  *
*        Note: The duct/tunnel is solved w.r.t. the propeller system.  *
*              When calculating the velocities on the duct/tunnel,     *
*              the rotational velocity should be subtracted from the   *
*              total velocity (w.r.t propeller system) calculated      *
*              by prsdif.f.                                            *
************************************************************************
      use m_BLADEM
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      DIMENSION U(3),NH1(3)
      REAL, ALLOCATABLE :: XVTMP(:), ATMP(:,:,:), BTMP(:,:,:)
      REAL, ALLOCATABLE :: UXTMP(:,:),UYTMP(:,:),UZTMP(:,:),TTMP(:)
      
      REAL RVT
      LOGICAL IRVT
      RVT=1.0
      INQUIRE(FILE='rvt.inp',EXIST=IRVT)
      IF (IRVT.EQ.(.TRUE.)) THEN
        RVT=-1.0
        IF (NTSTEP.EQ.0) WRITE(*,*) 'Left handed, reverse Ueff_swirl!'
      END IF

C-----------------------------------------------------------------------
C     Calculate the source strengths  ( DPDN = -V_inf*N, V_inf=u+wr )
C-----------------------------------------------------------------------
      IF (NTSTEP.EQ.0.OR.ISTEADY.EQ.0) THEN
        DO 10 IH=1,3
          NH1(IH)=1
   10   CONTINUE
      ELSE
        DO 15 IH=1,3
          NH1(IH)=NHARM(IH)
   15   CONTINUE
      END IF

      NNPAND1 = 0
      NNPAND2 = 0
      IF(IDUCT .EQ. 1) THEN
         NNPAND1 = NPANB + NPANH
         NNPAND2 = NPANB + NPANH + NPAND
      ENDIF

      NNPAN3 = 0
      NNPAN4 = 0
      IF(ITUN .NE. 0) THEN
         NNPAN3 = NPANB + NPANH + NPAND
         NNPAN4 = NPANB + NPANH + NPAND + NPANTN
      ENDIF

C #1

C ---------------
C Comments from Yiran Su:
C
C In PROPCAV, blade angle index increase in the opposite direction of \theta
C Here KK0 is the blade angle index
C KK is the sequence of blades stored in XCTP array only !
C ---------------

      DO 100 KK=1,NBLADE
         DTBLA=DELK*FLOAT(KK-1)
         DO 40 J=1,NPANEL

            RCP=SQRT(XCTP(J,2,KK)**2+XCTP(J,3,KK)**2)
            THP=ATAN2(XCTP(J,3,KK),XCTP(J,2,KK))

            DO IH=1,3
               U(IH)=0.0
               DO JH=1,NH1(IH)

                  CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
                  CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
                  T=FLOAT(JH-1)*(THP-TSTEP)
                  ST=SIN(T)
                  CT1=COS(T)
                  U(IH)=U(IH)+COEC*CT1+COES*ST
               ENDDO
            ENDDO

            VOX=U(1)
            VOR=U(2)
            VOT=U(3)

CYiranSu_PC2NS_read_effective wake
            IF (IPC2NS.EQ.(.TRUE.)) THEN
              IF (J.LE.(NC*MR)) THEN
                M=MR-FLOOR(REAL(J-1)/REAL(NC))
                KK0 = 1
                IF(KK .NE. 1) then
                  KK0 = NBLADE + 2 - KK
                END IF
                N=J-(MR-M)*NC
                K0=IDXREV+NINT(REAL(KK0-1)/REAL(NBLADE)*REAL(NTPREV))
                IF (K0.GT.NTPREV) K0=K0-NTPREV

                IF (ISTEADY.EQ.0) THEN  ! steady case
                  VOX=UEFX(N,M,1)
                  VOR=UEFR(N,M,1)
                  VOT=UEFT(N,M,1)*RVT

                ELSEIF (IDXREV.EQ.0) THEN  ! unsteady case steady step
                  VOX=0.0
                  VOR=0.0
                  VOT=0.0
                  DO II=1,NTPREV
                    VOX=VOX+UEFX(N,M,II)
                    VOR=VOR+UEFR(N,M,II)
                    VOT=VOT+UEFT(N,M,II)*RVT
                  END DO
                  VOX=VOX/REAL(NTPREV)
                  VOR=VOR/REAL(NTPREV)
                  VOT=VOT/REAL(NTPREV)

                ELSE                    ! unsteady steps

                  VOX=UEFX(N,M,K0)
                  VOR=UEFR(N,M,K0)
                  VOT=UEFT(N,M,K0)*RVT
                END IF
              ELSE
                VOX=1.0
                VOR=0.0
                VOT=0.0
              END IF
            END IF
CYiranSu_End

C#2

            COCP=XCTP(J,2,KK)/RCP
            SICP=XCTP(J,3,KK)/RCP
            WROVS=PI*RCP/ADVCO

C/s S. KIM| local friction coefficients
            IF (KK.EQ.1) THEN
               IF (J.LE.(NC*MR)) THEN
                  M = MR - FLOOR(REAL(J-1)/REAL(NC))
                  N = J - (MR-M)*NC 
                  VRR(N,M) = SQRT(VOX*VOX + WROVS*WROVS)
                  RER(N,M) = VRR(N,M)*CHORD(M)*XCDF ! In case XCDF > 10.0
                  XCDF1(N,M) = 0.075/(LOG10(RER(N,M))-2.0)**2.0 
                  IF (XCDF.LE.10.0) XCDF1(N,M) = XCDF ! In case XCDF <= 10.0, use a constant CF.
               ENDIF  
            ENDIF
C/e S. KIM| local friction coefficients

            IF(IDUCT.EQ.1.AND.IDOPT.EQ.1) THEN  !Solve duct only
               IF(J .GT. NNPAND1 .AND. J .LE. NNPAND2) THEN
                  WROVS = 0.0
               ENDIF
            END IF

            VOY=VOR*COCP-(WROVS+VOT)*SICP
            VOZ=VOR*SICP+(WROVS+VOT)*COCP

C-----------------------------------------------------------------------
C         If running hydrofoil case, set voy and voz to 0       CM121797
C-----------------------------------------------------------------------
            IF(ICON .EQ. 5 .OR. ICON .EQ. 6 .OR. ICON .EQ. 8) THEN
               VOX = COS(ALPHA * RAD)
               VOY = 0.0
               VOZ = SIN(ALPHA * RAD)
            ENDIF

            VELY=VEL(J,2)*COS(DTBLA)-VEL(J,3)*SIN(DTBLA)
            VELZ=VEL(J,2)*SIN(DTBLA)+VEL(J,3)*COS(DTBLA)
            BUG=VOX*VEL(J,1)+VOY*VELY+VOZ*VELZ

            IF(KK .EQ. 1) then
               KK0 = 1
               DPDN(J)=-BUG
               VOX1(J)=VOX
               VOY1(J)=VOY
               VOZ1(J)=VOZ
            ELSE
               KK0 = NBLADE + 2 - KK
            ENDIF
            STRGTH(J,KK0)=-BUG
 40      CONTINUE

C --- S.H.CHANG 04/20/2010 FOR CONTINUITY INSIDE OF A TUNNEL
         IF(ITUN.NE.0) THEN
           DO M = 1,MTUNEL
              DO N = 1,NAXT
                 J = INDEXTN(N,M)

                 RCP=SQRT(XCTP(J,2,KK)**2+XCTP(J,3,KK)**2)
                 COCP=XCTP(J,2,KK)/RCP
                 SICP=XCTP(J,3,KK)/RCP

                 WROVS = PI*RCP/ADVCO
                 VOX = 1.0
                 VOY = -WROVS*SICP
                 VOZ = WROVS*COCP

                 IF(N .GT. (NAXT-NSIDE)) THEN
                   BUG = 1.0-AREAR
                 ELSE
                   BUG = VOX*VEL(J,1)
                 END IF

                 IF(KK .EQ. 1) THEN
                   KK0 = 1
                   DPDN(J)=-BUG
                   VOX1(J)=VOX
                   VOY1(J)=VOY
                   VOZ1(J)=VOZ
                 ELSE
                   KK0 = NBLADE + 2 - KK
                 ENDIF

                 STRGTH(J,KK0) = -BUG

              END DO
           END DO
         END IF
C --- S.H.CHANG 04/20/2010 FOR CONTINUITY INSIDE OF A TUNNEL

 100  CONTINUE

C-----------------------------------------------------------------------
C           On-coming velocity in local coordinate system
C-----------------------------------------------------------------------

      DO N=1,NC
         DO M=1,MR
            L=INDEXB(N,M)
            VXIB(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
     *           +VOZ1(L)*DIR(L,1,3)
            VETAB(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
     *           +VOZ1(L)*DIR(L,2,3)
            VINFSB(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
         ENDDO
      ENDDO

      IF(IHUB .NE. 0) THEN
         DO N=1, NHBX
            DO M=1, MHBT
               L=INDEXH(N,M)
               VXIH(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
     *              +VOZ1(L)*DIR(L,1,3)
               VETAH(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
     *              +VOZ1(L)*DIR(L,2,3)
               VINFSH(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
            ENDDO
         ENDDO
      END IF

      IF(IDUCT .NE. 0) THEN
         DO M=1, MDUCT
            DO N=1, NDUCT
               L=INDEXD(N,M)
               VXID(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
     *              +VOZ1(L)*DIR(L,1,3)
               VETAD(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
     *              +VOZ1(L)*DIR(L,2,3)
               VINFSD(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
            ENDDO
         ENDDO
      END IF

      IF(ITUN .NE. 0) THEN
         DO N=1, NAXT
            DO M=1, MTUNEL
               L=INDEXTN(N,M)
               VXITN(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
     *              +VOZ1(L)*DIR(L,1,3)
               VETATN(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
     *              +VOZ1(L)*DIR(L,2,3)
               VINFSTN(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
            ENDDO
         ENDDO
      END IF

C**********************************************************************
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
C**********************************************************************
cC -- Begin Tip HSLEE(10/12/99)
c
c      IF(IAN .EQ. 2) THEN
c         DO N=1,NTHX
c            DO M=1, MCVT
c                 L=INDEXT(N,M)
c                 VXITH(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
c     *                     +VOZ1(L)*DIR(L,1,3)
c                 VETATH(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
c     *                      +VOZ1(L)*DIR(L,2,3)
c                 VINFSTH(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
c               ENDDO
c             ENDDO
c
c             DO N=1,NCVX
c               DO M=1, MCVT
c                 L=INDEXC(N,M)
c                 VXIC(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
c     *                     +VOZ1(L)*DIR(L,1,3)
c                 VETAC(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
c     *                      +VOZ1(L)*DIR(L,2,3)
c                 VINFSC(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
c               ENDDO
c             ENDDO
c      ENDIF
c
cC -- End Tip (10/12/99)
C**********************************************************************
C/e S.N.KIM | Aug. 2018.
C**********************************************************************
! T.WU plot the inflow wake
      IF (NTSTEP.EQ.0.AND.ICAVT.EQ.1) THEN
!
         ! ALLOCATE(XVTMP(MRP), ATMP(MRP,20,3), BTMP(MRP,20,3))
         ! ALLOCATE(UXTMP(MRP,NTPREV+1),UYTMP(MRP,NTPREV+1),UZTMP(MRP,NTPREV+1))
         ! ALLOCATE(TTMP(NTPREV+1))
         ALLOCATE(XVTMP(301), ATMP(201,20,3), BTMP(201,20,3))
         ALLOCATE(UXTMP(201,301),UYTMP(201,301),UZTMP(201,301))
         ALLOCATE(TTMP(301))
!
         OPEN(091725,FILE='UE_WAK.plt',STATUS='UNKNOWN')
         WRITE(091725,*) 'VARIABLES=Y,Z,UX,UY,UZ'
         WRITE(091725,*) 'ZONE T = "BLD",I=',MRP,',J=',NTPREV+1,',K=',1

         TTMP(1) = 0.0
         DO K = 1 , NTPREV
            TTMP(K+1) = TTMP(K) + DELTAT
         ENDDO
!
         DO I = 1 , 3
            DO K = 1 , 2
               DO J = 1 , NHARM(I)
                  CALL EVALDK(NWKCOE,MRP,XRW,RZ,XVTMP,XWCUB(1,J,K,I))
                  DO M = 1 , MRP
                     IF(K .EQ. 1) THEN
                        ATMP(M,J,I) = XVTMP(M)
                     ELSEIF(K .EQ. 2) THEN
                        BTMP(M,J,I) = XVTMP(M)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!
         DO M = 1 , MRP
            DO K = 1 , NTPREV+1
               DO J = 1 , NHARM(1)
                  UXTMP(M,K) = UXTMP(M,K) + 
     &                   ATMP(M,J,1) * COS( (J-1) * TTMP(K) ) +
     &                   BTMP(M,J,1) * SIN( (J-1) * TTMP(K) )
               ENDDO
               DO J = 1 , NHARM(2)
                  UYTMP(M,K) = UYTMP(M,K) + 
     &                   ATMP(M,J,2) * COS( (J-1) * TTMP(K) ) +
     &                   BTMP(M,J,2) * SIN( (J-1) * TTMP(K) )
               ENDDO
               DO J = 1 , NHARM(3)
                  UZTMP(M,K) = UZTMP(M,K) + 
     &                   ATMP(M,J,3) * COS( (J-1) * TTMP(K) ) +
     &                   BTMP(M,J,3) * SIN( (J-1) * TTMP(K) )
               ENDDO
            ENDDO
         ENDDO
!
         DO K = 1 , NTPREV+1
            DO M = 1 , MRP
               UR = UYTMP(M,K)
               UT = UZTMP(M,K)
               
               UYTMP(M,K) = UR * COS(TTMP(K)) 
     &              - UT * SIN(TTMP(K))
               UZTMP(M,K) = UR * SIN(TTMP(K)) 
     &              + UT * COS(TTMP(K))
            ENDDO
         ENDDO
         DO K = 1 , NTPREV+1
            DO M = 1 , MRP
               Y = RZ(M) * COS(TTMP(K))
               Z = RZ(M) * SIN(TTMP(K))
               WRITE(091725,*) Y, Z, UXTMP(M,K), UYTMP(M,K), UZTMP(M,K)
            ENDDO
         ENDDO
         DEALLOCATE(XVTMP, ATMP, BTMP)
         DEALLOCATE(UXTMP, UYTMP, UZTMP, TTMP)
      ENDIF
      RETURN
C>>>>>>>>>>>>>>>>>>>>>>End of subroutine INFLOW>>>>>>>>>>>>>>>>>>>>>>>>>
      END

C#1
*      REAL,DIMENSION(20,18,60)::UEX6,UEY6,UEZ6,UER6,UET6
*      REAL,DIMENSION(NH,MR,60)::UEX2,UER2,UET2
*      LOGICAL ISUTEM
*CYiranSu read wake from PF2NS effective wake file (NS2MPUF.dat)
*      INQUIRE(FILE='NS2MPUF.dat',EXIST=ISUTEM)
*      IF (ISUTEM.EQ.(.TRUE.)) THEN
*        WRITE(*,*) 'Read pf2ns wake: ISTEP=',IDXREV,',NREV=',NREV,NTSTEP
*        OPEN(991,FILE='NS2MPUF.dat',STATUS='OLD')
*        OPEN(992,FILE='EFF_XYZ.dat',STATUS='OLD')
*        READ(992,*) I,J
*
*        DO K=1,60
*          WRITE(993,*) 'ZONE T="AGL",I=',21,',J=',19,',K=',1,','
*
*          DO J=1,18
*            DO I=1,20
*              READ(991,*) TXX,UEX6(I,J,K),UEY6(I,J,K),UEZ6(I,J,K),TPP
*              READ(992,*) T_XX,T_YY,T_ZZ
*              T_TT=ATAN2(T_ZZ,T_YY)
*              UER6(I,J,K)=UEY6(I,J,K)*COS(T_TT)+UEZ6(I,J,K)*SIN(T_TT)
*              UET6(I,J,K)=UEZ6(I,J,K)*COS(T_TT)-UEY6(I,J,K)*SIN(T_TT)
*            END DO
*          END DO
*        END DO
*
*        CLOSE(991)
*        CLOSE(992)
*
*
*        !Interpolation begin
*        DO J=1,MR
*          DO I=1,NH
*            TII=(REAL(I)-0.5)/REAL(NH)*20.0+0.5
*            II=FLOOR(TII)
*            IF (II.LT.1) II=1
*            IF (II.GT.19) II=19
*            TII=TII-REAL(II)
*
*            TJJ=(REAL(J)-0.5)/REAL(MR)*18.0+0.5
*            JJ=FLOOR(TJJ)
*            IF (JJ.LT.1) JJ=1
*            IF (JJ.GT.17) JJ=17
*            TJJ=TJJ-REAL(JJ)
*
*            UEFF(I,J,1)=0.0
*            UEFF(I,J,2)=0.0
*            UEFF(I,J,3)=0.0
*
*            DO K=1,60
*              UEX2(I,J,K)=UEX6(II,JJ,K)*(1.0-TII)*(1.0-TJJ)+
*     &                    UEX6(II+1,JJ,K)*(TII)*(1.0-TJJ)+
*     &                    UEX6(II,JJ+1,K)*(1.0-TII)*(TJJ)+
*     &                    UEX6(II+1,JJ+1,K)*(TII)*(TJJ)
*              UER2(I,J,K)=UER6(II,JJ,K)*(1.0-TII)*(1.0-TJJ)+
*     &                    UER6(II+1,JJ,K)*(TII)*(1.0-TJJ)+
*     &                    UER6(II,JJ+1,K)*(1.0-TII)*(TJJ)+
*     &                    UER6(II+1,JJ+1,K)*(TII)*(TJJ)
*              UET2(I,J,K)=UET6(II,JJ,K)*(1.0-TII)*(1.0-TJJ)+
*     &                    UET6(II+1,JJ,K)*(TII)*(1.0-TJJ)+
*     &                    UET6(II,JJ+1,K)*(1.0-TII)*(TJJ)+
*     &                    UET6(II+1,JJ+1,K)*(TII)*(TJJ)
*          ! avoid extreme values
*              IF (UEX2(I,J,K).LT.0.6) UEX2(I,J,K)=0.6
*              IF (UEX2(I,J,K).GT.1.2) UEX2(I,J,K)=1.2
*              IF (UET2(I,J,K).LT.(-0.15)) UET2(I,J,K)=-0.15
*              IF (UET2(I,J,K).GT.(-0.15)) UET2(I,J,K)=0.15
*              UER2(I,J,K)=0.0
*
*              UEFF(I,J,1)=UEFF(I,J,1)+UEX2(I,J,K)
*              UEFF(I,J,2)=UEFF(I,J,2)+UER2(I,J,K)
*              UEFF(I,J,3)=UEFF(I,J,3)+UET2(I,J,K)
*            END DO
*            UEFF(I,J,1)=UEFF(I,J,1)/60.0
*            UEFF(I,J,2)=UEFF(I,J,2)/60.0
*            UEFF(I,J,3)=UEFF(I,J,3)/60.0
*          END DO
*        END DO
*        ! Finish interpolation
*
*        ! TEST PLOTS
**          OPEN(993,FILE='Test2_prop_panel_60.plt', STATUS='UNKNOWN')
**          WRITE(993,*) 'VARIABLES="X","Y","Z","VX","VR","VT"'
**          DO K=1,60
**            WRITE(993,*) 'ZONE T="AGL",I=',NH+1,',J=',MR+1,',K=',1,','
**            WRITE(993,*) 'DATAPACKING=BLOCK, VARLOCATION=(4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED)'
**            WRITE(993,*) ((XBM(I,J),I=1,(NH+1)),J=1,(MR+1))
**            WRITE(993,*) ((YBM(I,J),I=1,(NH+1)),J=1,(MR+1))
**            WRITE(993,*) ((ZBM(I,J),I=1,(NH+1)),J=1,(MR+1))
**            WRITE(993,*) ((UEX2(I,J,K),I=1,NH),J=1,MR)
**            WRITE(993,*) ((UER2(I,J,K),I=1,NH),J=1,MR)
**            WRITE(993,*) ((UET2(I,J,K),I=1,NH),J=1,MR)
**          END DO
**          CLOSE(993)
**          STOP
*
*      END IF
*CYiranSu_end

C#2
CYiranSu read wake from PF2NS effective wake file (NS2MPUF.dat)
*            IF (ISUTEM.EQ.(.TRUE.)) THEN
*              IF (J.LE.(NC*MR)) THEN
*                M=MR-FLOOR(REAL(J-1)/REAL(NC))
*                N=J-(MR-M)*NC
*                IF (N.GT.NH) THEN
*                  N2=N-NH
*                ELSE
*                  N2=NH-N+1
*                END IF
*
*                K=IDXREV-(KK-1)*15
*                IF (K.LE.0) K=K+60
*
*                K1=K+1
*                IF (K1.GT.60) K1=K1-60
*                K2=K-1
*                IF (K2.LT.1) K2=K2+60
*
*                IF (NTSTEP.EQ.0) THEN
*                  VOX=UEFF(N2,M,1)
*                  VOR=UEFF(N2,M,2)
*                  VOT=UEFF(N2,M,3)
*                ELSE
*                  VOX=(UEX2(N2,M,K)*2.0+UEX2(N2,M,K1)+UEX2(N2,M,K2))/4.0
*                  VOR=(UER2(N2,M,K)*2.0+UER2(N2,M,K1)+UER2(N2,M,K2))/4.0
*                  VOT=(UET2(N2,M,K)*2.0+UET2(N2,M,K1)+UET2(N2,M,K2))/4.0
*                END IF
*              END IF
*            END IF
CYiranSu_end




