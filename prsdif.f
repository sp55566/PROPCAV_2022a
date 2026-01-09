      SUBROUTINE PRSDIF(IT)

************************************************************************
*     PRSDIF: PReSsure calculations by DIFferentiating the potentials  *
*      --- Calculate the pressure coefficients on the blade and the    *
*          hub                                                         *
*                                                                      *
*      --- IT is the numer of iteration inside pressure kutta condition*
*             however, if IT = 999, PRSDIF will skip the calculation   *
*             of DPhi.                                                 *
*                                                                      *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     PARAMETER (MBZCUB=4*MBZ-4)

!     DIMENSION DPCUB(MBZCUB),RSQ1(MBZ)
      DIMENSION DPCUB(4*MBZ-4),RSQ1(MBZ)
      DIMENSION ICAV(NBZ),NCAV(NBZ),NDIR(NBZ),EPRES(NC,MR)
C      DIMENSION DELPDUM(MBZ)
      integer MBZCUB
      MBZCUB = 4*MBZ-4

C....DTN is the nondimensionalized delta(t).............................
      DTN=DELTAT*ADVCO/PI

      ISR=1
      IF(IFACE.EQ.2) ISR=2

C-----------------------------------------------------------------------
C     U velocity component along curvilinear coordinate (u,v,w) on the
C       blade by piecewise quadratic potentials
C     Notice that DPDXQ here gives -D(phi)/Dx because the
C       direction of u in RPAN is opposite to the direction of
C       'x' along the section. ('x' goes from the T.E. and goes
C       thru the lower surface, L.E., upper surface, to the T.E.
C       again
C-----------------------------------------------------------------------
      NCM=NC-1
      NCM2=NC-2
      DO 10 M=1,MR
         DO 20 N=2,NCM
            L0= INDEXB(N,M)
            L2= INDEXB(N-1,M)
            L1= INDEXB(N+1,M)
            X1=HALF*(DELU(L1)+DELU(L0))
            X2=HALF*(DELU(L2)+DELU(L0))
            DPDUB(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
 20      CONTINUE
         DPDUB(NC,M)=DPDUB(NCM,M)+(DPDUB(NCM,M)-DPDUB(NCM2,M))*X1/X2

         L0=INDEXB(2,M)
         L2=INDEXB(1,M)
         L1=INDEXB(3,M)
         DPDUB(1,M)=DPDUB(2,M)+(DPDUB(2,M)-DPDUB(3,M))
     *        *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))
 10   CONTINUE

C-----------------------------------------------------------------------
C     V velocity component along curvilinear coordinate (u,v,w) on the
C       blade by piecewise quadratic potentials
C-----------------------------------------------------------------------
      MRM=MR-1

      DO 30 N=1,NC
         DO 40 M=2,MRM
            L0= INDEXB(N,M)
            L1= INDEXB(N,M-1)
            L2= INDEXB(N,M+1)
            X1=0.5*(DELV(L1)+DELV(L0))
            X2=0.5*(DELV(L2)+DELV(L0))
            DPDVB(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
 40      CONTINUE
         DPDVB(N,MR)=DPDVB(N,MRM)
     *        +(DPDVB(N,MRM)-DPDVB(N,MR-2))*X2/X1

         L0= INDEXB(N,2)
         L1= INDEXB(N,1)
         L2= INDEXB(N,3)
         DPDVB(N,1)=DPDVB(N,2)+(DPDVB(N,2)-DPDVB(N,3))
     *        *(DELV(L1)+DELV(L0))/(DELV(L2)+DELV(L0))
 30   CONTINUE

C-----------------------------------------------------------------------
C     Calculate the pressure coefficients on the blade
C-----------------------------------------------------------------------
      THETAT=-FLOAT(NTSTEP-1)*DELTAT
      ADVCO2=ADVCO*ADVCO
      PI2=PI*PI

      ! Yiran Su calculate the effective pressure 02/18/2018
      EPRES = 0.0
C      if (IUNS.eq.1 .and. NTSTEP_KEY.GT.140) CALL PRES_EFF_CALC(EPRES)

C      write(*,*) 'THETAT =',THETAT*180./PI,NTSTEP

      DO 50 M=1,MR

C.......Flag with cavitating panels (including the split panel w/ ICAV..
         DO 60 N=1,NC
            ICAV(N)=0
 60      CONTINUE

         IF(IWET.EQ.0) THEN
            DO 70 II=1,ISR
               IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN
                  IDR=1
                  K=1
                  ISF=0
               ELSE
                  IDR=2
                  K=-1
                  ISF=1
               END IF
               M0M1=M0(M,IDR)-1
               DO 80 N=1,JCV(M,IDR)+NSPP(M,IDR)
                  N1=M0M1+K*N+ISF
                  ICAV(N1)=1
                  NCAV(N1)=N
                  NDIR(N1)=IDR
 80            CONTINUE
 70         CONTINUE
         END IF

         DO 90 N=1,NC
            L=INDEXB(N,M)

C..........Perturbation velocities in local coordinate system...........
            UXI=DPDUB(N,M)
            UETA=(DPDVB(N,M)-DPDUB(N,M)*SINPHI(L)) / COSPHI(L)

C..........Total velocity components in global coordinate system........
            UXIT=VXIB(N,M)+UXI
            UETAT=VETAB(N,M)+UETA

C..........If ITAN=1, the crossflow will be zeroed. (JY053000)
            ITAN=0
            IF(IWET.EQ.0.AND.ICON.NE.5.AND.ICON.NE.6) THEN
c               IF(ISP.EQ.1.OR.IFACE.EQ.2) THEN
               IF(ISP.EQ.1) THEN
                  ITAN=1
               END IF
               IF(ICON.EQ.7) ITAN=1
            END IF

            IF(ITAN.EQ.1) THEN
               UETAT=ZERO
               UT=VOX1(L)*VL(L,1)+VOY1(L)*VL(L,2)+VOZ1(L)*VL(L,3)
               IF(N.LE.NC/2) THEN
                  DPDVB(N,M)=ABS(UXIT)*SINPHI(L)-UT
               ELSE
                  DPDVB(N,M)=-ABS(UXIT)*SINPHI(L)-UT
               END IF
            END IF

C..........Magnitude and direction of total velocity in local system....
            if (abs(UXIT).gt.1e10.or.abs(UETAT).gt.1e10) then
              write(*,*) L,UXIT, UETAT
            end if
            VTOTS(L)=UXIT**2.+UETAT**2.
            UTOTU=UXIT/SQRT(VTOTS(L))
            UTOTW=UETAT/SQRT(VTOTS(L))

            IF(ICAV(N).EQ.1) THEN
               UTS=DPHIDS(NCAV(N),M,NDIR(N))
               UTV=DPHIDV(NCAV(N),M,NDIR(N))
               IF(NDIR(N).EQ.2) THEN
                  UTW=(UTV-UTS*SINPHI(L))/COSPHI(L)
               ELSE IF(NDIR(N).EQ.1) THEN
                  UTW=(UTV+UTS*SINPHI(L))/COSPHI(L)
               END IF
               VTOTS(L)=UTS**2.+UTW**2.
               UTOTU=UTS/SQRT(VTOTS(L))
               IF(NDIR(N).EQ.1) UTOTU=-UTOTU
               UTOTW=UTW/SQRT(VTOTS(L))
               UXIT=SQRT(VTOTS(L))*UTOTU
               UETAT=SQRT(VTOTS(L))*UTOTW
            END IF

            UXTOT(N,M)=UXIT*DIR(L,1,1)+UETAT*DIR(L,2,1)
            UYTOT(N,M)=UXIT*DIR(L,1,2)+UETAT*DIR(L,2,2)
            UZTOT(N,M)=UXIT*DIR(L,1,3)+UETAT*DIR(L,2,3)

C.........graviational term.............................................
            RCP=SQRT(XCT(L,2)*XCT(L,2)+XCT(L,3)*XCT(L,3))

            THETA=THETAT+ACOS(XCT(L,2)/RCP)
            RCP2=RCP*RCP
            IF((NTSTEP.EQ.0).OR.(ISTEADY.EQ.0).OR.(ISTEADY.EQ.3))THEN
               TWOGY = 0.
               IF((ITERGB.EQ.1).OR.(IFLAP.EQ.1))THEN
                  ZTEMP = RCP*COS(THETA)
                  TWOGY=-ZTEMP*2./(FNRUDDER)
               ENDIF
            ELSE
               TWOGY=RCP*COS(THETA)/(FROUDE*ADVCO2)
            ENDIF

C..........time dericavtive of potential................................
C
            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NTSTEP.EQ.0) THEN
               DPDT=0.0
            ELSE
               NMIN1=2
               IF(IFACE.EQ.2.AND.IWET.EQ.0) NMIN1=4
               IF(NREV.GT.NMIN1) THEN
                  DPDT=DPDTPRE(IDXREV,L)
               ELSE
                  DPDT=0.
               END IF

C --- Use DPOTDT earlier in UNS case ---- Yiran Su 09/26/2017
      IF ((IUNS.EQ.1).AND.(NTSTEP_KEY.GT.30)) DPDT=DPDTPRE(IDXREV,L)

            END IF
            CPBN(N,M) = VINFSB(N,M)-VTOTS(L)-2.0*DPDT

            ! Yiran Su added the effective pressure term 08/18/2018
            CPBN(N,M) = CPBN(N,M) + EPRES(N,M)

            CPB(N,M)=CPBN(N,M)-TWOGY
 90      CONTINUE
 50   CONTINUE

C-----------------------------------------------------------------------
C     Pressure at split panel.                                  JY120598
C-----------------------------------------------------------------------
      IF(IWET.EQ.0) THEN
         DO 100 M=1,MR
            DO 110 II=1,ISR
               IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN
                  IDR=1
                  I=1
               ELSE
                  IDR=2
                  I=-1
               END IF
               SGN=FLOAT(I)
               IF(JCV(M,IDR).GT.0.AND.NSPP(M,IDR).EQ.1) THEN
                  IF(IDR.EQ.1) THEN
                     N=LCV(M,IDR)-1
                     NN1=NC-N
                  ELSE
                     N=M0(M,IDR)-JCV(M,IDR)-1
                     NN1=N-1
                  END IF

                  CPL=CPB(N,M)
                  IF(NN1.GE.3) THEN
                     S2=(SPZ(N+(I*3),M)-SPZ(N+(I*2),M))*SGN
                     S1=(SPZ(N+(I*3),M)-SPZ(N+(I*1),M))*SGN
                     S0=(SPZ(N+(I*3),M)-SPZ(N,M))*SGN
                     CALL QUADEXT(S0,S1,S2,Q1,Q2,Q3)
                     CPR=CPB(N+(I*1),M)*Q1+CPB(N+(I*2),M)*Q2+
     *                    CPB(N+(I*3),M)*Q3
                  ELSE IF(NN1.EQ.2) THEN
                     S1=(SPZ(N+(I*2),M)-SPZ(N+(I*1),M))*SGN
                     S0=(SPZ(N+(I*2),M)-SPZ(N,M))*SGN
                     Q1=S0/S1
                     Q2=ONE-Q1
                     CPR=CPB(N+(I*1),M)*Q1+CPB(N+(I*2),M)*Q2
                  ELSE
                     CPR=CPB(N+(I*1),M)
                  END IF
                  CPB(N,M)=CPL*FLP(M,IDR)+CPR*FRP(M,IDR)
                  IF(-CPB(N,M)*ADVCO2.GT.SIGMA) CPB(N,M)=-SIGMA/ADVCO2

                 CPL=CPBN(N,M)
                  IF(NN1.GE.3) THEN
                     S2=(SPZ(N+(I*3),M)-SPZ(N+(I*2),M))*SGN
                     S1=(SPZ(N+(I*3),M)-SPZ(N+(I*1),M))*SGN
                     S0=(SPZ(N+(I*3),M)-SPZ(N,M))*SGN
                     CALL QUADEXT(S0,S1,S2,Q1,Q2,Q3)
                     CPR=CPBN(N+(I*1),M)*Q1+CPBN(N+(I*2),M)*Q2+
     *                    CPBN(N+(I*3),M)*Q3
                  ELSE IF(NN1.EQ.2) THEN
                     S1=(SPZ(N+(I*2),M)-SPZ(N+(I*1),M))*SGN
                     S0=(SPZ(N+(I*2),M)-SPZ(N,M))*SGN
                     Q1=S0/S1
                     Q2=ONE-Q1
                     CPR=CPBN(N+(I*1),M)*Q1+CPBN(N+(I*2),M)*Q2
                  ELSE
                     CPR=CPBN(N+(I*1),M)
                  END IF
                  CPBN(N,M)=CPL*FLP(M,IDR)+CPR*FRP(M,IDR)
                  IF(-CPBN(N,M)*ADVCO2.GT.SIGMA) CPBN(N,M)=-SIGMA/ADVCO2
               END IF
 110        CONTINUE
 100     CONTINUE
      END IF

C-----------------------------------------------------------------------
C     The old method of extrapolating the pressures with RZPSQ has been
C     replaced with a new 2-D least square extrapolation.       JY080400
C-----------------------------------------------------------------------
C....Same for ICON=8 as well (JY110100)

c      IF(ITUN .EQ. 0 .AND. IDUCT .EQ. 0) THEN
      IF(ITUN .EQ. 0) THEN ! Lets do extrapolation with ducted props.
         IF(MR.NE.MRTIP.AND.IT.NE.2.AND.ICON.NE.5.AND.ICON.NE.8) THEN
            CALL VELEXTRAP
            CALL PRSEXTRAP
         END IF
      ENDIF

C-----------------------------------------------------------------------
C     Calculate difference in pressures at blade T.E.
C-----------------------------------------------------------------------
      DO M=1,MR
         DELCP(M)=CPB(1,M)-CPB(NC,M)
      END DO

      IF(ITUN .NE. 0) GO TO 1111
c      IF(IDUCT .NE. 0) GO TO 1111 ! Lets do extrapolation with ducted props.

      IF(MR.NE.MRTIP) THEN
         MRTIP1=MRTIP+1
         DO M=1,MRTIP
            RSQ1(M)=RZPSQ(M)
         END DO
         RSQ1(MRTIP1)=ONE
         DELCP(MRTIP1)=ZERO
         CALL UGLYDK(MRTIP1,1,1,RSQ1,DELCP,0.0,0.0,DPCUB)
         CALL EVALDK(MRTIP1,MR,RSQ1,RZPSQ,DELCP,DPCUB)
      END IF

 1111 CONTINUE

C-----------------------------------------------------------------------
C     U velocity component along curvilinear coordinate (u,v,w) on the
C       hub by piecewise quadratic potentials
C-----------------------------------------------------------------------

      IF(IHUB.NE.0) THEN
        NHBXM=NHBX-1
        NHBUB=NHBU+NC/2
        DO 180 M=1,MHBT
           DO 190 N=2,NHBUB-1
              L0= INDEXH(N,M)
              L1= INDEXH(N+1,M)
              L2= INDEXH(N-1,M)
              X1=0.5*(DELU(L1)+DELU(L0))
              X2=0.5*(DELU(L2)+DELU(L0))
              DPDUH(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
 190       CONTINUE
           DPDUH(NHBUB,M)=DPDUH(NHBUB-1,M)
     *          +(DPDUH(NHBUB-1,M)-DPDUH(NHBUB-2,M))*X1/X2
           L0= INDEXH(2,M)
           L1= INDEXH(3,M)
           L2= INDEXH(1,M)
           DPDUH(1,M)=DPDUH(2,M)+(DPDUH(2,M)-DPDUH(3,M))
     *              *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))

C..........Calculate U by 1st order differencing after blade............
           DO 200 N=NHBUB+2,NHBXM
              L0= INDEXH(N,M)
              L1= INDEXH(N+1,M)
              L2= INDEXH(N-1,M)
              IF(N.EQ.NHBUB+2) THEN
                 DPDU2= (POT(L2)-POT(L0))/(DELU(L2)+DELU(L0))
              ELSE
                 DPDU2=DPDU1
              END IF
              DPDU1= (POT(L0)-POT(L1))/(DELU(L0)+DELU(L1))
              DPDUH(N,M)=DPDU1+DPDU2
 200       CONTINUE
           DPDUH(NHBX,M)=TWO*DPDU1

C........linear extrapolate u-velocity at t.e...........................
           L0= INDEXH(NHBUB+2,M)
           L1= INDEXH(NHBUB+3,M)
           L2= INDEXH(NHBUB+1,M)
           DPDUH(NHBUB+1,M)=DPDUH(NHBUB+2,M)
     *          +(DPDUH(NHBUB+2,M)-DPDUH(NHBUB+3,M))
     *          *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))
 180    CONTINUE

C-----------------------------------------------------------------------
C        U velocity component along curvilinear coordinate (u,v,w) on
C          the hub 1-st order difference formular
C-----------------------------------------------------------------------

C........Upstream of the blade..........................................
        DO 210 N=1,NHBU
           DO 220 M=1,MHBT
              L0=INDEXH(N,M)
              L1=INDEXH(N,M-1)
              L2=INDEXH(N,M+1)
              IF(M.EQ.1) THEN
                 L1=INDEXH(N,MHBT)
              ELSE IF(M.EQ.MHBT) THEN
                 L2=INDEXH(N,1)
              END IF
              IF(M.EQ.1) THEN
                 DPDV1=(POT(L1)-POT(L0))/(DELV(L1)+DELV(L0))
              ELSE
                 DPDV1=DPDV2
              END IF
              DPDV2=(POT(L0)-POT(L2))/(DELV(L0)+DELV(L2))
              DPDVH(N,M)=DPDV1+DPDV2
 220       CONTINUE
 210    CONTINUE

C......Near the blade & downstream of the blade.........................
        MHBTM=MHBT-1
        DO 230 N=NHBU+1,NHBX
           DO 240 M=2,MHBTM
              L0= INDEXH(N,M)
              L1= INDEXH(N,M-1)
              L2= INDEXH(N,M+1)
              IF(M.EQ.3) THEN
                 DPDV1= (POT(L1)-POT(L0))/(DELV(L0)+DELV(L1))
              ELSE
                 DPDV1=DPDV2
              END IF
              DPDV2= (POT(L0)-POT(L2))/(DELV(L2)+DELV(L0))
              DPDVH(N,M)=DPDV1+DPDV2
 240       CONTINUE
           DPDVH(N,MHBT)=TWO*DPDV2
 230    CONTINUE

C-----------------------------------------------------------------------
C       Calculate the pressure coefficients on the hub
C-----------------------------------------------------------------------
        DO 250 M=1,MHBT
           DO 260 N=1,NHBX
              L=INDEXH(N,M)

C............Perturbation velocities in local coordinate system.........
              UXI=DPDUH(N,M)
              UETA=(DPDVH(N,M)-DPDUH(N,M)*SINPHI(L)) / COSPHI(L)

C............Total velocity components in global coordinate system......
              UXIT=VXIH(N,M)+UXI
              UETAT=VETAH(N,M)+UETA
              UXHTOT(N,M)=UXIT*DIR(L,1,1)+UETAT*DIR(L,2,1)
              UYHTOT(N,M)=UXIT*DIR(L,1,2)+UETAT*DIR(L,2,2)
              UZHTOT(N,M)=UXIT*DIR(L,1,3)+UETAT*DIR(L,2,3)

C............Total velocity square......................................
              VTOTS(L)=(VXIH(N,M)+UXI)**2+(VETAH(N,M)+UETA)**2

C............Graviational term..........................................
              RCP=SQRT(XCT(L,2)*XCT(L,2)+XCT(L,3)*XCT(L,3))
              THETA=THETAT+ACOS(XCT(L,2)/RCP)
              RCP2=RCP*RCP
              IF((NTSTEP.EQ.0).OR.(ISTEADY.EQ.0).OR.(ISTEADY.EQ.3))THEN
                 TWOGY = 0.
CSH---------------shree----------07/08/2002--------
                 IF((ITERGB.EQ.1).OR.(IFLAP.EQ.1))THEN
                    ZTEMP = RCP*COS(THETA)
                    TWOGY=-ZTEMP*2./(FNRUDDER)
                 ENDIF
CSH---------------shree----------07/08/2002--------
              ELSE
                 TWOGY=RCP*COS(THETA)/(FROUDE*ADVCO2)
              ENDIF

C............Time dericavtive of potential..............................
              IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NTSTEP.EQ.0) THEN
                 DPDT=0.0
              ELSE
                 NMIN1=2
                 IF(IFACE.EQ.2.AND.IWET.EQ.0) NMIN1=4
                 IF(NREV.GT.NMIN1) THEN
                    DPDT=DPDTPRE(IDXREV,L)
                 ELSE
                    DPDT=0.
                 END IF

C --- Use DPOTDT earlier in UNS case ---- Yiran Su 09/26/2017
      IF ((IUNS.EQ.1).AND.(NTSTEP_KEY.GT.30)) DPDT=DPDTPRE(IDXREV,L)

              END IF

              CPHN(N,M)=VINFSH(N,M)-VTOTS(L)-2.0*DPDT
              CPH(N,M)=CPHN(N,M)-TWOGY

 260       CONTINUE
 250    CONTINUE
      END IF


C-----------------------------------------------------------------------
C     U velocity component along curvilinear coordinate (u,v,w) on the
C     duct by piecewise quadratic potentials
C-----------------------------------------------------------------------

      IF(IDUCT .NE. 0) THEN
         DO M=1, MDUCT
            DO N=2, NDUCT-1

C --- S.H.CHANG 03/22/2010
C               IF(N.LE.NDUCT/2) THEN
C                 L0= INDEXD(N,M)
C                 L1= INDEXD(N+1,M) !INDEXD(N+1,M) S.H.CHANG 03/21/2010
C                 L2= INDEXD(N-1,M) !INDEXD(N-1,M) S.H.CHANG 03/21/2010
C                 X1=0.5*(DELU(L1)+DELU(L0))
C                 X2=0.5*(DELU(L2)+DELU(L0))
C                 DPDUD(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
C               ELSE IF(N.GT.NDUCT/2) THEN
C                 L0= INDEXD(N,M)
C                 L1= INDEXD(N-1,M) !INDEXD(N+1,M) S.H.CHANG 03/21/2010
C                 L2= INDEXD(N+1,M) !INDEXD(N-1,M) S.H.CHANG 03/21/2010
C                 X1=0.5*(DELU(L1)+DELU(L0))
C                 X2=0.5*(DELU(L2)+DELU(L0))
C                 DPDUD(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
C               END IF
C --- S.H.CHANG 03/22/2010

               L0= INDEXD(N,M)
               L1= INDEXD(N+1,M)
               L2= INDEXD(N-1,M)
               X1=0.5*(DELU(L1)+DELU(L0))
               X2=0.5*(DELU(L2)+DELU(L0))
               DPDUD(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
            END DO

            L0= INDEXD(NDUCT-1,M)
            L1= INDEXD(NDUCT-2,M)
            L2= INDEXD(NDUCT,M)
            DPDUD(NDUCT,M)=DPDUD(NDUCT-1,M)+
     *           (DPDUD(NDUCT-1,M)-DPDUD(NDUCT-2,M))*
     *           (DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))
            L0= INDEXD(2,M)
            L1= INDEXD(3,M)
            L2= INDEXD(1,M)
            DPDUD(1,M)=DPDUD(2,M)+(DPDUD(2,M)-DPDUD(3,M))
     *           *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))
         END DO

         DO N=1, NDUCT
            DO M=2, MDUCT-1
               L0=INDEXD(N,M)
               L1=INDEXD(N,M-1)
               L2=INDEXD(N,M+1)
               X1=0.5*(DELV(L1)+DELV(L0))
               X2=0.5*(DELV(L2)+DELV(L0))
               DPDVD(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
            ENDDO

            L0= INDEXD(N,M)
            L1= INDEXD(N,M-1)
            L2= INDEXD(N,M+1)
            DPDVD(N,MDUCT)=DPDVD(N,MDUCT-1)+
     *           (DPDVD(N,MDUCT-1)-DPDVD(N,MDUCT-2))*
     *           (DELV(L2)+DELV(L0))/(DELV(L1)+DELV(L0))  !modified from X1/X2
            L0= INDEXD(N,2)
            L1= INDEXD(N,3)
            L2= INDEXD(N,1)
            DPDVD(N,1)=DPDVD(N,2)+(DPDVD(N,2)-DPDVD(N,3))
     *           *(DELV(L2)+DELV(L0))/(DELV(L1)+DELV(L0))
C Yiran Added
            IF ((N.GT.NDAFT).and.(N.LE.NDAFT+NH)) THEN
C            if ((n.ge.n1).and.(n.le.n2)) then   !n1 and n2 are the beginning and ending duct chord-wise index of the blade part
              l1=indexd(n,2)
              l0=indexd(n,1)
              dpdvd(n,1)=(pot(l1)-pot(l0))*2.0/(delv(l1)+delv(l0))
              l1=indexd(n,mduct)
              l0=indexd(n,mduct-1)
              dpdvd(n,mduct)=(pot(l1)-pot(l0))*2.0/(delv(l1)+delv(l0))
            end if
C Yiran End

         ENDDO

C Yiran modified the duct pressure evaluation part
         if(ian.ne.1) then !if ian=1 (PSF2), we don't need to do this.
            call duct_pres
         endif
C Finish

C --- S.H.CHANG 03/22/2010
         DO N = 1,NDUCT
            DO M = 1,MDUCT
               DPDVD(N,M) = -DPDVD(N,M)
            END DO
         END DO
C --- S.H.CHANG 03/22/2010

C-----------------------------------------------------------------------
C       Calculate the pressure coefficients on the tunnel
C-----------------------------------------------------------------------

         DO M=1, MDUCT
            DO N=1,NDUCT

               L=INDEXD(N,M)

C............Perturbation velocities in local coordinate system.........

               UXI=DPDUD(N,M)
               UETA=(DPDVD(N,M)-DPDUD(N,M)*SINPHI(L)) / COSPHI(L)

C............Total velocity components in global coordinate system......

               UXIT=VXID(N,M)+UXI
               UETAT=VETAD(N,M)+UETA

               UXDTOT(N,M)=UXIT*DIR(L,1,1)+UETAT*DIR(L,2,1)
               UYDTOT(N,M)=UXIT*DIR(L,1,2)+UETAT*DIR(L,2,2)
               UZDTOT(N,M)=UXIT*DIR(L,1,3)+UETAT*DIR(L,2,3)

C............Total velocity square......................................

               VTOTS(L)=(VXID(N,M)+UXI)**2+(VETAD(N,M)+UETA)**2
               RCP=SQRT(XCT(L,2)*XCT(L,2)+XCT(L,3)*XCT(L,3))
               THETA=THETAT+ACOS(XCT(L,2)/RCP)

               RCP2=RCP*RCP

               TWOGY = 0.

C............Time derivative of potential..............................

               IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NTSTEP.EQ.0) THEN
                  DPDT=0.0
               ELSE
                  NMIN1=2
                 IF(IFACE.EQ.2.AND.IWET.EQ.0) NMIN1=4
                 IF(NREV.GT.NMIN1) THEN
                    DPDT=DPDTPRE(IDXREV,L)
                 ELSE
                    DPDT=0.
                 END IF
              END IF
              CPDN(N,M)=VINFSD(N,M)-VTOTS(L)-2.0*DPDT
              CPD(N,M)=CPDN(N,M)-TWOGY
           ENDDO
        ENDDO

C -- Calculate Mean potential (to cal. mean pressure on duct)

         DO N = 1, NDUCT
            DPOTMEAN(N) = 0.0
            DO M = 1, MDUCT
               I1 = INDEXD(N,M)
               DPOTMEAN(N) = DPOTMEAN(N) + POT(I1)
            ENDDO
            DPOTMEAN(N) = DPOTMEAN(N) / FLOAT(MDUCT)
         ENDDO

         MM1 = MDUCT/2+1

C         DO N=2, NDUCT-1
C            L0= INDEXD(N,MM1)
C            L1= INDEXD(N+1,MM1)
C            L2= INDEXD(N-1,MM1)
C            X1=0.5*(DELU(L1)+DELU(L0))
C            X2=0.5*(DELU(L2)+DELU(L0))
C            DPDUMEAN(N)=DPDXQ(X1,X2,DPOTMEAN(N+1),
C     &                  DPOTMEAN(N),DPOTMEAN(N-1))
C         ENDDO

C         DPDUMEAN(NDUCT)=DPDUMEAN(NDUCT-1)
C     &        +(DPDUMEAN(NDUCT-1)-DPDUMEAN(NDUCT-2))*X1/X2
C         L0= INDEXD(2,MM1)
C         L1= INDEXD(3,MM1)
C         L2= INDEXD(1,MM1)
C         DPDUMEAN(1)=DPDUMEAN(2)+(DPDUMEAN(2)-DPDUMEAN(3))
C     &        *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))

C         DO N = 1,NDUCT
C            L = INDEXD(N,MM1)
C            UETA=-DPDUMEAN(N)*SINPHI(L)/COSPHI(L)
C            UXIT=VXID(N,MM1)+DPDUMEAN(N)
C            UETAT=VETAD(N,MM1)+UETA
C            VTOTMEAN = UXIT**2+UETAT**2
C            DCPMEAN(N) = VINFSD(N,MM1) - VTOTMEAN
C         ENDDO

         DO N = 1,NDUCT
            DPDUMEAN(N) = 0.0
            DPDVMEAN(N) = 0.0
            DO M = 1,MDUCT
               DPDUMEAN(N) = DPDUMEAN(N)+DPDUD(N,M)
               DPDVMEAN(N) = DPDVMEAN(N)+DPDUD(N,M)
            END DO
            DPDUMEAN(N) = DPDUMEAN(N)/FLOAT(MDUCT)
            DPDVMEAN(N) = DPDVMEAN(N)/FLOAT(MDUCT)
         END DO

         DO N = 1,NDUCT
            L = INDEXD(N,MM1)
            UXI = DPDUMEAN(N)
            UETA=(DPDVMEAN(N)-DPDUMEAN(N)*SINPHI(L))/COSPHI(L)
            UXIT=VXID(N,MM1)+UXI
            UETAT=VETAD(N,MM1)+UETA
            VTOTMEAN = UXIT**2+UETAT**2
            DCPMEAN(N) = VINFSD(N,MM1) - VTOTMEAN
         ENDDO

      ENDIF

C      DO I=1,NPANEL
C       WRITE(*,*) 'COSPHI(I)=',COSPHI(I),'SINPHI(I)=',SINPHI(I)
C      ENDDO

C-----------------------------------------------------------------------
C     U velocity component along curvilinear coordinate (u,v,w) on the
C       tunnel by piecewise quadratic potentials
C-----------------------------------------------------------------------

      IF(ITUN .NE. 0) THEN

         NNT1 = 2
         NNT2 = NAXT-NSIDE-1
         NNT3 = NAXT-NSIDE+2
         NNT4 = NAXT-1

         DO M = 1,MTUNEL
            DO N = NNT1,NNT2
               L0= INDEXTN(N,M)
               L1= INDEXTN(N+1,M)
               L2= INDEXTN(N-1,M)
               X1=0.5*(DELU(L1)+DELU(L0))
               X2=0.5*(DELU(L2)+DELU(L0))
               DPDUTN(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
            END DO
            DPDUTN(NNT2+1,M)=DPDUTN(NNT2,M)
     &           +(DPDUTN(NNT2,M)-DPDUTN(NNT2-1,M))*X1/X2
            L0= INDEXTN(2,M)
            L1= INDEXTN(3,M)
            L2= INDEXTN(1,M)
            DPDUTN(1,M)=DPDUTN(2,M)+(DPDUTN(2,M)-DPDUTN(3,M))
     &           *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))

            DO N = NNT3,NNT4
               L0= INDEXTN(N,M)
               L1= INDEXTN(N+1,M)
               L2= INDEXTN(N-1,M)
               X1=0.5*(DELU(L1)+DELU(L0))
               X2=0.5*(DELU(L2)+DELU(L0))
               DPDUTN(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
            END DO
            DPDUTN(NNT4+1,M)=DPDUTN(NNT4,M)
     &           +(DPDUTN(NNT4,M)-DPDUTN(NNT4-1,M))*X1/X2
            L0= INDEXTN(NNT3-2,M)
            L1= INDEXTN(NNT3-3,M)
            L2= INDEXTN(NNT3-1,M)
            DPDUTN(NNT3-1,M)=DPDUTN(NNT3-2,M)+
     &                      (DPDUTN(NNT3-2,M)-DPDUTN(NNT3-3,M))
     &           *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))
         END DO

C --- S.H.CHANG 03/24/2010
C         DO N = 1,NAXT
C            DO M = 1,MTUNEL
C               L0=INDEXTN(N,M)
C               L1=INDEXTN(N,M-1)
C               L2=INDEXTN(N,M+1)
C               IF(M.EQ.1) THEN
C                  L1=INDEXTN(N,MTUNEL)
C               ELSE IF(M.EQ.MTUNEL) THEN
C                  L2=INDEXTN(N,1)
C               END IF
C               IF(M.EQ.1) THEN
C                  DPDV1=(POT(L1)-POT(L0))/(DELV(L1)+DELV(L0))
C               ELSE
C                  DPDV1=DPDV2
C               END IF
C               DPDV2= (POT(L0)-POT(L2))/(DELV(L0)+DELV(L2))
C               DPDVTN(N,M)=DPDV1+DPDV2
C            END DO
C         END DO

         DO N = 1,NAXT
            DO M = 2,MTUNEL-1
               L0=INDEXTN(N,M)
               L1=INDEXTN(N,M-1)
               L2=INDEXTN(N,M+1)
               X1=0.5*(DELV(L1)+DELV(L0))
               X2=0.5*(DELV(L2)+DELV(L0))
               DPDVTN(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
            END DO
            DPDVTN(N,MTUNEL)=DPDVTN(N,MTUNEL-1)
     *           +(DPDVTN(N,MTUNEL-1)-DPDVTN(N,MTUNEL-2))*X2/X1
            L0= INDEXTN(N,2)
            L1= INDEXTN(N,3)
            L2= INDEXTN(N,1)
            DPDVTN(N,1)=DPDVTN(N,2)+(DPDVTN(N,2)-DPDVTN(N,3))
     *           *(DELV(L2)+DELV(L0))/(DELV(L1)+DELV(L0))
         END DO

         DO N = 1,NAXT
            DO M = 1,MTUNEL
               DPDVTN(N,M) = -DPDVTN(N,M)
            END DO
         END DO
C --- S.H.CHANG 03/24/2010

C-----------------------------------------------------------------------
C       Calculate the pressure coefficients on the tunnel
C-----------------------------------------------------------------------

         DO M=1, MTUNEL
            DO N=1,NAXT

               L=INDEXTN(N,M)

C............Perturbation velocities in local coordinate system.........
               UXI=DPDUTN(N,M)
               UETA=(DPDVTN(N,M)-DPDUTN(N,M)*SINPHI(L)) / COSPHI(L)

C............Total velocity components in global coordinate system......
               UXIT=VXITN(N,M)+UXI
               UETAT=VETATN(N,M)+UETA
               UXTNTOT(N,M)=UXIT*DIR(L,1,1)+UETAT*DIR(L,2,1)
               UYTNTOT(N,M)=UXIT*DIR(L,1,2)+UETAT*DIR(L,2,2)
               UZTNTOT(N,M)=UXIT*DIR(L,1,3)+UETAT*DIR(L,2,3)

C............Total velocity square......................................
               VTOTS(L)=(VXITN(N,M)+UXI)**2+(VETATN(N,M)+UETA)**2

C............Graviational term..........................................

               RCP=SQRT(XCT(L,2)*XCT(L,2)+XCT(L,3)*XCT(L,3))
               THETA=THETAT+ACOS(XCT(L,2)/RCP)
               RCP2=RCP*RCP

               TWOGY = 0.0

C............Time dericavtive of potential..............................

               IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NTSTEP.EQ.0) THEN
                  DPDT=0.0
               ELSE
                  NMIN1=2
                 IF(IFACE.EQ.2.AND.IWET.EQ.0) NMIN1=4
                 IF(NREV.GT.NMIN1) THEN
                    DPDT=DPDTPRE(IDXREV,L)
                 ELSE
                    DPDT=0.0
                 END IF
              END IF

              CPTNN(N,M)=VINFSTN(N,M)-VTOTS(L)-2.0*DPDT
              CPTN(N,M)=CPTNN(N,M)-TWOGY

           ENDDO
        ENDDO
      ENDIF

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
c      IF(IAN .EQ. 2) THEN
c         CALL CALVELC
c         CALL CALVELT
c      ENDIF
C/e S.N.KIM | Aug, 2018.

C-----------------------------------------------------------------------
C     Compute the potential differences at the trailing edge
C-----------------------------------------------------------------------

      IF(IT.NE.999) THEN
         DO 270 M=1,MR
            L=INDEXB(1,M)
            DELP(M)=POT(L+NC-1)-POT(L)
 270     CONTINUE

         IF(IDUCT .NE. 0) THEN
            DO M = 1 , MDUCT
               L1 = INDEXD(1,M)
               L2 = INDEXD(NDUCT,M)
               DELPD(M) = POT(L2) - POT(L1)
            ENDDO
         ENDIF

c         IF(ITUN .EQ. 0 .AND. IDUCT .EQ. 0) THEN
         IF(ITUN .EQ. 0) THEN ! Lets do extrapolation with ducted props.
            IF(MR.NE.MRTIP) THEN
               MRTIP1=MRTIP+1
               DO 280 M=1,MRTIP
                  RSQ1(M)=RZPSQ(M)
 280           CONTINUE
               RSQ1(MRTIP1)=ONE
               DELP(MRTIP1)=ZERO
               CALL UGLYDK(MRTIP1,1,1,RSQ1,DELP,0.0,0.0,DPCUB)
               CALL EVALDK(MRTIP1,MR,RSQ1,RZPSQ,DELP,DPCUB)
            ENDIF
         END IF
      ENDIF

!s--YE TIAN 04/13/2012
      if(ITERKT .eq. 0) then
        DELPonW=DELP
        if(iduct.ne.0) DELPonDW=DELPD ! S.N.KIM - 1/30/2017 
      endif
!e--YE TIAN 04/13/2012

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine PRSDIF>>>>>>>>>>>>>>>>>>>>>>>>>
      END


      subroutine duct_pres
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'
      real,dimension(nsw(mrp)):: tpx_n,tpy_n,tpz_n,tpt_n
      real,dimension(nduct+1,mduct+1)::td_n
      integer,dimension(nduct,mduct)::iva_n
      real,dimension(2)::p1_n,p2_n,p3_n,p4_n,pp_n
      logical::IISD

      nd_n=nduct/2
      ntp_n=nsw(mrp)
      nbld_n=nblade-1
C get the tip blade-wake curve
      do n=1,ntp_n
        tpx_n(n)=xw(n,mrp)
        tpy_n(n)=yw(n,mrp)
        tpz_n(n)=zw(n,mrp)
      end do
      do n=1,ntp_n
        tpt_n(n)=atan2(tpz_n(n),tpy_n(n))
      end do

      do n=2,ntp_n
        if ((tpt_n(n)-tpt_n(n-1)).gt.(4.0)) then
          do i=n,ntp_n
            tpt_n(i)=tpt_n(i)-3.14159265359*2.0
          end do
        elseif ((tpt_n(n)-tpt_n(n-1)).lt.(-4.0)) then
          do i=n,ntp_n
            tpt_n(i)=tpt_n(i)+3.14159265359*2.0
          end do
        end if
      end do

C convert into x-theta domain
      do m=1,(mduct+1)
        do n=1,(nd_n+1)
          td_n(n,m)=atan2(zd(n,m),yd(n,m))
        end do
      end do

      do m=2,(mduct+1)
        if ((td_n(1,m)-td_n(1,m-1)).gt.(4.0e0)) then
          do i=m,(mduct+1)
            td_n(1,i)=td_n(1,i)-3.14159265359*2.0
          end do
        elseif ((td_n(1,m)-td_n(1,m-1)).lt.(-4.0e0)) then
          do i=m,(mduct+1)
            td_n(1,i)=td_n(1,i)+3.14159265359*2.0
          end do
        end if
      end do

      do m=1,(mduct+1)
        do n=2,(nd_n+1)
          if ((td_n(n,m)-td_n(n-1,m)).gt.(4.0e0)) then
            do i=n,(nd_n+1)
              td_n(i,m)=td_n(i,m)-3.14159265359*2.0
            end do
          elseif ((td_n(n,m)-td_n(n-1,m)).lt.(-4.0e0)) then
            do i=n,(nd_n+1)
              td_n(i,m)=td_n(i,m)+3.14159265359*2.0
            end do
          end if
        end do
      end do


c      open(399,file='duct_pan_test.plt',status='unknown')
c      write(399,*) ('VARIABLES="X","T"')
c      write(399,*) ('zone T="a"')
c      do i=1,ntp_n
c        write(399,*) tpx_n(i),tpt_n(i)
c      end do
c      write(399,*) 'zone T="b",i=',nd_n+1,'j=',mduct+1
c      do m=1,(mduct+1)
c        do n=1,(nd_n+1)
c          write(399,*) xd(n,m),td_n(n,m)
c        end do
c      end do
c      close(399)


C find intersecting location
      do m=1,mduct
        do n=1,nduct
          iva_n(n,m)=0
        end do
      end do

      do m=1,mduct
        do n=1,nd_n
          p1_n(1)=xd(n,m)
          p2_n(1)=xd(n+1,m)
          p3_n(1)=xd(n+1,m+1)
          p4_n(1)=xd(n,m+1)
          p1_n(2)=td_n(n,m)
          p2_n(2)=td_n(n+1,m)
          p3_n(2)=td_n(n+1,m+1)
          p4_n(2)=td_n(n,m+1)
          do k=1,5
            do i=2,ntp_n
              do j=1,101
               pp_n(1)=(tpx_n(i-1)*real(j-1)+tpx_n(i)*real(101-j))/100.0
               pp_n(2)=(tpt_n(i-1)*real(j-1)+tpt_n(i)*real(101-j))/100.0
               pp_n(2)=pp_n(2)+real(k-3)*DELK
          if (iisd(p1_n,p2_n,p3_n,p4_n,pp_n).eq.(.true.)) iva_n(n,m)=2
              end do
            end do
          end do
        end do
      end do


C find adjacent points
      do m=1,mduct
        m1=m+1
        m2=m-1
        if (m1.gt.mduct) m1=m1-mduct
        if (m2.lt.1) m2=m2+mduct
        do n=1,nd_n
          n1=n+1
          n2=n-1
          if (iva_n(n,m).eq.0) then
            if (iva_n(n,m1).eq.2) iva_n(n,m)=1
            if (iva_n(n,m2).eq.2) iva_n(n,m)=1
            if (n.eq.1) then
              if (iva_n(n1,m).eq.2) iva_n(n,m)=1
            end if
            if (iva_n(n2,m).eq.2) iva_n(n,m)=1
          end if
        end do
      end do

C calculate all dpdu on iva=1 points
      do m=1,mduct
        do n=2,nduct
          if (iva_n(n,m).eq.1) then
            if (iva_n(n-1,m).eq.2) then
              l0=indexd(n,m)
              l1=indexd(n+1,m)
              dpdud(n,m)=-1.0*(pot(l1)-pot(l0))*2.0/(delu(l1)+delu(l0))
            elseif (iva_n(n+1,m).eq.2) then
              l0=indexd(n-1,m)
              l1=indexd(n,m)
              dpdud(n,m)=-1.0*(pot(l1)-pot(l0))*2.0/(delu(l1)+delu(l0))
            end if
          end if
        end do
        if ((iva_n(1,m).eq.1).and.(iva_n(2,m).eq.2)) then
          m1=m+1
          m2=m-1
          if (m1.gt.mduct) m1=m1-mduct
          if (m2.lt.1) m2=m2+mduct
          if (iva_n(1,m1).eq.0) dpdud(n,m)=dpdud(n,m1)
          if (iva_n(1,m2).eq.0) dpdud(n,m)=dpdud(n,m2)
        end if
      end do

C calculate all dpdv on iva=1 points
      do m=1,mduct
        m1=m+1
        m2=m-1
        if (m1.gt.mduct) m1=m1-mduct
        if (m2.lt.1) m2=m2+mduct
        do n=1,nduct
          if (iva_n(n,m).eq.1) then
            if (iva_n(n,m1).eq.2) then
              l0=indexd(n,m2)
              l1=indexd(n,m)
              dpdvd(n,m)=(pot(l1)-pot(l0))*2.0/(delv(l1)+delv(l0))
            elseif (iva_n(n,m2).eq.2) then
              l0=indexd(n,m)
              l1=indexd(n,m1)
              dpdvd(n,m)=(pot(l1)-pot(l0))*2.0/(delv(l1)+delv(l0))
            end if
          end if
        end do
      end do


C calculate all dpdu and dpdv on iva=2 points
      do m=1,mduct
        m1=m+1
        m2=m-1
        if (m1.gt.mduct) m1=m1-mduct
        if (m2.lt.1) m2=m2+mduct
        do n=1,nduct
          n1=n+1
          n2=n-1
          if (iva_n(n,m).eq.2) then
            ict_n=0
            dpdud(n,m)=0.0
            dpdvd(n,m)=0.0
            if (iva_n(n,m1).ne.2) then
              ict_n=ict_n+1
              dpdud(n,m)=dpdud(n,m)+dpdud(n,m1)
              dpdvd(n,m)=dpdvd(n,m)+dpdvd(n,m1)
            end if
            if (iva_n(n,m2).ne.2) then
              ict_n=ict_n+1
              dpdud(n,m)=dpdud(n,m)+dpdud(n,m2)
              dpdvd(n,m)=dpdvd(n,m)+dpdvd(n,m2)
            end if
            if (iva_n(n1,m).ne.2) then
              ict_n=ict_n+1
              dpdud(n,m)=dpdud(n,m)+dpdud(n1,m)
              dpdvd(n,m)=dpdvd(n,m)+dpdvd(n1,m)
            end if
            if (n2.ge.1) then
              if (iva_n(n2,m).ne.2) then
                ict_n=ict_n+1
                dpdud(n,m)=dpdud(n,m)+dpdud(n2,m)
                dpdvd(n,m)=dpdvd(n,m)+dpdvd(n2,m)
              end if
            end if

            if (ict_n.eq.0) then
              iva_n(n,m)=3
            else
              dpdud(n,m)=dpdud(n,m)/ict_n
              dpdvd(n,m)=dpdvd(n,m)/ict_n
            end if

          end if
        end do
      end do

C calculate all dpdu and dpdv on iva=3 points
      do m=1,mduct
        m1=m+1
        m2=m-1
        if (m1.gt.mduct) m1=m1-mduct
        if (m2.lt.1) m2=m2+mduct
        do n=1,nduct
          n1=n+1
          n2=n-1
          if (iva_n(n,m).eq.3) then
            ict_n=0
            dpdud(n,m)=0.0
            dpdvd(n,m)=0.0
            if (iva_n(n,m1).ne.3) then
              ict_n=ict_n+1
              dpdud(n,m)=dpdud(n,m)+dpdud(n,m1)
              dpdvd(n,m)=dpdvd(n,m)+dpdvd(n,m1)
            end if
            if (iva_n(n,m2).ne.3) then
              ict_n=ict_n+1
              dpdud(n,m)=dpdud(n,m)+dpdud(n,m2)
              dpdvd(n,m)=dpdvd(n,m)+dpdvd(n,m2)
            end if
            if (iva_n(n1,m).ne.3) then
              ict_n=ict_n+1
              dpdud(n,m)=dpdud(n,m)+dpdud(n1,m)
              dpdvd(n,m)=dpdvd(n,m)+dpdvd(n1,m)
            end if
            if (n2.ge.1) then
              if (iva_n(n2,m).ne.3) then
                ict_n=ict_n+1
                dpdud(n,m)=dpdud(n,m)+dpdud(n2,m)
                dpdvd(n,m)=dpdvd(n,m)+dpdvd(n2,m)
              end if
            end if

            if (ict_n.eq.0) then
              iva_n(n,m)=4
              write (*,*) "Warning: iva_n=4 in duct_pres subroutine"
            else
              dpdud(n,m)=dpdud(n,m)/ict_n
              dpdvd(n,m)=dpdvd(n,m)/ict_n
            end if

          end if
        end do
      end do

cC output value
c      open(399,file='duct_pan_status.plt',status='unknown')
c      write(399,*) ('VARIABLES="X","Y","Z","stat"')
c      write(399,*) 'ZONE T="duct",I=',NDUCT+1,',J=',MDUCT+1,',DATAPACKI
c     *NG=BLOCK,VARLOCATION=','(4=CELLCENTERED)'
cC      write(399,*) (399,6891) nduct,mduct
c      write(399,*) ((XD(N,M),N=1,NDUCTP),M=1,MDUCTP)
c      write(399,*) ((YD(N,M),N=1,NDUCTP),M=1,MDUCTP)
c      write(399,*) ((ZD(N,M),N=1,NDUCTP),M=1,MDUCTP)
c      write(399,*) ((IVA_N(N,M),N=1,NDUCT),M=1,MDUCT)
c      close(399)

* 6891 FORMAT(1x,'ZONE T="duct",I=',I5,',J=',I5,
*     * ',DATAPACKING=BLOCK,VARLOCATION=(4=CELLCENTERED)')

      end
