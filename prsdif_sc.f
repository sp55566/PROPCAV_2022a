      SUBROUTINE PRSDIF_SC(IT)

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
      DIMENSION ICAV(NBZ),NCAV(NBZ),NDIR(NBZ)

      integer MBZCUB
      MBZCUB = 4*MBZ-4

C....DTN is the nondimensionalized delta(t).............................
      DTN=DELTAT*ADVCO/PI

      NF=1
      NL=NC

C-----------------------------------------------------------------------
C     U velocity component along curvilinear coordinate (u,v,w) on the 
C       blade by piecewise quadratic potentials
C     Notice that DPDXQ here gives -D(phi)/Dx because the 
C       direction of u in RPAN is opposite to the direction of 
C       'x' along the section. ('x' goes from the T.E. and goes
C       thru the lower surface, L.E., upper surface, to the T.E.
C       again
C-----------------------------------------------------------------------
      NCM=NL-1
      NCM2=NL-2
      DO 10 M=1,MR
         DO 20 N=NF+1,NCM
            L0= INDEXB(N,M)
            L2= INDEXB(N-1,M)
            L1= INDEXB(N+1,M)
            X1=HALF*(DELU(L1)+DELU(L0))
            X2=HALF*(DELU(L2)+DELU(L0))
            DPDUB(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))
 20      CONTINUE
         DPDUB(NL,M)=DPDUB(NCM,M)+(DPDUB(NCM,M)-DPDUB(NCM2,M))*X1/X2

         L0=INDEXB(NF+1,M)
         L2=INDEXB(NF,M)
         L1=INDEXB(NF+2,M)
         DPDUB(NF,M)=DPDUB(NF+1,M)+(DPDUB(NF+1,M)-DPDUB(NF+2,M))
     *        *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))
 10   CONTINUE

C-----------------------------------------------------------------------
C     V velocity component along curvilinear coordinate (u,v,w) on the 
C       blade by piecewise quadratic potentials
C-----------------------------------------------------------------------
      MRM=MR-1

      DO 30 N=NF,NL
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

      DO 50 M=1,MR

C.......Flag with cavitating panels (including the split panel w/ ICAV..
         DO 60 N=1,NC
            ICAV(N)=0
 60      CONTINUE

         IF(IWET.EQ.0) THEN
            DO 70 IDR=1,2
               IF(IDR.EQ.1) THEN
                  K=1
                  ISF=0
               ELSE
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

C..........time dericavtive of potential................................
            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.IDXREV.EQ.0) THEN
               DPDT=0.0
            ELSE
               IF(NREV.GT.4) THEN
                  DPDT=DPDTPRE(IDXREV,L)
               ELSE
                  DPDT=0.
               END IF
            END IF

            IF(N.GE.N0(2).AND.N.LT.N0(1)) THEN

C..........Perturbation velocities in local coordinate system...........
               UXI=DPDUB(N,M)
               UETA=(DPDVB(N,M)-DPDUB(N,M)*SINPHI(L)) / COSPHI(L)

C..........Total velocity components in global coordinate system........
               UXIT=VXIB(N,M)+UXI
               UETAT=VETAB(N,M)+UETA

C..........If ITAN=1, the crossflow will be zeroed. (JY053000)
               ITAN=1
c               IF(IWET.EQ.0.AND.ICON.NE.5.AND.ICON.NE.6) THEN
c                  IF(ISP.EQ.1) ITAN=1
c               END IF
            
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
               
            ELSE
               DPDT=ZERO

               IF(N.GE.N0(1)) THEN
                  IDR=1
                  N1=N-N0(1)+1
               ELSE IF(N.LT.N0(2)) THEN
                  IDR=2
                  N1=N0(2)-N
               END IF
               UTS=DPHIDS2(N1,M,IDR)
               UTV=DPHIDV2(N1,M,IDR)
               IF(IDR.EQ.2) THEN
                  UTW=(UTV-UTS*SINPHI(L))/COSPHI(L)
               ELSE IF(IDR.EQ.1) THEN
                  UTW=(UTV+UTS*SINPHI(L))/COSPHI(L)
               END IF
               VTOTS(L)=UTS**2.+UTW**2.
               UTOTU=UTS/SQRT(VTOTS(L))
               IF(IDR.EQ.1) UTOTU=-UTOTU
               UTOTW=UTW/SQRT(VTOTS(L))
               UXIT=SQRT(VTOTS(L))*UTOTU
               UETAT=SQRT(VTOTS(L))*UTOTW

               UT=VOX1(L)*VL(L,1)+VOY1(L)*VL(L,2)+VOZ1(L)*VL(L,3)
               IF(N.LE.NC/2) THEN
                  DPDVB(N,M)=ABS(UXIT)*SINPHI(L)-UT               
               ELSE
                  DPDVB(N,M)=-ABS(UXIT)*SINPHI(L)-UT
               END IF
            END IF

            UXTOT(N,M)=UXIT*DIR(L,1,1)+UETAT*DIR(L,2,1)
            UYTOT(N,M)=UXIT*DIR(L,1,2)+UETAT*DIR(L,2,2)
            UZTOT(N,M)=UXIT*DIR(L,1,3)+UETAT*DIR(L,2,3)

C.........graviational term.............................................
            RCP=SQRT(XCT(L,2)*XCT(L,2)+XCT(L,3)*XCT(L,3))
            THETA=THETAT+ACOS(XCT(L,2)/RCP)
            RCP2=RCP*RCP
            IF((ISTEADY.EQ.0).OR.(ISTEADY.EQ.3))THEN
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

            CPBN(N,M)=VINFSB(N,M)-VTOTS(L)-2.0*DPDT
            CPB(N,M)=CPBN(N,M)-TWOGY

 90      CONTINUE
 50   CONTINUE

C-----------------------------------------------------------------------
C     Pressure at split panel.                                  JY120598
C-----------------------------------------------------------------------
      IF(IWET.EQ.0) THEN
         DO 100 M=1,MR
            DO 110 IDR=1,2
               IF(IDR.EQ.1) THEN
                  I=1
               ELSE
                  I=-1
               END IF
               SGN=FLOAT(I)
               IF(JCV(M,IDR).GT.0.AND.NSPP(M,IDR).EQ.1) THEN
                  IF(IDR.EQ.1) THEN
                     N=LCV(M,IDR)-1
                     NN1=N0(1)-1-N
                  ELSE
                     N=M0(M,IDR)-JCV(M,IDR)-1
                     NN1=N-N0(2)
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
      IF(MR.NE.MRTIP.AND.IT.NE.2.AND.ICON.NE.5.AND.ICON.NE.8) THEN
C         CALL VELEXTRAP
C         CALL PRSEXTRAP
         DO M=MRTIP+1,MR
            DO N=1,NC
               CPB(N,M)=CPB(N,M-1)
            END DO
         END DO
      END IF

C-----------------------------------------------------------------------
C     Calculate difference in pressures at blade T.E.       
C-----------------------------------------------------------------------
      DO M=1,MR
         DELCP(M)=CPB(1,M)-CPB(NC,M)
      END DO
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
              DPDV2= (POT(L0)-POT(L2))/(DELV(L0)+DELV(L2))
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
              IF((ISTEADY.EQ.0).OR.(ISTEADY.EQ.3))THEN
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
              IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1) THEN
                 DPDT=0.0
              ELSE
                 IF(NREV.GT.4) THEN
                    DPDT=DPDTPRE(IDXREV,L)
                 ELSE
                    DPDT=0.
                 END IF
              END IF
              
              CPHN(N,M)=VINFSH(N,M)-VTOTS(L)-2.0*DPDT
              CPH(N,M)=CPHN(N,M)-TWOGY

 260       CONTINUE
 250    CONTINUE
      END IF

c      if(ian .eq. 2) then
c         call calvelc
c         call calvelt
c      endif

C-----------------------------------------------------------------------
C     Compute the potential differences at the trailing edge
C-----------------------------------------------------------------------
      IF(IT.NE.999) THEN
         DO 270 M=1,MR
            L=INDEXB(1,M)
            DELP(M)=POT(L+NC-1)-POT(L)
 270     CONTINUE
         IF(MR.NE.MRTIP) THEN
            MRTIP1=MRTIP+1
            DO 280 M=1,MRTIP
               RSQ1(M)=RZPSQ(M)
 280        CONTINUE
            RSQ1(MRTIP1)=ONE
            DELP(MRTIP1)=ZERO
            CALL UGLYDK(MRTIP1,1,1,RSQ1,DELP,0.0,0.0,DPCUB)
            CALL EVALDK(MRTIP1,MR,RSQ1,RZPSQ,DELP,DPCUB)
         END IF
      END IF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine PRSDIF>>>>>>>>>>>>>>>>>>>>>>>>>
      END



