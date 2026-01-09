      SUBROUTINE INFWAK_DUCT
************************************************************************
*     INFWAK: INFluence coefficients due to the DUCT wake              *
*      --- Calculate the influence coefficients due to the transition  *
*          wake                                                        *
*                                                                      *
*  Date of last Revision       Revision                                *
*  ---------------------       --------                                *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)

      IAAA1 = NPANB+NPANH
      IAAA2 = NPANB + NPANH + NPAND

C-----------------------------------------------------------------------
C
C        Prepare panel geometries of     (N,M+1) 2 *-----* 3 (N+1,M+1)
C          the wake, clockwise viewing             |     |
C          from the suction side of the            |     |
C          blade                           (N,M) 1 *-----* 4 (N+1,M)
C
C-----------------------------------------------------------------------

      DO N=1 , NDWK
         DO M = 1, MDUCT
            L = INDEXWD(N,M)

            XGWD(L,1,1)=XDW(N,M+1)
            XGWD(L,1,2)=YDW(N,M+1)
            XGWD(L,1,3)=ZDW(N,M+1)
            XGWD(L,2,1)=XDW(N,M)
            XGWD(L,2,2)=YDW(N,M)
            XGWD(L,2,3)=ZDW(N,M)
            XGWD(L,3,1)=XDW(N+1,M)
            XGWD(L,3,2)=YDW(N+1,M)
            XGWD(L,3,3)=ZDW(N+1,M)
            XGWD(L,4,1)=XDW(N+1,M+1)
            XGWD(L,4,2)=YDW(N+1,M+1)
            XGWD(L,4,3)=ZDW(N+1,M+1)

         ENDDO
      ENDDO
      
      NPWAKED = MDUCT * NDWK

      CALL GEO3DWD(NPWAKED,XGWD,CHRLEWSD,IER)

      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAK'
         STOP
      END IF

      
C******  Geom info on the wake surface   by Hong Sun 06/04/2006  ******
      DO 60 J = 1, NPWAKED

        RCPWD = SQRT(XCTWD(J,2)**2+XCTWD(J,3)**2)
        THPWD = ATAN2(XCTWD(J,3),XCTWD(J,2))
        XCTPDW(J,1,1)=XCTWD(J,1)
        XCTPDW(J,2,1)=XCTWD(J,2)
        XCTPDW(J,3,1)=XCTWD(J,3)

        DO K = 1, 3
         DO K1 = 1, 3
           DIRDWCP(J,K1,K,1) = DIRWD(J,K1,K)     
         ENDDO
        ENDDO

      IF(IHULL .EQ. 1 .and. ICON .EQ. 5) THEN
         XCTPDW(J,1,2) = XCTWD(J,1)
         XCTPDW(J,2,2) = 2. - XCTWD(J,2)
         XCTPDW(J,3,2) = XCTWD(J,3)
      ELSE
          IF(NBLADE.GT.1) THEN
           DO 55 KK=2,NBLADE
            XCTPDW(J,1,KK)=XCTWD(J,1)
            XCTPDW(J,2,KK)=RCPWD*COS(THPWD+DELK*(KK-1))
            XCTPDW(J,3,KK)=RCPWD*SIN(THPWD+DELK*(KK-1))
 55        CONTINUE    
          ENDIF 
        ENDIF

 60   CONTINUE
 
C*********************************************************************

      DO M=1,MDUCT
         DO I=1,NPANEL
            WSTINFD(I,M)=0.0
         ENDDO
      ENDDO

      DO 290 N = 1 , NDWK
         DO 270 M= 1 , MDUCT
            L = INDEXWD(N,M)

            DO 110 K=1,4
               XV(K)=XVPWD(L,K)
               YV(K)=YVPWD(L,K)
               SIDE(K)=SIDWD(L,K)
110         CONTINUE

            DO 130 K=1,15
               S(K)=SSWD(L,K)
130         CONTINUE

            XM1(1)=XDW(N,M+1)
            XM1(2)=YDW(N,M+1)
            XM1(3)=ZDW(N,M+1)
            XM2(1)=XDW(N,M)
            XM2(2)=YDW(N,M)
            XM2(3)=ZDW(N,M)
            XM3(1)=XDW(N+1,M)
            XM3(2)=YDW(N+1,M)
            XM3(3)=ZDW(N+1,M)
            XM4(1)=XDW(N+1,M+1)             
            XM4(2)=YDW(N+1,M+1)             
            XM4(3)=ZDW(N+1,M+1)

            IMR0=0            

            DO 230 KK=1,NBLADE

              DO I=1,NPANEL
                  STRGTHD(I,KK)=ZERO
               ENDDO

               DO 210 I=1,NPANEL
                  
C.....Transfer control points to the local coordinate...................
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 170 K=1,3
                     XLOC=XLOC+(XCTP(I,K,KK)-XCTWD(L,K))*DIRWD(L,1,K)
                     YLOC=YLOC+(XCTP(I,K,KK)-XCTWD(L,K))*DIRWD(L,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,KK)-XCTWD(L,K))*DIRWD(L,3,K)
 170              CONTINUE

                  IMR=IMR0
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWSD(L),FS,FD,
     *                 FSX,FSY,FDX,FDY,FDZ,0,IMR)
                  
                  IF(IMR.EQ.2) THEN
                     DO 190 IXYZ=1,3
                        XMC(IXYZ)=XCTP(I,IXYZ,KK)
 190                 CONTINUE
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF

                  IF(ABS(FD).GT.6.28) THEN
                     FD=0.0
                  END IF

                  IF(IDOPT .EQ. 1) THEN
                     IF(I .LE. IAAA1 .OR. I .GT. IAAA2) THEN
                        FD = 0.0
                     ENDIF

                  ENDIF

                  STRGTHD(I,KK) = FD

                  WSTINFD(I,M)=WSTINFD(I,M)+FD

C -- Infl. of first duct wake for Unsteady Problem

C                  IF(N .EQ. 1) WUTINFD(I,M) = FD
 
 210           CONTINUE
 230        CONTINUE
            
            DO KK=1,NBLADE
               IO=500+KK
               CALL WRITE1(IO,STRGTHD(1,KK),NPANEL)
            ENDDO

 270     CONTINUE
 290  CONTINUE

      RETURN
      END

      MODULE m_rpan_infwak_duct
        type t_rpan_infwak_duct
          real:: XVrpan(4)
          real:: YVrpan(4) 
          real:: SQrpan(15)
          real:: SIDErpan(4)
        end type t_rpan_infwak_duct
      END MODULE m_rpan_infwak_duct

      MODULE m_hypot_infwak_duct
        type t_hypot_infwak_duct
          real:: XM1hypot(3)
          real:: XM2hypot(3)
          real:: XM3hypot(3)
          real:: XM4hypot(3) 
          real:: XMChypot(3)
        end type t_hypot_infwak_duct
      END MODULE m_hypot_infwak_duct

      SUBROUTINE INFWAK_DUCT_FAST
*************************************************************************
*     INFWAK: INFluence coefficients due to the DUCT wake               *
*      --- Calculate the influence coefficients due to the transition   *
*          wake                                                         *
*                                                                       *
*  Date of last Revision       Revision                                 *
*  ---------------------       --------                                 *
*   21-08-2018 S.N.KIM         Rearranged the order of DO loops with    *
*                              OpenMP parallel code in order to enhance *
*                              computing efficiency of calculating wake *
*                              influenct coeff.                         *
*************************************************************************
      use omp_lib
      use m_rpan_infwak_duct
      use m_hypot_infwak_duct
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      TYPE(t_rpan_infwak_duct):: rpan_geo_data
      TYPE(t_hypot_infwak_duct):: hypot_geo_data
      dimension strgth_tmp(NPANZ,KZ,NDWK,MDUCT)

      IAAA1 = NPANB+NPANH
      IAAA2 = NPANB + NPANH + NPAND

C-----------------------------------------------------------------------
C
C        Prepare panel geometries of     (N,M+1) 2 *-----* 3 (N+1,M+1)
C          the wake, clockwise viewing             |     |
C          from the suction side of the            |     |
C          blade                           (N,M) 1 *-----* 4 (N+1,M)
C
C-----------------------------------------------------------------------

      DO N=1 , NDWK
         DO M = 1, MDUCT
            L = INDEXWD(N,M)

            XGWD(L,1,1)=XDW(N,M+1)
            XGWD(L,1,2)=YDW(N,M+1)
            XGWD(L,1,3)=ZDW(N,M+1)
            XGWD(L,2,1)=XDW(N,M)
            XGWD(L,2,2)=YDW(N,M)
            XGWD(L,2,3)=ZDW(N,M)
            XGWD(L,3,1)=XDW(N+1,M)
            XGWD(L,3,2)=YDW(N+1,M)
            XGWD(L,3,3)=ZDW(N+1,M)
            XGWD(L,4,1)=XDW(N+1,M+1)
            XGWD(L,4,2)=YDW(N+1,M+1)
            XGWD(L,4,3)=ZDW(N+1,M+1)

         ENDDO
      ENDDO
      
      NPWAKED = MDUCT * NDWK

      CALL GEO3DWD(NPWAKED,XGWD,CHRLEWSD,IER)

      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAK'
         STOP
      END IF

      
C******  Geom info on the wake surface   by Hong Sun 06/04/2006  ******
      DO 60 J = 1, NPWAKED

        RCPWD = SQRT(XCTWD(J,2)**2+XCTWD(J,3)**2)
        THPWD = ATAN2(XCTWD(J,3),XCTWD(J,2))
        XCTPDW(J,1,1)=XCTWD(J,1)
        XCTPDW(J,2,1)=XCTWD(J,2)
        XCTPDW(J,3,1)=XCTWD(J,3)

        DO K = 1, 3
         DO K1 = 1, 3
           DIRDWCP(J,K1,K,1) = DIRWD(J,K1,K)     
         ENDDO
        ENDDO

      IF(IHULL .EQ. 1 .and. ICON .EQ. 5) THEN
         XCTPDW(J,1,2) = XCTWD(J,1)
         XCTPDW(J,2,2) = 2. - XCTWD(J,2)
         XCTPDW(J,3,2) = XCTWD(J,3)
      ELSE
          IF(NBLADE.GT.1) THEN
           DO 55 KK=2,NBLADE
            XCTPDW(J,1,KK)=XCTWD(J,1)
            XCTPDW(J,2,KK)=RCPWD*COS(THPWD+DELK*(KK-1))
            XCTPDW(J,3,KK)=RCPWD*SIN(THPWD+DELK*(KK-1))
 55        CONTINUE    
          ENDIF 
        ENDIF

 60   CONTINUE
 
C*********************************************************************

      DO M=1,MDUCT
         DO I=1,NPANEL
            WSTINFD(I,M)=0.0
         ENDDO
      ENDDO

!$omp parallel
!$omp&private(N,M,L,K,rpan_geo_data,hypot_geo_data,CTLT,
!$omp&        IMR0,KK,I,XLOC,YLOC,ZLOC,IMR,FS,FD,
!$omp&        FSX,FSY,FDX,FDY,FDZ,IXYZ)
!$omp&shared(NDWK,MDUCT,XVPWD,YVPWD,SIDWD,SSWD,
!$omp&       XDW,YDW,ZDW,CHRLEWSD,NBLADE,NPANEL,ZERO,
!$omp&       STRGTH_TMP,XCTWD,DIRWD,XCTP,IDOPT,IAAA1,IAAA2,
!$omp&       WSTINFD)
!$omp do
      DO 290 M = 1, MDUCT
         DO 270 N = 1, NDWK
            L = INDEXWD(N,M)

            DO 110 K=1,4
               rpan_geo_data%XVrpan(K)=XVPWD(L,K)
               rpan_geo_data%YVrpan(K)=YVPWD(L,K)
               rpan_geo_data%SIDErpan(K)=SIDWD(L,K)
110         CONTINUE

            DO 130 K=1,15
               rpan_geo_data%SQrpan(K)=SSWD(L,K)
130         CONTINUE
            CTLT = CHRLEWSD(L) 

            hypot_geo_data%XM1hypot(1)=XDW(N,M+1)
            hypot_geo_data%XM1hypot(2)=YDW(N,M+1)
            hypot_geo_data%XM1hypot(3)=ZDW(N,M+1)
            hypot_geo_data%XM2hypot(1)=XDW(N,M)
            hypot_geo_data%XM2hypot(2)=YDW(N,M)
            hypot_geo_data%XM2hypot(3)=ZDW(N,M)
            hypot_geo_data%XM3hypot(1)=XDW(N+1,M)
            hypot_geo_data%XM3hypot(2)=YDW(N+1,M)
            hypot_geo_data%XM3hypot(3)=ZDW(N+1,M)
            hypot_geo_data%XM4hypot(1)=XDW(N+1,M+1)             
            hypot_geo_data%XM4hypot(2)=YDW(N+1,M+1)             
            hypot_geo_data%XM4hypot(3)=ZDW(N+1,M+1)

            IMR0=0            

            DO 230 KK=1,NBLADE

              DO I=1,NPANEL
                  STRGTH_TMP(I,KK,N,M)=ZERO
               ENDDO

               DO 210 I=1,NPANEL
                  
C.....Transfer control points to the local coordinate...................
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 170 K=1,3
                     XLOC=XLOC+(XCTP(I,K,KK)-XCTWD(L,K))*DIRWD(L,1,K)
                     YLOC=YLOC+(XCTP(I,K,KK)-XCTWD(L,K))*DIRWD(L,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,KK)-XCTWD(L,K))*DIRWD(L,3,K)
 170              CONTINUE

                  IMR=IMR0
                  CALL RPAN_PAR4(XLOC,YLOC,ZLOC,CTLT,FS,FD,
     *                 FSX,FSY,FDX,FDY,FDZ,0,IMR,rpan_geo_data)
                  
                  IF(IMR.EQ.2) THEN
                     DO 190 IXYZ=1,3
                        hypot_geo_data%XMChypot(IXYZ)=XCTP(I,IXYZ,KK)
 190                 CONTINUE
                     CALL HYPOT_PAR3(FD,FS,hypot_geo_data)
                  END IF

                  IF(ABS(FD).GT.6.28) THEN
                     FD=0.0
                  END IF

                  IF(IDOPT .EQ. 1) THEN
                     IF(I .LE. IAAA1 .OR. I .GT. IAAA2) THEN
                        FD = 0.0
                     ENDIF
                  ENDIF

                  STRGTH_TMP(I,KK,N,M) = FD
                  WSTINFD(I,M)=WSTINFD(I,M)+FD
C -- Infl. of first duct wake for Unsteady Problem

 210           CONTINUE
 230        CONTINUE
 270     CONTINUE
 290  CONTINUE
!$omp end do
!$omp end parallel

      DO N = 1, NDWK
        DO M = 1, MDUCT
          DO KK= 1, NBLADE
            IO=500+KK
            CALL WRITE1(IO,STRGTH_TMP(1,KK,N,M),NPANEL)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END 


      SUBROUTINE RPAN_PAR4(Xi,Yi,Zi,CHRLENSi, FSs,FDd,FSXi,FSYi,
     %                  FDXi,FDYi,FDZi,ID,IFLAG,rpan_geo_data)
      use m_rpan_infwak_duct
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      type(t_rpan_infwak_duct)::rpan_geo_data
      REAL XI,YI,ZI,CHRLENSI,FSS,FDD,FSXI,FSYI,FDXI,FDYI,FDZI
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
      xv(i) = dble(rpan_geo_data%XVrpan(i))
      yv(i) = dble(rpan_geo_data%YVrpan(i))
      side(i) = dble(rpan_geo_data%SIDErpan(i))
      enddo
      do i = 1 , 15
        s(i) = dble(rpan_geo_data%SQrpan(i))
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
  11  FD=ZERO
      FS=ZERO
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
      FSX=ZERO
      FSY=ZERO
      FDX=ZERO
      FDY=ZERO
      FDZ=ZERO
C----------------------------------------------------------------------
C
C    LOOP FOR CORNER FUNCTIONS
C
C----------------------------------------------------------------------
      DO 13 N=1,4
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
   13 CONTINUE
C----------------------------------------------------------------------
C
C    LOOP FOR SIDE FUNCTIONS AND SUMS OVER FOUR SIDES
C
C----------------------------------------------------------------------
      DO 33 N=1,4
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
  33  CONTINUE
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


C**************************** MIT PSF10 ********************************
C
C                PROPELLER STEADY FLOW ANALYSIS PROGRAM
C                        PANEL METHOD SOLUTION         
C
C                             Version 1.0
C
C       Copyright (c) Massachusetts Institute of Technology 1997
C
C                     Release Date: 1 January 1997
C
C***********************************************************************
      SUBROUTINE HYPOT_PAR3(FDR,FSR,hypot_geo_data)
C**********************************************************************
C     DOUBLE PRECISION
C     Compute the potential due to sources and dipoles based on 
C      Morino's formula
C     -- x1, x2, x3 and x4 should be input in the clock wise direction
C     -- This subroutine CAN NOT calculate the influence functions of 
C        a TRIANGULAR panel
C     
C     10-04-88 C.Y.HSIN @MHL
C     14-08-17 S.N.KIM  @MHL (UT-AUSTIN)
C
      USE m_hypot_infwak_duct
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TYPE(t_hypot_infwak_duct)::hypot_geo_data
      REAL FDR,FSR 
      DIMENSION X1(3),X2(3),X3(3),X4(3),P(3),PC(3),P1(3),P2(3),P3(3),
     *          XCC(3),RC(3),XI0(4),ETA0(4),FDP(4),FSP(4),  
     *          R(3),A1(3),A2(3),U(3),XNC(3),XN(3),RXA1(3),RXA2(3),
     *          NNEG(4)
      DATA XI0 /-1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA ETA0/-1.0D0,1.0D0,1.0D0,-1.0D0/
      DATA ZERO,QUAD,HALF,TWO,FOUR/0.0D0,0.25D0,0.50D0,
     *                                 2.0D0,4.0D0/
C
C.....Transfer to double precision
      DO 10 I=1,3
         X1(I)=DBLE(hypot_geo_data%XM1hypot(I))         
         X2(I)=DBLE(hypot_geo_data%XM2hypot(I))         
         X3(I)=DBLE(hypot_geo_data%XM3hypot(I))         
         X4(I)=DBLE(hypot_geo_data%XM4hypot(I))         
         XCC(I)=DBLE(hypot_geo_data%XMChypot(I))         
 10   CONTINUE
C      PI=3.14159265D0
      PI = DACOS(-1.D0)
      TWOPI=TWO*PI 
      FOURPI=FOUR*PI 
      PI2=HALF*PI
      DO 20 I=1,3
         PC(I)=QUAD*( X1(I)+X2(I)+X3(I)+X4(I) )
         P1(I)=QUAD*( X3(I)+X4(I)-X2(I)-X1(I) )           
         P2(I)=QUAD*( X3(I)-X4(I)+X2(I)-X1(I) )           
         P3(I)=QUAD*( X3(I)-X4(I)-X2(I)+X1(I) )           
 20   CONTINUE  
C
C.....normal vector at the center point
C
      XI=ZERO
      ETA=ZERO
      DO 30 I=1,3
         P(I)=PC(I)+XI*P1(I)+ETA*P2(I)+XI*ETA*P3(I)
         RC(I)=P(I)-XCC(I)            
         A1(I)=P1(I)+ETA*P3(I)
         A2(I)=P2(I)+XI*P3(I)
 30   CONTINUE
      CALL EXPROD(A1,A2,U)
      UL=DSQRT( U(1)*U(1)+U(2)*U(2)+U(3)*U(3) )
      RCL=DSQRT( RC(1)*RC(1)+RC(2)*RC(2)+RC(3)*RC(3) )
      DO 40 I=1,3 
         XNC(I)=U(I)/UL
 40   CONTINUE        
C
C.....Compute Is, Id at four corner points       
      DO 50 J=1,4
         XI=XI0(J)           
         ETA=ETA0(J)
         DO 60 I=1,3
            P(I)=PC(I)+XI*P1(I)+ETA*P2(I)+XI*ETA*P3(I)
            R(I)=P(I)-XCC(I)
            A1(I)=P1(I)+ETA*P3(I)
            A2(I)=P2(I)+XI*P3(I)
 60      CONTINUE
         CALL EXPROD(A1,A2,U)
         UL=DSQRT( U(1)*U(1)+U(2)*U(2)+U(3)*U(3) )
         DO 70 I=1,3 
            XN(I)=U(I)/UL
 70      CONTINUE
         CALL EXPROD(R,A1,RXA1)
         CALL EXPROD(R,A2,RXA2)
         CALL ENPROD(RXA1,RXA2,RA1RA2)
         CALL ENPROD(R,U,RU)
         CALL ENPROD(R,A1,RA1)
         CALL ENPROD(R,A2,RA2)
         CALL ENPROD(RXA1,XNC,RXA1NC)
         CALL ENPROD(RXA2,XNC,RXA2NC)
         CALL ENPROD(RC,XNC,RNC)
         RL=DSQRT( R(1)*R(1)+R(2)*R(2)+R(3)*R(3) )
         RXA1L=DSQRT( RXA1(1)*RXA1(1)+RXA1(2)*RXA1(2)+RXA1(3)*RXA1(3) ) 
         RXA2L=DSQRT( RXA2(1)*RXA2(1)+RXA2(2)*RXA2(2)+RXA2(3)*RXA2(3) ) 
         A1L=DSQRT( A1(1)*A1(1)+A1(2)*A1(2)+A1(3)*A1(3) )        
         A2L=DSQRT( A2(1)*A2(1)+A2(2)*A2(2)+A2(3)*A2(3) )        
         DEN=RL*RU
         FDP(J)=ATAN3( DEN,RA1RA2 ) 
C........Determine the quardrants
         IF(FDP(J).GE.ZERO) THEN
            NNEG(J)=0
         ELSE
            NNEG(J)=1
         END IF
         FS1=RNC*FDP(J)
         FS2=RXA1NC*DASINH(RA1/RXA1L)/A1L
         FS3=RXA2NC*DASINH(RA2/RXA2L)/A2L
         FSP(J)=-FS2+FS3
 50   CONTINUE
C
C     Compute potential of source and dipole
      FD=FDP(1)-FDP(2)+FDP(3)-FDP(4)
      MNEG=NNEG(1)+NNEG(2)+NNEG(3)+NNEG(4)
      IF(MNEG.EQ.1) THEN
         IF(NNEG(1).EQ.1.OR.NNEG(3).EQ.1) THEN
            FD=FD+TWOPI 
         END IF
      ELSE IF(MNEG.EQ.3) THEN
         IF(NNEG(1).EQ.0.OR.NNEG(3).EQ.0) THEN
            FD=FD-TWOPI 
         END IF
      END IF
      IF(FD.LT.-6.28318531) THEN
         FD=FOURPI+FD
      END IF
      FS1=FSP(1)-FSP(2)+FSP(3)-FSP(4)
      FS=-FD*RNC+FSP(1)-FSP(2)+FSP(3)-FSP(4)
      FDR=SNGL(FD)
      FSR=SNGL(FS)
      RETURN
C))))))))))))))))))))))) End of subroutine HYPOT  (((((((((((((((((((((
      END

