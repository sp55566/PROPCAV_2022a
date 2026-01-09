      SUBROUTINE GPANEL2(V,CHRLENS,IER)
C-----------------------------------------------------------------------
C
C  EVALUATES GEOMETRICAL DATA FOR THE QUADRILATERAL GENERATED
C  BY THE COORDINATES OF THE INPUT VERTICES.
C
C  INPUT :
C  -----
C     V(I,J), I=1,..,4  VERTEX INDEX
C             J=1,..,3  REFERENCE x,y,z VERTEX COORDINATES
C
C  OUTPUT:  PANEL GEOMERTICAL DATA STORED IN THE COMMON BLOCK / PANEL /
C  ------   SHARED WITH THE CALLING SUBROUTINE GEOMT.
C
C    Date of last revision           Revision
C    ---------------------         -----------
C    01-30-91 Neal Fine        -Put UL,VL,WL in common /PANEL/
C    03-04-91 Neal Fine        -removed UL,VL,WL from common /PANEL/
C-----------------------------------------------------------------------

      USE GEOCAMB, ONLY : VPP,XCTP,DIRP,S,SIDE,GAUSSX,GAUSSW,DEU,DEV,
     &                    COSPH,SINPH,ULC,VLC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION N1(4),XVXV(4),YVYV(4),XV(4),YV(4),DX(4),DY(4),XP(6,4),
     *          YP(6,4),SI(15),CS(15)
      DIMENSION V(4,3),IX(4),IY(4),GXL(4,3),GSUM(4),XORP(3),IVER(4),
     *          UL(3),VL(3),WL(3)
      DATA QUART / 0.25D0 /, GND / 0.5773503D0 /, IVER / 1, 4, 3, 2 /
      DATA IX / 1, -1, -1, 1 /, IY / -1, -1, 1, 1 /
      DATA CS / 1.D0, 1.D0, 1.5D0, 1.5D0, 3.75D0, 1.D0, 3.0D0, 1.5D0,
     *          7.5D0, 1.5D0, 1.5D0, 3.75D0, 1.5D0, 7.5D0, 3.75D0 /
      DATA HALF/ .5D+00/, ONE/ 1D+00/
      DATA ONE6/.1666667D+00/
      DATA TOL/ 1.0D-06 /, TOLS/ 1D-12 /
      DATA N1/ 2, 3, 4, 1 /
      ULS=0.D0
      VLS=0.D0
      WLS=0.D0

C-----------------------------------------------------------------------
C     EVALUATION OF REFERENCE COORDINATES OF THE ORIGIN OF THE LOCAL
C     COORDINATE SYSTEM.
C-----------------------------------------------------------------------
C---------Ul, Vl, and Wl are unit vectors in the local (s,v,n)----------
C---------coordinate system.  V1-V4 are panel corners----------CM111097-

      DO 10 J=1,3
         V1=V(1,J)
         V2=V(2,J)
         V3=V(3,J)
         V4=V(4,J)
         V12=V1+V2
         V23=V2+V3
         V34=V3+V4
         V41=V4+V1
         XORP(J)=(V12+V34)*QUART
         ULOC=(V12-V34)*HALF
         VLOC=(V23-V41)*HALF
         ULS=ULS+ULOC*ULOC
         VLS=VLS+VLOC*VLOC
         UL(J)=ULOC
         VL(J)=VLOC
10    CONTINUE
      CHRLENS=DMAX1(ULS,VLS)
      ULS=DSQRT(ULS)
      VLS=DSQRT(VLS)
      DEU=ULS
      DEV=VLS
      IF(ULS.LT.TOL.OR.VLS.LT.TOL) GO TO 98
      DO 20 J=1,3
         UL(J)=UL(J)/ULS
         VL(J)=VL(J)/VLS
         VLC(J)=VL(J)
         ULC(J)=UL(J)
20    CONTINUE

C-----------------------------------------------------------------------
C     EVALUATION OF THE DIRECTION COSINES OF LOCAL FRAME
C-----------------------------------------------------------------------

      CALL EXPROD(UL,VL,WL)
      DO 30 J=1,3
         WLOC=WL(J)
         WLS=WLS+WLOC*WLOC
30    CONTINUE
      WLS=DSQRT(WLS)
      IF(WLS.LT.TOL) GO TO 98
      DO 40 J=1,3
         WL(J)=WL(J)/WLS
40    CONTINUE

      CALL EXPROD(WL,UL,VL)

      COSPH=0.D0
      SINPH=0.D0
      DO 50 J=1,3
         DIRP(1,J)=UL(J)
         DIRP(2,J)=VL(J)
         DIRP(3,J)=WL(J)
         COSPH=COSPH+VLC(J)*VL(J)
         SINPH=SINPH+VLC(J)*UL(J)
50    CONTINUE

C-----------------------------------------------------------------------
C     LOCAL COORDINATES OF VERTICES
C-----------------------------------------------------------------------

      DO 80 I=1,4
         IV=IVER(I)
         DO 70 J=1,3
            DSUM=0.0D0
            DO 60 K=1,3
               DDIR=DIRP(J,K)
               DIF=V(I,K)-XORP(K)
               DSUM=DSUM+DDIR*DIF
60          CONTINUE
            VPP(IV,J)=DSUM
70       CONTINUE
80    CONTINUE

C-----------------------------------------------------------------------
C     LOCAL COORDINATES OF GAUSSIAN NODES RELATIVE TO CENTROID - WEIGHTS
C-----------------------------------------------------------------------

      A1=(VPP(4,1)+VPP(1,1))*HALF
      A2=(VPP(4,1)+VPP(3,1))*HALF
      A3=(VPP(4,1)+VPP(2,1))*HALF
      B1=(VPP(4,2)+VPP(1,2))*HALF
      B2=(VPP(4,2)+VPP(3,2))*HALF
      B3=(VPP(4,2)+VPP(2,2))*HALF
      D0=A1*B2-A2*B1
      D1=A1*B3-A3*B1
      D2=A3*B2-A2*B3
      D0I=1.D0/(3.D0*D0)
      UL(1)=(D1*A1+D2*A2)*D0I
      UL(2)=(D1*B1+D2*B2)*D0I
      UL(3)=0.0D0
      DO 90 J=1,4
         RG=IX(J)*GND
         SG=IY(J)*GND
         GAUSSW(J)=D0  +D1*RG+D2*SG
         GXL(J,1)=A1*RG+A2*SG+A3*RG*SG - UL(1)
         GXL(J,2)=B1*RG+B2*SG+B3*RG*SG - UL(2)
         GXL(J,3)= 0.0D0
         DO 85 JJ=1,2
            VPP(J,JJ)=VPP(J,JJ) - UL(JJ)
85       CONTINUE
90    CONTINUE

C-----------------------------------------------------------------------
C     COMPUTE MOMENTS I(M,N) OF PANEL.  I(0,0)=AREA, I(1,0)=XC, I(0,1)=
C     YC. REMAINING MOMENTS I(M,N) ARE NORMALIZED BY THE ARRAY CS(K).
C     THEN I(M,N) IS REPLACED BY ARRAY S(K) WHERE K=N+1+.5*M*(11-M).
C-----------------------------------------------------------------------

      DO 100 N=1,4
         XV(N)=VPP(N,1)
         YV(N)=VPP(N,2)
         XVXV(N)=XV(N)*XV(N)
         YVYV(N)=YV(N)*YV(N)
100   CONTINUE
      DO 110 N=1,15
        S(N)=0.0D0
110   CONTINUE
      DO 120 N=1,4
         DX(N)=XV(N1(N))-XV(N)
         DY(N)=YV(N1(N))-YV(N)
         SIDE(N)=DSQRT(DX(N)*DX(N)+DY(N)*DY(N))
         S(1)=S(1)+HALF*DX(N)*(YV(N)+YV(N1(N)))
         S(2)=S(2)+ONE6*DX(N)*(YVYV(N)+YV(N)*YV(N1(N))+YVYV(N1(N)))
         S(6)=S(6)-ONE6*DY(N)*(XVXV(N)+XV(N)*XV(N1(N))+XVXV(N1(N)))
120   CONTINUE
      IF (S(1).LT.TOLS) THEN
         WRITE(*,*) 'GPANEL: S(1).LT.TOLS'
         STOP
      END IF
      AI=ONE/S(1)
      S(6)=S(6)*AI
      S(2)=S(2)*AI
      DO 140 N=1,4
         XP(1,N)=XV(N)-S(6)
         YP(1,N)=YV(N)-S(2)
         DO 125 M=2,6
            XP(M,N)=XP(1,N)*XP(M-1,N)
            YP(M,N)=YP(1,N)*YP(M-1,N)
125      CONTINUE
140   CONTINUE
      DO 200 N=1,4
         NXT=N1(N)
         DDX=DX(N)*DX(N)
         DDY=DY(N)*DY(N)
         IF (DDX+DDY.LT.1D-12) GO TO 200
         IF (DDY.GT.DDX) GOTO 1

C-----------------------------------------------------------------------
C        FOLLOWING SECTION IS FOR TANGENT LESS THAN ONE
C-----------------------------------------------------------------------

         T=DY(N)/DX(N)
         DO 160 M=2,4
            K=1+.5D0*M*(11-M)
            M1=M+1
            M2=M+2
            SI(K)=(XP(M1,NXT)*YP(1,NXT)-XP(M1,N)*YP(1,N))/(M1)
     *           +T*(XP(M2,N)-XP(M2,NXT))/(M1*M2)
            S(K)=S(K)+SI(K)
            DO 150 NN=1,M
               MM=M-NN
               NN1=NN+1
               MM1=MM+1
               K=NN1+HALF*MM*(11-MM)
               L=NN+HALF*MM1*(10-MM)
               SI(K)=(XP(MM1,NXT)*YP(NN1,NXT)-XP(MM1,N)*YP(NN1,N))
     *               /(MM1*NN1)-T*NN/MM1*SI(L)
               S(K)=S(K)+SI(K)
150         CONTINUE
160      CONTINUE
         GO TO 200

C-----------------------------------------------------------------------
C        FOLLOWING SECTION IS FOR COTANGENT LESS THAN ONE
C-----------------------------------------------------------------------

1        C=DX(N)/DY(N)
         DO 180 NN=2,4
            NN1=NN+1
            NN2=NN+2
            K=NN1
            SI(K)=C/(NN1*NN2)*(YP(NN2,NXT)-YP(NN2,N))
            S(K)=S(K)+SI(K)
            DO 170 M=1,NN
               NM=NN-M
               NM1=NM+1
               NM2=NM+2
               K=NM1+HALF*M*(11-M)
               L=NM2+HALF*(M-1)*(12-M)
               SI(K)=C*((XP(M,NXT)*YP(NM2,NXT)-XP(M,N)*YP(NM2,N))
     1              /(NM1*NM2)-M*SI(L)/NM1)
               S(K)=S(K)+SI(K)
170         CONTINUE
180      CONTINUE

200   CONTINUE

C----------------------------------------------------------------------
C     RENORMALIZE ARRAY OF MOMENTS
C----------------------------------------------------------------------

      DO 210 N=3,15
        S(N)=S(N)*CS(N)
210   CONTINUE

C-----------------------------------------------------------------------
C     EVALUATION OF REFERENCE COORDINATES OF CENTROID AND GAUSSIAN NODES
C-----------------------------------------------------------------------

      DO 250 J=1,3
         DSUM=0.0D0
         DO 220 L2=1,4
            GSUM(L2)=0.0D0
220      CONTINUE
         DO 230 K=1,3
            DDIR=DIRP(K,J)
            DSUM=DSUM+DDIR*UL(K)
            DO 225 L2=1,4
               GSUM(L2)=GSUM(L2)+DDIR*GXL(L2,K)
225         CONTINUE
230      CONTINUE
         XCTP(J)=DSUM+XORP(J)
         DO 240 L2=1,4
            GAUSSX(L2,J)=GSUM(L2)+XCTP(J)
240      CONTINUE
250   CONTINUE
      IER=1
      GO TO 99

98    IER=0
      WRITE(*,*) 'IER = 0 in gpanel'
!Allen Du 06/16/2019 add the option of coupling with the extension
      open(70401,file="ex_check.log")
      close(70401)
      STOP

99    CONTINUE
      RETURN
      END



