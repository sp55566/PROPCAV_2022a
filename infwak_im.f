      SUBROUTINE INFWAK_IM
************************************************************************
*     INFWAK: INFluence coefficients due to the transition WAKe IMages *
*      --- Calculate the influence coefficients due to the transition  *
*          wake images.                                                *
*                                                                      *
*     Date      Comments                                               *
*     --------- -----------------------------------------------------  *
************************************************************************

      use m_WKNP
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

!     COMMON /WKNP/ NWIDX(MBZ),NWSEC(MBZ)
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)

C-----------------------------------------------------------------------
C     Prepare parameters
C-----------------------------------------------------------------------
C.....NWMIN is the number of wake panels at each radius that the........
C.....vortices shedding.................................................
      NWMIN=NWMINFW

C.....Calculate geometric characteristics...............................
      NMIN=0
      NWIDX0=0
      DO 50 M=MR,1,-1

C-----------------------------------------------------------------------
C
C        Prepare panel geometries of     (N,M+1) 2 *-----* 3 (N+1,M+1)
C          the wake, clockwise viewing             |     |
C          from the suction side of the            |     |
C          blade                           (N,M) 1 *-----* 4 (N+1,M)
C
C-----------------------------------------------------------------------

C........NWIDX0 is the beginning panel index at section M, NMIN here....
C........is the total number of wake panels of last section.............
         NWIDX0=NWIDX0+NMIN
         NWIDX(M)=NWIDX0

C........NMIN: initial number of wake panels of this section............
         NMIN=NWSEC(M)

C........Use triangular panels in the wake near the trailing edge.......
         DO 30 N=1,NMIN-1
            L=N+NWIDX0
            XGW(L,1,1)=XW(N,M)
            XGW(L,1,2)=YW(N,M)
            XGW(L,1,3)=ZW(N,M)
            XGW(L,2,1)=XW(N,M+1)
            XGW(L,2,2)=YW(N,M+1)
            XGW(L,2,3)=ZW(N,M+1)
            XGW(L,3,1)=XW(N+1,M+1)
            XGW(L,3,2)=YW(N+1,M+1)
            XGW(L,3,3)=ZW(N+1,M+1)
            XGW(L,4,1)=XW(N+1,M)
            XGW(L,4,2)=YW(N+1,M)
            XGW(L,4,3)=ZW(N+1,M)
 30      CONTINUE
         IF(NSW(M).EQ.NSW(M+1)) THEN
            NMIN=NMIN-1
         ELSE 
            IF(NMIN.EQ.NSW(M)) THEN
               NMIN=NMIN
               NMIN1=NMIN+NWIDX0
               XGW(NMIN1,1,1)=XW(NSW(M),M)
               XGW(NMIN1,1,2)=YW(NSW(M),M)
               XGW(NMIN1,1,3)=ZW(NSW(M),M)
               XGW(NMIN1,2,1)=XW(NSW(M),M+1)
               XGW(NMIN1,2,2)=YW(NSW(M),M+1)
               XGW(NMIN1,2,3)=ZW(NSW(M),M+1)
               XGW(NMIN1,3,1)=XW(NSW(M+1),M+1)
               XGW(NMIN1,3,2)=YW(NSW(M+1),M+1)
               XGW(NMIN1,3,3)=ZW(NSW(M+1),M+1)
               XGW(NMIN1,4,1)=XW(NSW(M),M)
               XGW(NMIN1,4,2)=YW(NSW(M),M)
               XGW(NMIN1,4,3)=ZW(NSW(M),M)
            ELSE
               NMIN=NMIN
               NMIN1=NMIN+NWIDX0
               XGW(NMIN1,1,1)=XW(NSW(M+1),M)
               XGW(NMIN1,1,2)=YW(NSW(M+1),M)
               XGW(NMIN1,1,3)=ZW(NSW(M+1),M)
               XGW(NMIN1,2,1)=XW(NSW(M+1),M+1)
               XGW(NMIN1,2,2)=YW(NSW(M+1),M+1)
               XGW(NMIN1,2,3)=ZW(NSW(M+1),M+1)
               XGW(NMIN1,3,1)=XW(NSW(M),M)
               XGW(NMIN1,3,2)=YW(NSW(M),M)
               XGW(NMIN1,3,3)=ZW(NSW(M),M)
               XGW(NMIN1,4,1)=XW(NSW(M+1),M)
               XGW(NMIN1,4,2)=YW(NSW(M+1),M)
               XGW(NMIN1,4,3)=ZW(NSW(M+1),M)
            END IF
         END IF
         NWPAN(M)=NMIN
 50   CONTINUE

C.....NPWAKE is the total number of panels in the wake..................
      NPWAKE=NWIDX0+NMIN
      CALL GEO3DW(NPWAKE,XGW,CHRLEWS,IER)
      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAK_IM'
         STOP
      END IF

C-----------------------------------------------------------------------
C     Calculate influence coefficients of the wake images.
C-----------------------------------------------------------------------
      DO KK=1,NBLADE
         IO=180+KK
         REWIND IO
      END DO

      DO 290 L1=1,NWMIN
         DO 270 M=MR,1,-1
            L=IDXWAK(L1,M)

C..........Transfer data to the common block /GEOM/.....................
            DO 110 K=1,4
               XV(K)=XVPW(L,K)
               YV(K)=YVPW(L,K)
               SIDE(K)=SIDW(L,K)
 110        CONTINUE
            DO 130 K=1,15
               S(K)=SSW(L,K)
 130        CONTINUE

C..........XM1,XM2,XM3,XM4 for Morino's formulation.....................
            XM1(1)=XW(L1,M)
            XM1(2)=YW(L1,M)
            XM1(3)=ZW(L1,M)
            XM2(1)=XW(L1,M+1)
            XM2(2)=YW(L1,M+1)
            XM2(3)=ZW(L1,M+1)
            XM3(1)=XW(L1+1,M+1)
            XM3(2)=YW(L1+1,M+1)
            XM3(3)=ZW(L1+1,M+1)
            XM4(1)=XW(L1+1,M)             
            XM4(2)=YW(L1+1,M)             
            XM4(3)=ZW(L1+1,M) 
            IMR0=0            

            DO 230 KK=1,NBLADE

               IREC1=NTPOS(KK)

               DO 210 I=1,NPANEL

                  IF(ISUBM(I,IDXREV).EQ.1.AND.L1.LE.MSW(M,IREC1)) THEN

C................Transfer control points to the local coordinate........
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 170 K=1,3
                     XLOC=XLOC+(XCTP_IM(I,K,KK)-XCTW(L,K))*DIRW(L,1,K)
                     YLOC=YLOC+(XCTP_IM(I,K,KK)-XCTW(L,K))*DIRW(L,2,K)
                     ZLOC=ZLOC+(XCTP_IM(I,K,KK)-XCTW(L,K))*DIRW(L,3,K)
 170              CONTINUE

C................Computethe induced potentials..........................
                  IMR=IMR0
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     *                 FSX,FSY,FDX,FDY,FDZ,0,IMR)
                  
C................Near field use Morino's formulation....................
                  IF(IMR.EQ.2) THEN
                     DO 190 IXYZ=1,3
                        XMC(IXYZ)=XCTP_IM(I,IXYZ,KK)
 190                 CONTINUE
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF
                
                  IF(ABS(FD).GT.6.28) THEN
C                     WRITE(*,*) 'INFWAK_IM-A',NTSTEP,KK,I,M,L1,IMR,FD
                     FD=0.0
                  END IF

                  ELSE
                  FD=ZERO
                  END IF

                  STRGTH(I,KK)=FD
 210           CONTINUE
 230        CONTINUE

C..........Write wake inf. functions to files...........................
            DO 250 KK=1,NBLADE
               IO=180+KK
               CALL WRITE1(IO,STRGTH(1,KK),NPANEL)
 250        CONTINUE
 270     CONTINUE
 290  CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine INFWAK_IM>>>>>>>>>>>>>>>>>>>>>>
      END
