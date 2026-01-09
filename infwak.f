      SUBROUTINE INFWAK
************************************************************************
*     INFWAK: INFluence coefficients due to the transition WAKe        *
*      --- Calculate the influence coefficients due to the transition  *
*          wake                                                        *
*                                                                      *
*  Date of last Revision       Revision                                *
*  ---------------------       --------                                *
*     05-14-89   Multi-bladed                                          *
*                The change of tip seperated wake will be made later   *
*     05-24-90   Correect one index error                              *
*     JY040398   I replaced the triangular panels near the blade with  *
*                hyperboloidal panels.                                 *
*                                                                      *
************************************************************************
      use m_WKNP
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
c XM YU 10/2011
      INCLUDE 'PUFCAVC.INC'
c XM YU 10/2011
!     COMMON /WKNP/ NWIDX(MBZ),NWSEC(MBZ)
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)
CSH--Calculate the influence coeff from the image, so NBLADE=2
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2
CSH--------------------------------------------
C-----------------------------------------------------------------------
C     Prepare parameters
C-----------------------------------------------------------------------
C.....NWMIN is the number of wake panels at each radius that the........
C.....vortices shedding.................................................

      NWMIN=999

      IDDUCT2 = NPANB + NPANH
      IDDUCT3 = NPANB + NPANH + NPAND

      IAAA1 = NPANB
      IF(IHUB .NE. 0) IAAA1 = IAAA1+NPANH
      IF(IDUCT .NE. 0) IAAA1 = IAAA1 + NPAND

C      WRITE(*,'(5X,((10I5)))') (NSW(M),M=1,MRP)
c      open(371,file='nsw.dat')
c      do m=1,mrp
c         write(371,*) nsw(m)
c      enddo 
 
      DO 5 M=1,MRP
         NWMIN=MIN(NWMIN,NSW(M))
    5 CONTINUE
      NWMIN=NWMIN-1

C.....Determine NWSEC: number of wake panels in each section............
      DO 10 M=1,MR
        NWSEC(M)=MIN(NSW(M),NSW(M+1))
   10 CONTINUE
      
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
 30     CONTINUE
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
c      open(438,file='nwidx.dat')
c      do m=mr,1,-1
c         write(438,*) nwidx(m)
c      enddo
c      write(*,*) '2', npwake

      CALL GEO3DW(NPWAKE,XGW,CHRLEWS,IER)

      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAK'
         STOP
      END IF

c XM YU add if
      if (ian.ne.6) then
C******  Geom info on the wake surface   by Hong Sun 01/20/2005  *******
      DO 60 J = 1, NPWAKE

        RCPW = SQRT(XCTW(J,2)**2+XCTW(J,3)**2)
        THPW = ATAN2(XCTW(J,3),XCTW(J,2))
        XCTPW(J,1,1)=XCTW(J,1)
        XCTPW(J,2,1)=XCTW(J,2)
        XCTPW(J,3,1)=XCTW(J,3)

        DO K = 1, 3
         DO K1 = 1, 3
           DIRWCP(J,K1,K,1) = DIRW(J,K1,K)     
         END DO
        END DO

      IF(IHULL .EQ. 1 .and. ICON .EQ. 5) THEN
         XCTPW(J,1,2) = XCTW(J,1)
         XCTPW(J,2,2) = 2. - XCTW(J,2)
         XCTPW(J,3,2) = XCTW(J,3)
      ELSE
          IF(NBLADE.GT.1) THEN
           DO 55 KK=2,NBLADE
            XCTPW(J,1,KK)=XCTW(J,1)
            XCTPW(J,2,KK)=RCPW*COS(THPW+DELK*(KK-1))
            XCTPW(J,3,KK)=RCPW*SIN(THPW+DELK*(KK-1))
 55        CONTINUE    
          END IF 
      END IF

 60   CONTINUE
C***********************************************************************
      endif
c XM YU add if

cC------------------------------(S.H.CHANG 02/25/2010)-------------------------
cC    OUTPUT GEOMETRIES IN THE WAKE FOR HULLFPP
c      WRITE(736,*) NPWAKE
c      DO L = 1,NPWAKE
c         DO M = 1,4
c            DO N = 1,3
c               WRITE(736,*) XGW(L,M,N)
c            END DO
c         END DO
c      END DO
c 
c      WRITE(736,*) NWMIN
c      DO I = 1,MR
c         WRITE(736,*) NWIDX(I),NWPAN(I)
c      END DO
c      WRITE(736,'(1X,((40I3)))') (NSW(M),M=1,MRP)
cC------------------------------(S.H.CHANG 02/25/2010)------------------------- 

C-----------------------------------------------------------------------
C        Calculate influence coefficients of the wake.
C        W(I,M):  induced potentials at I due to Mth helical strip of 
C                 dipoles
C-----------------------------------------------------------------------

C.....Initialization....................................................
C.....WSTINF is for the steady problem,WUSINF is the wake inf. functions
C.....wich always have steady soln's strength...........................
      DO 90 M=1,MR
        DO 70 I=1,NPANEL
          WSTINF(I,M)=0.0
          WUSINF(I,M)=0.0
   70   CONTINUE 
 90   CONTINUE

!     DO 290 m=mr,1,-1
!        DO 270 L1=1,nwmin

      DO 290 L1=1, NWMIN
        DO 270 M=MR,1,-1

            L=IDXWAK(L1,M)
C.....Transfer data to the common block /GEOM/..........................
            DO 110 K=1,4
               XV(K)=XVPW(L,K)
               YV(K)=YVPW(L,K)
               SIDE(K)=SIDW(L,K)
110         CONTINUE
            DO 130 K=1,15
               S(K)=SSW(L,K)
130         CONTINUE
C
C.....XM1,XM2,XM3,XM4 for Morino's formulation..........................
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


C.....Initialization STRGTH: wake inf. function at i of blade k.........
            DO 230 KK=1,NBLADE
               DO 150 I=1,NPANEL
                  STRGTH(I,KK)=ZERO
 150           CONTINUE

C.....Add the influence coefficients from all panels....................
               DO 210 I=1,NPANEL
                  
C.....Transfer control points to the local coordinate...................
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 170 K=1,3
                     XLOC=XLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,1,K)
                     YLOC=YLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,3,K)
 170              CONTINUE

C.....Computethe induced potentials.....................................
                  IMR=IMR0
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     *                 FSX,FSY,FDX,FDY,FDZ,0,IMR)

                 
                  
C.....Near field use Morino's formulation...............................
                  IF(IMR.EQ.2) THEN
                     DO 190 IXYZ=1,3
                        XMC(IXYZ)=XCTP(I,IXYZ,KK)
 190                 CONTINUE
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF

c.....if not self-inf. functions, and fd > 6.28, then...................
c.....I will correct to be 0............................................
CJY                  IF(ABS(FD).GT.6.27) THEN
                  IF(ABS(FD).GT.6.28) THEN
C                     WRITE(*,*) 'INFWAK-A',NTSTEP,KK,I,M,L1,IMR,FD
                     FD=0.0
                  END IF
c XM YU IMAGE MODEL
                  if ((icon.eq.5).and.(wingimag.ne.0)) then

                     call rudim1(L,i,kk,imr,xm1,xm2,xm3,xm4,xv12,xv13)
                     fd=fd+xv12
                     fs=fs+xv13
                  endif
c XM YU IMAGE MODEL

                  IF(ITUN .NE. 0 .AND. I .GT. IAAA1) FD = -FD

                  IF(IDUCT .EQ. 1) THEN 
                     IF(IDOPT .EQ. 1) THEN     
                        IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                           FD = 0.0
                        ENDIF
                     ENDIF
                    

                  ENDIF

                  STRGTH(I,KK)=FD
                  WSTINF(I,M)=WSTINF(I,M)+FD
   
 210           CONTINUE
 230        CONTINUE

C.....Write wake inf. functions to files................................
            DO 250 KK=1,NBLADE
               IO=80+KK
               CALL WRITE1(IO,STRGTH(1,KK),NPANEL)
 250        CONTINUE
 270     CONTINUE
 290  CONTINUE

C.....non-triangular panels (i.e. before the wake "tail")...............
      DO 610 M=MR,1,-1
         DO 590 L1=NWMIN+1,NWPAN(M)
            L=IDXWAK(L1,M)
C.....Transfer data to the common block /GEOM/..........................
            DO 490 K=1,4
               XV(K)=XVPW(L,K)
               YV(K)=YVPW(L,K)
               SIDE(K)=SIDW(L,K)
 490        CONTINUE
            DO 510 K=1,15
               S(K)=SSW(L,K)
 510        CONTINUE

C.....Addthe influence coefficients from all panels.....................
CSH--ADD the influence coeff from the image, so NBLADE=2
            DO 570 KK=1,NBLADE
               DO 550 I=1,NPANEL

C.....Transfer control points to the local coordinate...................
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 530 K=1,3
                     XLOC=XLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,1,K)
                     YLOC=YLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,3,K)
 530              CONTINUE
                        
C.....Compute the induced potentials....................................
                  IMR=1
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     *                 FSX,FSY,FDX,FDY,FDZ,0,IMR)

                 IF(ITUN .NE. 0 .AND. I .GT. IAAA1) FD = -FD
   
                 IF(IDUCT .EQ. 1) THEN 
                    IF(IDOPT .EQ. 1) THEN     
                       IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                          FD = 0.0
                       ENDIF

                     ENDIF 
                  ENDIF
                 
                  WSTINF(I,M)=WSTINF(I,M)+FD
                  WUSINF(I,M)=WUSINF(I,M)+FD
 550           CONTINUE
 570        CONTINUE
 590     CONTINUE
 610  CONTINUE
C-----------------------------------------------------------------------
C     Total influence coefficient at 0.7R
C-----------------------------------------------------------------------
      DO I = 1, MR
         IF(RZP(I) .GE. 0.7) THEN
            MPDK = I
            GO TO 1100
         ENDIF
      ENDDO
 1100 CONTINUE

      IF(ICON .EQ. 5) GO TO 1232

      DO 1230 I=1,NPANEL
         IF(IDUCT .EQ. 1 .AND. IDOPT .EQ. 1 
     %        .AND. I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) GO TO 1230 
         WSTINF(I,MPDK)=WSTINF(I,MPDK)- (0.8/RULT/DELK/TANBUW)
     *        * WDK(I)
         WUSINF(I,MPDK)=WUSINF(I,MPDK)- (0.8/RULT/DELK/TANBUW)
     *        * WDK(I)
 1230 CONTINUE
      
 1232 CONTINUE

      WRITE(2,'(5X,'' MPDK='',I3,''    RZP(MPDK)='',F10.3)')
     *      MPDK,RZP(MPDK)

C.....Write the influence functions of the real wake....................
      DO 1240 M=MR,1,-1
        CALL WRITE1(81+NBLADE,WUSINF(1,M),NPANEL)
 1240 CONTINUE

C-----store the wake IC due to all blades (including the effect of the 
C-----far wake) into WINF(I,M).  This is necessary because WSTINF is 
C-----written over later on in INDPOT.F (JY061399)

      DO M = 1, MR
         DO I = 1, NPANEL
            WINF(I,M) = WSTINF(I,M)
         ENDDO
      ENDDO

CSH--REPLACE, NBLADE=1-------------------------
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine INFWAK>>>>>>>>>>>>>>>>>>>>>>>>>
      END
