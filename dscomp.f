      SUBROUTINE DSCOMP
************************************************************************
*  COMPute the arclength DS of each panel by discretizing each strip   *
*  with nodes for each original node and then adding the lengths of    *
*  each resulting sub-panel. Uses subroutine GBLADE as a base.         *
*  06-26-92 NF                                                         *
*  moved to PUFCAV from PSFCAV on 06-26-92 NF                          *
*                                                                      *
*  Date     Revision or Comment                                        *
*  -------- -----------------------                                    *
*  CM061297 Added NACA16 Thickness Distribution                        *
*  CM110797 Working on a read in thickness distribution                *
*  JY010798 Added a parabolic camber and thickness option for          *
*           blade geometry definition.  Here are the new               *
*           additions:                                                 *
*           ITHK = 10   Parabolic thickness form                       *
*           ICAM = 0    NACA A=0.8 mean line                           *
*           ICAM = 1    Parabolic camber line                          *
*  JY112198 Added a new camber and a new thickness option:             *
*                     ICAM=99  User input camber distribution          *
*                     ITHK=99  User input thickness distribution       *
*                                                                      *
************************************************************************
      use m_INPGEO2
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'INPGEO.INC'

      COMMON /IFORT2/ NPOINTS
      COMMON /MODGEO1/ DUMXC(MXINP),PC(16),
     %                 YP(MXINP,NXMAX),YS(MXINP,NXMAX)
      DIMENSION SBS(10*NBHZ),YTCS(10*NBHZ),YCCS(10*NHBZ),TCCS(10*NHBZ)
     *    ,     XIS(10*NBPZ),ETAS(10*NBPZ),XBS(10*NBPZ,MBPZ)
     *    ,     YBS(10*NBPZ,MBPZ),ZBS(10*NBPZ,MBPZ)
C      DIMENSION SPY1(10*NHBZ,NXMAX),SSY1(10*NHBZ,NXMAX)
C      DIMENSION TMPP1(NXMAX),TMPP2(NXMAX)
C      DIMENSION SPY11(MBPZ,10*NHBZ),SSY11(MBPZ,10*NHBZ)
C      DIMENSION CUBIC(MAXSPL),CUBIC2(MBPZ*4)

!s-- YE TIAN --- 07/07/2013----
      real vtemp(15)
!e-- YE TIAN --- 07/07/2013----

      NCS=10*NC
      NHS=10*NH
      NHSP=NHS+1
      NCSP=NCS+1

C-----------------------------------------------------------------------
C     Set up chordwise spacing
C-----------------------------------------------------------------------

      CALL SPACE(NCS,ICSPAC,RLET,SBS)

C-----------------------------------------------------------------------
C     Generate blade geometry
C-----------------------------------------------------------------------

      DO 50 M=1,MRP

C.......Midchord centered local coordinate xi and eta (nondimen. by R)

C.......CHORD=C/D,THK=T/D,CAMBR=F/C,YTC=Y/D,YCC=Y/C,SB=S/C,RAKE=XM/D

C-----------------------------------------------------------------------
C     ITHK=99, user input thickness distribution.               JY112798
C-----------------------------------------------------------------------
         IF(ITHK.EQ.99) THEN
            CALL THKINP(NHS,XTS1(1,M),SBS,YTCS)
            GO TO 1000
         END IF

C-----------------------------------------------------------------------
C     ITHK=0 option is no longer valid.                         JY051300
C-----------------------------------------------------------------------

         IF(THICK(M).EQ.0.0) THEN
            DO N=1,NHS
               YTCS(N)=ZERO
            END DO
         ELSE
            IF(ITHK.EQ.1) THEN
               CALL NACA66(NHS,THICK(M),RLE,SBS,YTCS)
            ELSE IF(ITHK.EQ.2) THEN
               CALL RAE(NHS,THICK(M),SBS,YTCS)
            ELSE IF(ITHK.EQ.3) THEN
               CALL NACA65A(NHS,THICK(M),RLE,SBS,YTCS)
            ELSE IF(ITHK.EQ.4) THEN 
               CALL NACA64A(NHS,THICK(M),RLE,SBS,YTCS)
            ELSE IF(ITHK.EQ.5) THEN
               CALL NACA00(NHS,THICK(M),RLE,SBS,YTCS)
            ELSE IF(ITHK.EQ.6) THEN
               CALL ELIPSE(NHS,THICK(M),RLE,SBS,YTCS)
C--------------------------------------------------------------------------
C ITHK=8, NAVSEA TYPE_ONE (NOTE: NACA16 is NAVSEA TYPE_TWO) SNKIM
C 08232017
C--------------------------------------------------------------------------
            ELSE IF(ITHK.EQ.8) THEN
               CALL NAVSEA1(NHS,THICK(M),RLE,SBS,YTCS)
            ELSE IF(ITHK.EQ.9) THEN
               CALL NACA16(NHS,THICK(M),RLE,SBS,YTCS)
            ELSE IF(ITHK.EQ.10) THEN
               CALL TKPB(NHS,THICK(M),RLE,SBS,YTCS)
            END IF
         END IF

 1000    CONTINUE

C
C.......0.8 camber mean line
         IF(ICAM.EQ.0) THEN
            CALL A8ML(NHS,CAMBR(M),SBS,YCCS,TCCS)
C.......parabolic camber mean line (jy 010798)
         ELSE IF(ICAM.EQ.1) THEN
            CALL CBPB(NHS,CAMBR(M),SBS,YCCS,TCCS)
         ELSE IF(ICAM.EQ.2) THEN                    !SNKIM082317
            CALL NACA65(NHS,CAMBR(M),SBS,YCCS,TCCS)
C-----------------------------------------------------------------------
C       ICAM=99, user input thickness distribution.            JY112798
C-----------------------------------------------------------------------
         ELSE IF(ICAM.EQ.99) THEN
            CALL CAMINP(NHS,XCS1(1,M),SBS,YCCS,TCCS)
         END IF

         IF(ICON .EQ. 5) GO TO 1100

C.......local coordinate xi and eta
         XIS(NHSP)=-CHORD(M)
         ETAS(NHSP)=ZERO
         DO 20 N=1,NHS
            IF(ITHK .EQ. 99 .AND. ICAM .EQ. 99 
     %           .AND. IFORMAT .EQ. 2) THEN
               DXT = 0.0
               DYT = YTCS(N)
            ELSE
               DXT=YTCS(N)*SIN(TCCS(N))
               DYT=YTCS(N)*COS(TCCS(N))
            ENDIF

            NBOT=NHSP-N
            NTOP=NHSP+N
            XIS(NTOP)=((SBS(N)-.5)*CHORD(M)-DXT)*TWO
            ETAS(NTOP)=(YCCS(N)*CHORD(M)+DYT)*TWO
            XIS(NBOT)=((SBS(N)-.5)*CHORD(M)+DXT)*TWO
            ETAS(NBOT)=(YCCS(N)*CHORD(M)-DYT)*TWO
 20    CONTINUE

C.......Reference coordinate of blade element
       IF(RZ(M).EQ.0.)THEN
          RZ(M)=1
       ENDIF
       TANP=PITCH(M)/PI/RZ(M)
       COSP=1.0/SQRT(1+TANP*TANP)
       SINP=TANP*COSP
       DO 30 N=1,NCSP
          DX= XIS(N)*SINP-ETAS(N)*COSP
          XBS(N,M)=RAKE(M)*TWO +DX
          THETA=SKEW(M)*RAD+(XIS(N)*COSP+ETAS(N)*SINP)/RZ(M)
          YBS(N,M)=RZ(M)*COS(THETA)
          ZBS(N,M)=RZ(M)*SIN(THETA)
 30    CONTINUE

       GO TO 50

C -- For Hydrofoil case.

 1100  CONTINUE

       COST=COS(PITCH(M))
       SINT=SIN(PITCH(M))
       XIS(NHSP)=-CHORD(M) * COST
       ETAS(NHSP)=CHORD(M) * SINT

       DO N = 1 , NHS
          NBOT=NHSP-N
          NTOP=NHSP+N
          NMEAN=N+1
          
          XXX = (SBS(N) - 0.5) * CHORD(M) * TWO
          YYY = YCCS(N) * CHORD(M) * TWO

               
          IF(ITHK .EQ. 99 .AND. ICAM .EQ. 99 
     %         .AND. IFORMAT .EQ. 2) THEN
             DXT = 0.0
             DYT = YTCS(N)
          ELSE
             DXT = YTCS(N) * SIN(TCCS(N)) * TWO
             DYT = YTCS(N) * COS(TCCS(N)) * TWO
          ENDIF

          XU = XXX - DXT
          XL = XXX + DXT
          
          YU = YYY + DYT
          YL = YYY - DYT
          
          XIS(NTOP) = XU * COST + YU * SINT
          ETAS(NTOP) = -XU * SINT + YU * COST
          
          XIS(NBOT) = XL * COST + YL * SINT
          ETAS(NBOT) = -XL * SINT + YL * COST
       ENDDO
       
       DO N=1,NCSP
          XBS(N,M)=RAKE(M) * TWO + XIS(N)  
          YBS(N,M)=RZ(M)
          ZBS(N,M)=ETAS(N) + SKEW(M) * TWO
       ENDDO

 50   CONTINUE

C-----------------------------------------------------------------------
C     compute the arclength on the blade, ds
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      DO 80 M=1,MR
         I1=0
         DO 70 N=1,NC
            SUM=ZERO
            DO 60 I=1,10
               I1=I1+1
               DU=HALF*(XBS(I1+1,M)+XBS(I1+1,M+1)-XBS(I1,M)-XBS(I1,M+1))
               DV=HALF*(YBS(I1+1,M)+YBS(I1+1,M+1)-YBS(I1,M)-YBS(I1,M+1))
               DW=HALF*(ZBS(I1+1,M)+ZBS(I1+1,M+1)-ZBS(I1,M)-ZBS(I1,M+1))
               SUM=SUM+SQRT(DU*DU+DV*DV+DW*DW)
 60         CONTINUE
            DS(N,M)=SUM
 70      CONTINUE
 80   CONTINUE
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine DSCOMP>>>>>>>>>>>>>>>>>>>>>>>>>
      END

      SUBROUTINE DSCOMP_IDRT
************************************************************************
*  COMPute the arclength DS of each panel by discretizing each strip   *
*  with nodes for each original node and then adding the lengths of    *
*  each resulting sub-panel. This subroutine is special version of     *
*  DSCOMP when IDRT=1.                                                 *
*                                                                      *
*  10-15-2018 (birth day of M.LEE)  S.KIM                              *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'INPGEO.INC'

      COMMON /IFORT2/ NPOINTS
      DIMENSION XBS_TMP(10*NC+1),YBS_TMP(10*NC+1),ZBS_TMP(10*NC+1)
      DIMENSION XBS(10*NC+1,MRP),YBS(10*NC+1,MRP),ZBS(10*NC+1,MRP)
      DIMENSION SPLINEL(NCP),SPLINELS(10*NC+1)
      DIMENSION P1(3),P2(3),XP3(NCP),YP3(NCP),ZP3(NCP)
      REAL,DIMENSION(4*(NCP-1)) :: CUBXB,CUBYB,CUBZB

      NCS=10*NC
      NCSP=NCS+1

C-----------------------------------------------------------------------
C     interpolate blade panels into 10 chordwise subdivisions. 
C-----------------------------------------------------------------------
      DO M = 1, MRP
        I1 = 0
        SPLINEL(1) = ZERO
        SPLINELS(1) = ZERO
        DO N = 1,NC
          P1(1) = XB(N,M)
          P1(2) = YB(N,M)
          P1(3) = ZB(N,M)
          P2(1) = XB(N+1,M)
          P2(2) = YB(N+1,M)
          P2(3) = ZB(N+1,M)
          DELL = SQRT(DOT_PRODUCT(P2-P1,P2-P1))
          SPLINEL(N+1) = SPLINEL(N) + DELL
          DELL2 = DELL/10
          DO I = 1, 10
            I1 = I1 + 1
            SPLINELS(I1+1) = SPLINELS(I1) + DELL2
          ENDDO
        ENDDO
        DO N = 1, NCP
          XP3(N) = XB(N,M)
          YP3(N) = YB(N,M)
          ZP3(N) = ZB(N,M)
        ENDDO
        CALL UGLYDK(NCP,1,1,SPLINEL,XP3,0.0,0.0,CUBXB)
        CALL UGLYDK(NCP,1,1,SPLINEL,YP3,0.0,0.0,CUBYB)
        CALL UGLYDK(NCP,1,1,SPLINEL,ZP3,0.0,0.0,CUBZB)
        CALL EVALDK(NCP,NCSP,SPLINEL,SPLINELS,XBS_TMP,CUBXB)
        CALL EVALDK(NCP,NCSP,SPLINEL,SPLINELS,YBS_TMP,CUBYB)
        CALL EVALDK(NCP,NCSP,SPLINEL,SPLINELS,ZBS_TMP,CUBZB)
        DO N = 1, NCSP
          XBS(N,M) = XBS_TMP(N)
          YBS(N,M) = YBS_TMP(N)
          ZBS(N,M) = ZBS_TMP(N)           
        ENDDO
        XBS(NCSP,M) = XBS(1,M) ! Close the trailing edge.
        YBS(NCSP,M) = YBS(1,M)
        ZBS(NCSP,M) = ZBS(1,M)
      ENDDO

C-----------------------------------------------------------------------
C     compute the arclength on the blade, ds
C-----------------------------------------------------------------------
      DO 80 M=1,MR
         I1=0
         DO 70 N=1,NC
            SUM=ZERO
            DO 60 I=1,10
               I1=I1+1
               DU=HALF*(XBS(I1+1,M)+XBS(I1+1,M+1)-XBS(I1,M)-XBS(I1,M+1))
               DV=HALF*(YBS(I1+1,M)+YBS(I1+1,M+1)-YBS(I1,M)-YBS(I1,M+1))
               DW=HALF*(ZBS(I1+1,M)+ZBS(I1+1,M+1)-ZBS(I1,M)-ZBS(I1,M+1))
               SUM=SUM+SQRT(DU*DU+DV*DV+DW*DW)
 60         CONTINUE
            DS(N,M)=SUM
 70      CONTINUE
 80   CONTINUE
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine DSCOMP>>>>>>>>>>>>>>>>>>>>>>>>>
      END
















