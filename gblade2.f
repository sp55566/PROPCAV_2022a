      SUBROUTINE GBLADE2
************************************************************************
*     GBLADE: Geometry of the wing BLADE                               *
*      --- Generate a blade surface geometry                           *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     CM010598 This subroutine added for the hydrofoil geometry case   *
*              by H.S. Lee                                             *
*     JY010798 Added a parabolic camber and thickness option for       *
*           blade geometry definition.  Here are the new               *
*           additions:                                                 *
*           ITHK = 10   Parabolic thickness form                       *
*           ICAM = 0    NACA A=0.8 mean line                           *
*           ICAM = 1    Parabolic camber line                          *
*     JY112198 Added a new camber and a new thickness option:          *
*                     ICAM=99  User input camber distribution          *
*                     ITHK=99  User input thickness distribution       *
*                                                                      *
************************************************************************

      use m_INPGEO2
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'INPGEO.INC'

      DIMENSION XIM(NBHZP),ETAM(NBHZP)
      DIMENSION XIOLD(NBPZ,MBPZ),ETAOLD(NBPZ,MBPZ),
     *     XIMOLD(NBHZP,MBPZ),ETAMOLD(NBHZP,MBPZ)

C-----------------------------------------------------------------------
C     Set up chordwise spacing
C-----------------------------------------------------------------------
      IF(ICSPAC.EQ.2.OR.ICSPAC.EQ.3) THEN
         RLET=RADLE(XTI(1),XCHD(1),ITHK)
      END IF
      CALL SPACE(NC,ICSPAC,RLET,SB)

C-----------------------------------------------------------------------
C     Generate blade geometry
C-----------------------------------------------------------------------
      DO 50 M=1,MRP
          COST=1
C         COST=COS(PITCH(M))
C         SINT=SIN(PITCH(M))
          SINT=0

C........Midchord centered local coordinate xi and eta (nondimensionaliz
C........  by R
C........CHORD=C/D,THK=T/D,CAMBR=F/C,YTC=Y/D,YCC=Y/C,SB=S/C,RAKE=XM/D
C
C-----------------------------------------------------------------------
C     ITHK=99, user input thickness distribution.               JY112798
C-----------------------------------------------------------------------
         IF(ITHK.EQ.99) THEN
            CALL THKINP(NH,XTS1(1,M),SB,YTC)
            GO TO 1000
         END IF

         IF(THICK(M).EQ.0.0) THEN
            DO 10 N=1,NH
               YTC(N)=ZERO
10          CONTINUE 
         ELSE
            IF(ITHK.EQ.1) THEN
               CALL NACA66(NH,THICK(M),RLE,SB,YTC)
            ELSE IF(ITHK.EQ.2) THEN
               CALL RAE(NH,THICK(M),SB,YTC)
            ELSE IF(ITHK.EQ.3) THEN
               CALL NACA65A(NH,THICK(M),RLE,SB,YTC)
            ELSE IF(ITHK.EQ.4) THEN 
               CALL NACA64A(NH,THICK(M),RLE,SB,YTC)
            ELSE IF(ITHK.EQ.5) THEN
               CALL NACA00(NH,THICK(M),RLE,SB,YTC)
            ELSE IF(ITHK.EQ.6) THEN
               CALL ELIPSE(NH,THICK(M),RLE,SB,YTC)
            ELSE IF(ITHK.EQ.9) THEN
               CALL NACA16(NH,THICK(M),RLE,SB,YTC)
            ELSE IF(ITHK.EQ.10) THEN
               CALL TKPB(NH,THICK(M),RLE,SB,YTC)
            END IF

         END IF

 1000    CONTINUE

C-----------------------------------------------------------------------
C        New if statement allows user to select camber distribution
C        using the variable ICAM.  The available cambers are:
C        
C        ICAM = 0      A=0.8 Mean Line Camber
C        ICAM = 1      Parabolic Camber Distribution           JY010898
C-----------------------------------------------------------------------
         IF(ICAM.EQ.0) THEN
            CALL A8ML(NH,CAMBR(M),SB,YCC,TCC)
         ELSE IF(ICAM.EQ.1) THEN
            CALL CBPB(NH,CAMBR(M),SB,YCC,TCC)
C-----------------------------------------------------------------------
C       ICAM=99, user input thickness distribution.            JY112798
C-----------------------------------------------------------------------
         ELSE IF(ICAM.EQ.99) THEN
            CALL CAMINP(NH,XCS1(1,M),SB,YCC,TCC)
         END IF

        XI(NHP) = -CHORD(M) * COST
        ETA(NHP) = CHORD(M) * SINT
          XIM(1)=XI(NHP)
          ETAM(1)=ZERO

          DO 20 N=1,NH

            NBOT=NHP-N
            NTOP=NHP+N
            NMEAN=N+1

          XXX = (SB(N) - 0.5) * CHORD(M) * TWO
          YYY = YCC(N) * CHORD(M) * TWO

               
            IF(ITHK .EQ. 99 .AND. ICAM .EQ. 99 
     %           .AND. IFORMAT .EQ. 2) THEN
               DXT = 0.0
               DYT = YTC(N)
            ELSE
               DXT = YTC(N) * SIN(TCC(N)) * TWO
               DYT = YTC(N) * COS(TCC(N)) * TWO
            ENDIF

          XU = XXX - DXT
          XL = XXX + DXT

          YU = YYY + DYT
          YL = YYY - DYT
         
          XI(NTOP) = XU * COST + YU * SINT
          ETA(NTOP) = -XU * SINT + YU * COST

          XI(NBOT) = XL * COST + YL * SINT
          ETA(NBOT) = -XL * SINT + YL * COST

            XIM(NMEAN)=(SB(N)-.5)*CHORD(M)*TWO
            ETAM(NMEAN)=YCC(N)*CHORD(M)*TWO
20       CONTINUE

         IF(ISC.EQ.1) THEN
C..........Store coordinates of original blades.
            DO N=1,NCP
               XIOLD(N,M)=XI(N)
               ETAOLD(N,M)=ETA(N)
               IF(N.LE.NHP) THEN
                  XIMOLD(N,M)=XIM(N)
                  ETAMOLD(N,M)=ETAM(N)
               END IF
            END DO
         END IF

         DO 31 N=1,NCP
            XB(N,M)=RAKE(M) * TWO + XI(N)  
            YB(N,M)=RZ(M)
            ZB(N,M)=ETA(N) + SKEW(M) * TWO
 31      CONTINUE
CSH------------------------Shreenaath Natarajan------04/02/2003----------
CSH-------Changes for FLAP
         IF((M.EQ.1).AND.(IFLAP.EQ.1)) CALL XFLAPPT
CSH
C.......Plot blade sectsions (JY071201)
         IF(M.EQ.1) THEN
            OPEN(55,FILE='bldsec.plt',STATUS='UNKNOWN')
            WRITE(55,*) 'VARIABLES="c/R","y/R"'
         END IF
         WRITE(55,*) 'ZONE T="M=',M,'"'
         DO N=1,NCP
            WRITE(55,*) XI(N)+CHORD(M),ETA(N)+RZ(M)
C...... Added XIV, ETAV for 2-D strip node coordinates calculation 
C                    By  Hong Sun            11/04/05
            XIV(N,M) = XI(N)
            ETAV(N,M) = ETA(N)
         END DO
         IF(M.EQ.MRP) CLOSE(55)

50    CONTINUE
CSH------------------------Shreenaath Natarajan------04/02/2003----------
CSH-------Changes for FLAP
       IF(IFLAP.EQ.1)   CALL GFLAP
CSH
C-----------------------------------------------------------------------
C     Create geometry for the separated region.                 JY081601
C-----------------------------------------------------------------------
      IF(ISC.EQ.1) CALL SRGEO(XIOLD,ETAOLD,XIMOLD,ETAMOLD)

      RETURN
      END










