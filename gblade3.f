      SUBROUTINE GBLADE3
************************************************************************
*     GBLADE: Geometry of the special surface-piercing BLADE           *
*      --- Generate a blade surface geometry                           *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     JY083000 This subroutine added for a special type of surface-    *
*              piercing propeller for testing of the code.             *
*                                                                      *
************************************************************************

      use m_INPGEO2
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'INPGEO.INC'
    
      DIMENSION XX1(NBHZP),YY1(NBHZP),XICUB(4*NBHZP-4),
     *     XIM(NBHZP),ETAM(NBHZP),CHORDOLD(MBPZ),XBM(NBHZP,MBPZ),
     *     YBM(NBHZP,MBPZ),ZBM(NBHZP,MBPZ)

      ANGB=FLOAT(IANG)*DELTAT
      ANGB2=ANGB/2.
      DANG=DELTAT/FLOAT(INH)

      INH1=INH*2
      DANG1=DELTAT/FLOAT(INH1)

      NH=2*INH1+(IANG-2)*INH
      NC=NH*2

      WRITE(*,*) 'IANG,INH,NC,DANG',IANG,INH,NC,DANG

      NHP=NH+1
      NHM=NH-1
      NCP=NC+1
      MRP=MR+1
      NPANB=NC*MR
      
      IF(NC.GT.NBZ) THEN
         WRITE(*,*) 'Need to increase NBZ in PARAM.INC!'
         WRITE(*,*) 'NC,NBZ:',NC,NBZ
         STOP
      END IF

C-----------------------------------------------------------------------
C     Store old chord length
C-----------------------------------------------------------------------
      DO M=1,MRP
         CHORDOLD(M)=CHORD(M)
      END DO

C-----------------------------------------------------------------------
C     Calculate chord length: CHORD=C/D, CHORDZ=CHORD*COS(PITCH)
C
C     Calculate Y/R & Z/R as well as LE coordinates.
C-----------------------------------------------------------------------
      DO M=1,MRP
         TANP=PITCH(M)/PI/RZ(M)
         COSP=1.0/SQRT(1+TANP*TANP)
         SINP=TANP*COSP
         YB(NHP,M)=RZ(M)*COS(-ANGB2)
         ZB(NHP,M)=RZ(M)*SIN(-ANGB2)
         YB(1,M)=RZ(M)*COS(ANGB2)
         ZB(1,M)=RZ(M)*SIN(ANGB2)
         YB(NCP,M)=YB(1,M)
         ZB(NCP,M)=ZB(1,M)

         CHORDZ=SQRT((YB(NHP,M)-YB(1,M))**2.+
     *        (ZB(NHP,M)-ZB(1,M))**2.)
         CHORD(M)=CHORDZ/COSP/TWO

         XB(NHP,M)=RAKE(M)*TWO-CHORD(M)*SINP
         XB(1,M)=RAKE(M)*TWO+CHORD(M)*SINP
         XB(NCP,M)=XB(1,M)


         DO N=2,NH
            IF(N.LE.INH1+1) THEN
               YB(NH+N,M)=RZ(M)*COS(-ANGB2+DANG1*(N-1))
               ZB(NH+N,M)=RZ(M)*SIN(-ANGB2+DANG1*(N-1))
            ELSE IF(N.GE.NHP-INH1) THEN
               YB(NH+N,M)=RZ(M)*COS(ANGB2-DANG1*(NHP-N))
               ZB(NH+N,M)=RZ(M)*SIN(ANGB2-DANG1*(NHP-N))
            ELSE
               YB(NH+N,M)=RZ(M)*COS(-ANGB2+DANG*(N-INH1+INH-1))
               ZB(NH+N,M)=RZ(M)*SIN(-ANGB2+DANG*(N-INH1+INH-1))
            END IF
            YB(NHP-N+1,M)=YB(NH+N,M)
            ZB(NHP-N+1,M)=ZB(NH+N,M)
         END DO

         DO N=1,NHP
            YBM(N,M)=YB(NH+N,M)
            ZBM(N,M)=ZB(NH+N,M)
            IF(N.EQ.1.OR.N.EQ.NHP) XBM(N,M)=XB(NH+N,M)
         END DO
      END DO

C-----------------------------------------------------------------------
C     Temporary chordwise spacing for 2D section: 
C-----------------------------------------------------------------------
      DO N=1,NH
         SB(N)=FLOAT(N)/FLOAT(NH)
      END DO
      
C-----------------------------------------------------------------------
C     Generate blade geometry
C-----------------------------------------------------------------------
      DO M=1,MRP

C........New F/C ratio due to new chord length
         CAMBR(M)=CAMBR(M)*CHORD(M)/CHORDOLD(M)

C........Midchord centered local coordinate xi and eta 
C........(nondimensionalized by R)
C........CHORD=C/D,THK=T/D,CAMBR=F/C,YTC=Y/D,YCC=Y/C,SB=S/C,RAKE=XM/D
C
         IF(ITHK.EQ.99) THEN
            CALL THKINP(NH,XTS1(1,M),SB,YTC)
            GO TO 1000
         END IF

CT.......For this special case, the last line is t/C, not t/D (JY091000)
         THICK(M)=THICK(M)*CHORD(M)

         IF(THICK(M).EQ.0.0) THEN
            DO N=1,NH
               YTC(N)=ZERO
            END DO
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

         IF(ICAM.EQ.0) THEN
            CALL A8ML(NH,CAMBR(M),SB,YCC,TCC)
         ELSE IF(ICAM.EQ.1) THEN
            CALL CBPB(NH,CAMBR(M),SB,YCC,TCC)
         ELSE IF(ICAM.EQ.99) THEN
            CALL CAMINP(NH,XCS1(1,M),SB,YCC,TCC)
         END IF

C........local coordinate xi and eta
         XI(NHP)=-CHORD(M)
         ETA(NHP)=ZERO
         XIM(1)=XI(NHP)
         ETAM(1)=ZERO
         DO 20 N=1,NH
            DXT=YTC(N)*SIN(TCC(N))
            DYT=YTC(N)*COS(TCC(N))
            NBOT=NHP-N
            NTOP=NHP+N
            NMEAN=N+1
            XI(NTOP)=( (SB(N)-.5)*CHORD(M)-DXT )*TWO
            ETA(NTOP)=(YCC(N)*CHORD(M) +DYT)*TWO
            XI(NBOT)=( (SB(N)-.5)*CHORD(M)+DXT )*TWO
            ETA(NBOT)=(YCC(N)*CHORD(M) -DYT)*TWO
            XIM(NMEAN)=(SB(N)-.5)*CHORD(M)*TWO
            ETAM(NMEAN)=YCC(N)*CHORD(M)*TWO
20       CONTINUE

         TANP=PITCH(M)/PI/RZ(M)
         COSP=1.0/SQRT(1+TANP*TANP)
         SINP=TANP*COSP

C.......X-coord on pressure side
         DO N=NH,2,-1
            XX1(NHP-N+1)=XI(N)*COSP+ETA(N)*SINP
            YY1(NHP-N+1)=XI(N)*SINP-ETA(N)*COSP
         END DO
         XX1(1)=ATAN(ZB(NHP,M)/YB(NHP,M))*RZ(M)
         YY1(1)=XB(NHP,M)-RAKE(M)*TWO
         XX1(NHP)=ATAN(ZB(1,M)/YB(1,M))*RZ(M)
         YY1(NHP)=XB(1,M)-RAKE(M)*TWO

         CALL UGLYDK(NHP,1,1,XX1,YY1,0.0,0.0,XICUB)

         DO N=NH,2,-1
            THETA=ATAN(ZB(N,M)/YB(N,M))
            XX2=THETA*RZ(M)
            CALL EVALDKs(NHP,1,XX1,XX2,YY2,XICUB)
            XB(N,M)=YY2+RAKE(M)*TWO
         END DO

C.......X-coord on suction side
         DO N=NHP+1,NC
            XX1(N-NH)=XI(N)*COSP+ETA(N)*SINP
            YY1(N-NH)=XI(N)*SINP-ETA(N)*COSP
         END DO
         XX1(1)=ATAN(ZB(NHP,M)/YB(NHP,M))*RZ(M)
         YY1(1)=XB(NHP,M)-RAKE(M)*TWO
         XX1(NHP)=ATAN(ZB(1,M)/YB(1,M))*RZ(M)
         YY1(NHP)=XB(1,M)-RAKE(M)*TWO

         CALL UGLYDK(NHP,1,1,XX1,YY1,0.0,0.0,XICUB)

         DO N=NHP+1,NC
            THETA=ATAN(ZB(N,M)/YB(N,M))
            XX2=THETA*RZ(M)            
            CALL EVALDKs(NHP,1,XX1,XX2,YY2,XICUB)
            XB(N,M)=YY2+RAKE(M)*TWO
         END DO         

C.......X-coord on mean camber surface
         DO N=1,NHP
            XX1(N)=XIM(N)*COSP+ETAM(N)*SINP
            YY1(N)=XIM(N)*SINP-ETAM(N)*COSP
         END DO

         CALL UGLYDK(NHP,1,1,XX1,YY1,0.0,0.0,XICUB)

         DO N=2,NH
            THETA=ATAN(ZBM(N,M)/YBM(N,M))
            XX2=THETA*RZ(M)            
            CALL EVALDKs(NHP,1,XX1,XX2,YY2,XICUB)
            XBM(N,M)=YY2+RAKE(M)*TWO
         END DO         

      END DO

c      WRITE(999,*) 'ZONE T="1", I=',NCP,', J=',MRP
c      DO M=1,MRP
c         DO N=1,NCP
c            WRITE(999,*) XB(N,M),YB(N,M),ZB(N,M)
c         END DO
c      END DO
c      WRITE(999,*) 'ZONE T="2", I=',NHP,', J=',MRP
c      DO M=1,MRP
c         DO N=1,NHP
c            WRITE(999,*) XBM(N,M),YBM(N,M),ZBM(N,M)
c         END DO
c      END DO

c      DO M=1,MRP
c         WRITE(999,*) 'ZONE T="M=',M,'"'
c         DO N=1,NCP
c            WRITE(999,*) XB(N,M),ZB(N,M)
c         END DO
c      END DO

c      STOP

C-----------------------------------------------------------------------
C     Compute the normal vectors of the camber surface
C-----------------------------------------------------------------------
      DO N=1,NH
         DO M=1,MR
            XG(1,1,1)=XBM(N,M)
            XG(1,1,2)=YBM(N,M)
            XG(1,1,3)=ZBM(N,M)
            XG(1,2,1)=XBM(N,M+1)
            XG(1,2,2)=YBM(N,M+1)
            XG(1,2,3)=ZBM(N,M+1)
            XG(1,3,1)=XBM(N+1,M+1)
            XG(1,3,2)=YBM(N+1,M+1)
            XG(1,3,3)=ZBM(N+1,M+1)
            XG(1,4,1)=XBM(N+1,M)
            XG(1,4,2)=YBM(N+1,M)
            XG(1,4,3)=ZBM(N+1,M)

            CALL GEOM3D(1,XG,CHRLEPS,IER)
            IF(IER.EQ.0) THEN
               WRITE(*,'(A)') ' UNACCEPTABLE PANELS ON CAMBER SURFACE'
               STOP
            END IF
            XCON(N,M)=VEL(1,1)
            YCON(N,M)=VEL(1,2)
            ZCON(N,M)=VEL(1,3)
         END DO
      END DO

      RETURN
      END










