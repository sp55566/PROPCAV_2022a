       SUBROUTINE MODGEO(IFORMAT,ICAM,ICON,NX,NC)
C-----------------------------------------------------------------------
C     This subroutine modifies the suction side geometry aft of the 
C     midchord to give zero trailing edge thickness.
C
C     Date                 Comments
C     ------------------   ----------------------------------------
C     JY111800             Subroutine created.
C-----------------------------------------------------------------------
      INCLUDE 'PARAM.INC'
       COMMON /MODGEO1/ XC(MXINP),XTMP(16),
     %                  YP(MXINP,NXMAX),YS(MXINP,NXMAX)
       COMMON/GINP/XR(NXMAX),XRSQ(NXMAX),XPI(NXMAX),XRAKE(NXMAX),
     *             XSKEW(NXMAX),XCHD(NXMAX),XCI(NXMAX),XTI(NXMAX),
     *             XVA(NXMAX),XVR(NXMAX),XVT(NXMAX),XUA(NXMAX),
     *             XUAU(NXMAX),XUT(NXMAX),XUTU(NXMAX)
       COMMON /INPGEO1/ XCS(NXMAX,15),XTS(NXMAX,16)
       DIMENSION YTMP(MXINP),CUBTMP(4*(MXINP-1))
       DIMENSION CAM(MXINP),THK(MXINP),XCB(MXINP),XCT(MXINP),
     *           XXX(MXINP),YYY(MXINP),YYY1(MXINP),CUBTMP1(4*(MXINP-1)),
     *           XCOLD(MXINP),XKEEP(MXINP,NXMAX),YKEEP(MXINP,NXMAX)

       DO I=1,NX
          DO N=1,NC
             XKEEP(N,I)=XC(N)*XCHD(I)*2.
             YKEEP(N,I)=YS(N,I)*XCHD(I)*2.+XR(I)
          END DO
       END DO

       PI=DATAN(1D0)*4

C.....Only modify the last 33% of the suction side geometry
       NKEEP=MIN(NC*2/3,NC-5)
       NKEEP1=NKEEP+1

C-----------------------------------------------------------------------
C      Plot original blade section geometry.
C-----------------------------------------------------------------------
       OPEN(25,FILE='bld1.plt',STATUS='UNKNOWN')
       WRITE(25,*) 'VARIABLES="c/R","y/R"'
 100   FORMAT(1X,'ZONE T="r/R=',F8.3,'"')
       IF(IFORMAT.NE.2) THEN
          DO I=1,NX
             WRITE(25,100) XR(I)
             WRITE(25,*) 0.,XR(I)
             DO N=1,NC
                WRITE(25,*) XC(N)*XCHD(I)*2.,YP(N,I)*XCHD(I)*2.+XR(I)
             END DO
             DO N=NC,1,-1
                WRITE(25,*) XC(N)*XCHD(I)*2.,YS(N,I)*XCHD(I)*2.+XR(I)
             END DO
             WRITE(25,*) 0.,XR(I)
          END DO
          CLOSE(25)
       END IF

       OPEN(25,FILE='bld-mod.plt',STATUS='UNKNOWN')
       WRITE(25,*) 'VARIABLES="c/R","y/R"'

       DO I=1,NX

C-----------------------------------------------------------------------
C      Modify suction-side geometry if IFORMAT.NE.2
C-----------------------------------------------------------------------
          IF(IFORMAT.NE.2) THEN

C........modify suction-side geometry
             DO N=1,NKEEP
                XXX(N)=XC(N)
                YYY(N)=YS(N,I)
             END DO
             XXX(NKEEP1)=XC(NC)
             YYY(NKEEP1)=YP(NC,I)

             CALL UGLYDK(NKEEP1,1,1,XXX,YYY,0.0,0.0,CUBTMP1)
             CALL EVALDK(NKEEP1,NC,XXX,XC,YYY1,CUBTMP1)
             DO N=1,NC-1
                YS(N,I)=YYY1(N)
             END DO
             YS(NC,I)=YP(NC,I)

C...........plot new blade section geometry
             WRITE(25,100) XR(I)
             WRITE(25,*) 0.,XR(I)
             DO N=1,NC
                WRITE(25,*) XC(N)*XCHD(I)*2.,YP(N,I)*XCHD(I)*2.+XR(I)
             END DO
             DO N=NC,1,-1
                WRITE(25,*) XC(N)*XCHD(I)*2.,YS(N,I)*XCHD(I)*2.+XR(I)
             END DO
             WRITE(25,*) 0.,XR(I)

C........calculate the new chord, skew, rake, and pitch due to the
C........change in the suction side
             XCHDNEW=SQRT(XCHD(I)**2.+(XTS(I,NC)*XCHD(I))**2./4.)
             PITCH=ATAN(XPI(I)/PI/XR(I))

             IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) PITCH=0.

             DPITCH=ACOS(XCHD(I)/XCHDNEW)
             PITCHNEW=DPITCH+PITCH
             XRAKE(I)=XRAKE(I)+0.5*(XCHDNEW*SIN(PITCHNEW)-
     *            XCHD(I)*SIN(PITCH))
             SKEWOLD=XSKEW(I)/180.*PI
             SKEWNEW=SKEWOLD-(XCHD(I)*COS(PITCH)-XCHDNEW*COS(PITCHNEW))
     *            /XR(I)
             XSKEW(I)=SKEWNEW/PI*180.
             XPI(I)=PI*XR(I)*TAN(PITCHNEW)

             DO N=1,NC
                XCOLD1=XC(N)
                XCB(N)=XCOLD1*COS(DPITCH)-YP(N,I)*SIN(DPITCH)
                XCT(N)=XCOLD1*COS(DPITCH)-YS(N,I)*SIN(DPITCH)
                YP(N,I)=XCOLD1*SIN(DPITCH)+YP(N,I)*COS(DPITCH)
                YS(N,I)=XCOLD1*SIN(DPITCH)+YS(N,I)*COS(DPITCH)
             END DO

             DCHD=XCHD(I)/XCHDNEW
             DO N=1,NC
                XCB(N)=XCB(N)*DCHD
                XCT(N)=XCT(N)*DCHD
                YP(N,I)=YP(N,I)*DCHD
                YS(N,I)=YS(N,I)*DCHD
             END DO
             XCHD(I)=XCHDNEW

             CALL UGLYDK(NC,1,1,XCB,YP(1,I),0.0,0.0,CUBTMP)
             CALL EVALDK(NC,NC-1,XCB,XTMP,YTMP,CUBTMP)
             DO N=1,NC-1
                YP(N,I)=YTMP(N)
             END DO
             CALL UGLYDK(NC,1,1,XCT,YS(1,I),0.0,0.0,CUBTMP)
             CALL EVALDK(NC,NC-1,XCT,XTMP,YTMP,CUBTMP)
             DO N=1,NC-1
                YS(N,I)=YTMP(N)
             END DO
             YP(NC,I)=0.0
             YS(NC,I)=0.0

C-----------------------------------------------------------------------
C     Modify suction-side geometry if IFORMAT.EQ.2.  
C
C     Note:  1) This option assumes the chord is defined as the distance
C               between the L.E. and the T.E. of the PRESSURE side.
C            2) Special treament is used when the T.E. of the pressure 
C               side extendS below the chord-line (like the one given 
C               in pg 51 of Dr. Olofsson's thesis.  In this case, the 
C               chord will be redifined (so will the associated 
C               variables) so that the chord is the distance between 
C               the L.E. and the point where the pressure side 
C               intersects the origiinal chordline.
C-----------------------------------------------------------------------
          ELSE

             IFLAG=0
             IF(YP(NC,I).LT.0.0) THEN
                IFLAG=1

C..............redefine the new chord length
                DO N=NC/2,NC
                   IF(YP(N,I).LT.0.0.AND.YP(N-1,I).GT.0.0) THEN
                      N2=N
                      N1=N-1
                   ELSE IF(YP(N,I).EQ.0) THEN
                      N2=N
                      N1=N
                   END IF
                END DO
                IF(N2.EQ.N1) THEN
                   XCNEW=XC(N1)
                ELSE
                   XCNEW=XC(N1)+(XC(N2)-XC(N1))/(YP(N2,I)-YP(N1,I))*
     *                  (-YP(N1,I))
                END IF
                NCOLD=NC
                DO N=1,NC
                   XCOLD(N)=XC(N)
                END DO
                NC=N2
                DO N=1,NC-1
                   XC(N)=XC(N)/XCNEW
                   YP(N,I)=YP(N,I)/XCNEW
                   YS(N,I)=YS(N,I)/XCNEW
                END DO
                XC(NC)=1.0
                YP(NC,I)=0.0
                YS(NC,I)=0.0
                DCHD=XCHD(I)*(1.-XCNEW)
                XCHD(I)=XCNEW*XCHD(I)

C..............redefine new rake and skew
                PITCH=ATAN(XPI(I)/PI/XR(I))
                SKEWOLD=XSKEW(I)/180.*XR(I)/2.*PI/COS(PITCH)
                SKEWNEW=SKEWOLD-DCHD/2.
                XSKEW(I)=SKEWNEW*2./XR(I)*COS(PITCH)*180./PI
                XRAKE(I)=XRAKE(I)+SKEWNEW*SIN(PITCH)

             END IF

C...........modify suction-side geometry
             DO N=1,NKEEP
                XXX(N)=XC(N)
                YYY(N)=YS(N,I)
             END DO
             XXX(NKEEP1)=XC(NC)
             YYY(NKEEP1)=YP(NC,I)
             
             CALL UGLYDK(NKEEP1,1,1,XXX,YYY,0.0,0.0,CUBTMP1)
             CALL EVALDK(NKEEP1,NC,XXX,XC,YYY1,CUBTMP1)
             DO N=1,NC
                YS(N,I)=YYY1(N)
             END DO 

C...........plot original blade section geometry
             WRITE(15,100) XR(I)
             WRITE(15,*) 0.,XR(I)
             DO N=1,NC
                WRITE(15,*) XC(N)*XCHD(I)*2.,
     *               YP(N,I)*XCHD(I)*2.+XR(I)
             END DO
             DO N=NCOLD,1,-1
                WRITE(15,*) XKEEP(N,I),YKEEP(N,I)
             END DO
             WRITE(15,*) 0.,XR(I)

C...........plot new blade section geometry
             WRITE(25,100) XR(I)
             WRITE(25,*) 0.,XR(I)
             DO N=1,NC
                WRITE(25,*) XC(N)*XCHD(I)*2.,
     *               YP(N,I)*XCHD(I)*2.+XR(I)
             END DO
             DO N=NC,1,-1
                WRITE(25,*) XC(N)*XCHD(I)*2.,
     *               YS(N,I)*XCHD(I)*2.+XR(I)
             END DO
             WRITE(25,*) 0.,XR(I)
                
          END IF

C-----------------------------------------------------------------------
C     Calculate new XTS, XCS, XCI, and XTI.
C-----------------------------------------------------------------------
C........calculating new camber and thickness distribution
          DO N=1,NC
             CAM(N)=(YP(N,I)+YS(N,I))/2.
             THK(N)=YS(N,I)-YP(N,I)
          END DO

C........calculating max camber and max thickness.  This information
C........is not used in PROPCAV if ICAM=99 and ITHK=99.
          FC=0.
          TC=0.
          DO N=1,NC
             FC=AMAX1(FC,CAM(N))
             TC=AMAX1(TC,THK(N))
          END DO
          IF(I.LT.NX) THEN
             TD=TC*XCHD(I)
          ELSE
             TD=0.0
          END IF
          XCI(I)=FC
          XTI(I)=TD

C........extrapolate camber and thickness for the X/C coordinates 
C........required in PROPCAV.
          IF(IFORMAT.EQ.2) THEN
             CALL UGLYDK(NC,1,1,XC,CAM,0.0,0.0,CUBTMP1)
             CALL EVALDK(NC,15,XC,XTMP,YTMP,CUBTMP1)
             DO N=1,15
                XCS(I,N)=YTMP(N)
             END DO
             CALL UGLYDK(NC,1,1,XC,THK,0.0,0.0,CUBTMP1)
             CALL EVALDK(NC,16,XC,XTMP,YTMP,CUBTMP1)
             DO N=1,16
                XTS(I,N)=YTMP(N)
             END DO 
          ELSE
             DO N=1,16
                IF(N.LE.15) XCS(I,N)=CAM(N)
                XTS(I,N)=THK(N)
             END DO
          END IF

C........Restore orginal x/C location for next strip.             
          IF(IFLAG.EQ.1.AND.IFORMAT.EQ.2) THEN
             NC=NCOLD
             DO N=1,NC
                XC(N)=XCOLD(N)
             END DO
          END IF

       END DO
       
       IF(IFORMAT.NE.2.AND.ICAM.NE.99) ICAM=99

       IF(IFORMAT.EQ.2) CLOSE(15)
       CLOSE(25)
C-----------------------------------------------------------------------
C     Write new geometry to bld-mod.geo
C-----------------------------------------------------------------------
       OPEN(20,FILE='bld-mod.geo',STATUS='UNKNOWN')
 200   FORMAT(15(2X,F8.4))
       WRITE(20,200) (XR(I),I=1,NX)
       WRITE(20,200) (XPI(I),I=1,NX)
       WRITE(20,200) (XRAKE(I),I=1,NX)
       WRITE(20,200) (XSKEW(I),I=1,NX)
       WRITE(20,200) (XCHD(I),I=1,NX)
       WRITE(20,200) (XCI(I),I=1,NX)
       WRITE(20,200) (XTI(I),I=1,NX)
       DO N=1,15
          WRITE(20,200) (XCS(I,N),I=1,NX)
       END DO
       DO N=1,16
          WRITE(20,200) (XTS(I,N),I=1,NX)
       END DO
       CLOSE(20)

       RETURN
       END
