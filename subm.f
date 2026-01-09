       SUBROUTINE SUBM
************************************************************************
*      This subroutine calculates submergence of key blade panels at   *
*      all time steps.                                                 *
************************************************************************

       use m_PSXYZ
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       DIMENSION XN(4,3),XSN(5,3)
!      COMMON/PSXYZ/XPS(NSTEP,NPBZ,5),YPS(NSTEP,NPBZ,5),
!    *     ZPS(NSTEP,NPBZ,5),IPSN(NPBZ,NSTEP),NEL(NSTEP)


       if (.NOT. allocated(xps)) then
         allocate(XPS(NSTEP,NPBZ,5),YPS(NSTEP,NPBZ,5))
         allocate(ZPS(NSTEP,NPBZ,5),IPSN(NPBZ,NSTEP),NEL(NSTEP))
       end if

 1000  FORMAT(1x,'VARIABLES="X","Y","Z"')
 1005  FORMAT(1x,'ZONE T="Time=',F8.1,'", N=',I5,', E=',I5,
     *      ' F=FEPOINT, ET=QUADRILATERAL')
 1015  FORMAT(3(1X,F14.8))
 1020  FORMAT(4(1X,I6))

C-----------------------------------------------------------------------
C      IPS=0: neglect effects of partially submerged panels
C         =1: include effects of partially submerged panels
C-----------------------------------------------------------------------
       IPS=1

C-----------------------------------------------------------------------
C     Loop over all time steps
C-----------------------------------------------------------------------
       DO KK=1,NTPREV

          DT1=-DELTAT*FLOAT(KK-1)

C-----------------------------------------------------------------------
C     Count number of number of panels 
C-----------------------------------------------------------------------
          NEL(KK)=0

C-----------------------------------------------------------------------
C     Tag submerged panels
C     ISUBM=0    :  panel is above F.S.
C     ISUBM=1    :  panel is submerged
C     ISUBM=2    :  panel is partially submerged
C     
C     IPSN       :  no. of nodes on the partially submerged panel
C     NPS        :  total no. of partially submerged panels
C     XSN,YSN,ZSN:  global nodal points of partially submerged panel
C-----------------------------------------------------------------------
          NPS=0
 
          DO M=MR,1,-1

C-----------------------------------------------------------------------
C     Find submerged panels on the blade
C-----------------------------------------------------------------------
             DO N=1,NC

C..............Determine the nodal coordinate of panel
                CALL NODE(1,N,M,XN)

                I=INDEXB(N,M)

C..............Determine if the four sides of the panels are cut by 
C..............the free surface
                CALL FSCUT(XN,XSN,DT1,IS,ID,ICNT)

C..............tag the panels
                IF(IS.EQ.4) THEN
                   ISUBM(I,KK)=1
                   IF(IPS.EQ.1) NEL(KK)=NEL(KK)+1
                ELSE IF(ID.EQ.4) THEN
                   ISUBM(I,KK)=0
                ELSE
                   NPS=NPS+1
                   ISUBM(I,KK)=2
                 
                   IPSN(NPS,KK)=ICNT
                   DO NN=1,ICNT
                      XPS(KK,NPS,NN)=XSN(NN,1)
                      YPS(KK,NPS,NN)=XSN(NN,2)
                      ZPS(KK,NPS,NN)=XSN(NN,3)
                   END DO
                   
                   IF(IPS.EQ.1) THEN
                      NEL(KK)=NEL(KK)+1
                      IF(ICNT.GT.4) NEL(KK)=NEL(KK)+1
                   END IF
                END IF          

             END DO

C-----------------------------------------------------------------------
C     Calculate parameters for blade
C-----------------------------------------------------------------------
C             CALL SBLD(M,KK,IPS,NEL(KK))

C-----------------------------------------------------------------------
C     Find submerged panels on the wake
C-----------------------------------------------------------------------
             ICW(M,KK)=0
             MSW(M,KK)=0

             I1=INDEXB(1,M)
             INC=INDEXB(NC,M)
             IF(IPS.EQ.0) THEN
C                IF(ICB(M,1,KK).EQ.0.OR.ICB(M,2,KK).EQ.0) GO TO 100
                IF(ISUBM(I1,KK).NE.1.OR.ISUBM(INC,KK).NE.1) GO TO 100
             ELSE
                IF(ISUBM(I1,KK).EQ.0.OR.ISUBM(INC,KK).EQ.0) GO TO 100
             END IF

             IFLAG=0
             DO N=1,NWMIN

                I=IDXWAK(N,M)
                ISUWM(I,KK)=0

                IF(IFLAG.EQ.1) GO TO 100

C..............Determine the nodal coordinate of panel
                CALL NODE(3,N,M,XN)

                   
C..............Determine if the four sides of the panels are cut by 
C..............the free surface
                CALL FSCUT(XN,XSN,DT1,IS,ID,ICNT)

C..............tag the panels
                IF(IS.EQ.4) THEN
                   ISUWM(I,KK)=1
                   ICW(M,KK)=1
                   NEL(KK)=NEL(KK)+1
                ELSE IF(ID.EQ.4) THEN
                   ISUWM(I,KK)=0
                ELSE
                   NPS=NPS+1
                   ISUWM(I,KK)=2
                 
                   IPSN(NPS,KK)=ICNT
                   DO NN=1,ICNT
                      XPS(KK,NPS,NN)=XSN(NN,1)
                      YPS(KK,NPS,NN)=XSN(NN,2)
                      ZPS(KK,NPS,NN)=XSN(NN,3)
                   END DO

                   IF(IPS.EQ.1) THEN
                      NEL(KK)=NEL(KK)+1
                      IF(ICNT.GT.4) NEL(KK)=NEL(KK)+1
                   END IF
                END IF          

                IM1=IDXWAK(N-1,M)
                IF(IPS.EQ.1) THEN
                   IF(ISUWM(I,KK).EQ.0.AND.ISUWM(IM1,KK).GT.0) IFLAG=1
                ELSE
                   IF(ISUWM(I,KK).NE.1.AND.ISUWM(IM1,KK).EQ.1) IFLAG=1
                END IF

             END DO

C-----------------------------------------------------------------------
C     Calculate parameters for wake
C-----------------------------------------------------------------------
C             CALL SWAKE(M,KK)

 100         CONTINUE

          END DO

C          GO TO 3000

C-----------------------------------------------------------------------
C     Find submerged panels on the hub
C-----------------------------------------------------------------------
          DO N=1,NHBX
             DO M=1,MHBT

C..............Determine the nodal coordinate of panel
                CALL NODE(2,N,M,XN)

                I=INDEXH(N,M)

C..............Determine if the four sides of the panels are cut by 
C..............the free surface
                CALL FSCUT(XN,XSN,DT1,IS,ID,ICNT)

C..............tag the panels
                IF(IS.EQ.4) THEN
                   ISUBM(I,KK)=1
                   IF(IPS.EQ.1) NEL(KK)=NEL(KK)+1
                ELSE IF(ID.EQ.4) THEN
                   ISUBM(I,KK)=0
                ELSE
                   NPS=NPS+1
                   ISUBM(I,KK)=2
                   
                   IPSN(NPS,KK)=ICNT
                   DO NN=1,ICNT
                      XPS(KK,NPS,NN)=XSN(NN,1)
                      YPS(KK,NPS,NN)=XSN(NN,2)
                      ZPS(KK,NPS,NN)=XSN(NN,3)
                   END DO

                   NEL(KK)=NEL(KK)+1
                   IF(ICNT.GT.4) NEL(KK)=NEL(KK)+1
                END IF          

             END DO
          END DO

 3000     CONTINUE

       END DO

C.....Plot submerged panels
       CALL SUBMPLT2

       STOP

       RETURN
       END


       SUBROUTINE FSCUT(XN,XSN,DT1,IS,ID,ICNT)
       INCLUDE 'PUFCAV.INC'
       DIMENSION XX(2),YY(2),ZZ(2),XN(4,3),XSN(5,3)

       ICNT=0
       IS=0
       ID=0

C.....Loop over all the edges of the panel
       DO ISIDE=1,4

          N1=ISIDE
          IF(ISIDE.LT.4) THEN
             N2=ISIDE+1
          ELSE
             N2=1
          END IF

C........Rotate to global coordinate
          XX(1)=XN(N1,1)
          XX(2)=XN(N2,1)                   
          CALL ROTATE(XN(N1,2),XN(N1,3),YY(1),ZZ(1),DT1)
          CALL ROTATE(XN(N2,2),XN(N2,3),YY(2),ZZ(2),DT1)

C........Check if edge ISIDE is cut by the free surface
          IF(YY(1).LE.YY(2)) THEN
             IB=1
             IT=2
          ELSE
             IB=2
             IT=1
          END IF
          DEN=YY(IT)-YY(IB)
          IF(DEN.EQ.ZERO) THEN
             DUM=YFS-YY(IB)
             IF(DUM.GT.ZERO) CUT=MAX(ONE,DUM)
          ELSE
             CUT=(YFS-YY(IB))/DEN
          END IF
          
          IF(CUT.GE.ONE) THEN
             IS=IS+1
             ICNT=ICNT+1
             XSN(ICNT,1)=XX(1)
             XSN(ICNT,2)=YY(1)
             XSN(ICNT,3)=ZZ(1)
          ELSE IF(CUT.LE.ZERO) THEN
             ID=ID+1
          ELSE 
             IF(IB.EQ.1) THEN
                ICNT=ICNT+1
                XSN(ICNT,1)=XX(1)
                XSN(ICNT,2)=YY(1)
                XSN(ICNT,3)=ZZ(1)
                ICNT=ICNT+1
                XSN(ICNT,1)=XX(IB)+CUT*(XX(IT)-XX(IB))
                XSN(ICNT,2)=YFS
                XSN(ICNT,3)=ZZ(IB)+CUT*(ZZ(IT)-ZZ(IB))
             ELSE
                ICNT=ICNT+1
                XSN(ICNT,1)=XX(IB)+CUT*(XX(IT)-XX(IB))
                XSN(ICNT,2)=YFS
                XSN(ICNT,3)=ZZ(IB)+CUT*(ZZ(IT)-ZZ(IB))
             END IF
          END IF
       END DO

       RETURN
       END
