       SUBROUTINE M841BFLOW
************************************************************************
*     M841BFLOW: INFLOW velocities using distribution shown on page 56 *
*                of Dr. Olofsson's dissertation.                       *
*      --- Compute the source strength of each blade at current time   *
*            step                                                      *
*      --- Compute the total inflow velocities at the blade panels and *
*          wake panels                                                 *
*                                                                      *
* Date        Revision                                                 *
* ----        --------                                                 *
* JY091000    Subroutine Created.                                      *
*                                                                      *
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'
       COMMON/M841BVEL/ XXCUB(40),XXY(11)
       DIMENSION XXV(11)

       DATA XXV/0.000E+00,6.000E-01,7.300E-01,8.800E-01,9.400E-01, 
     *      9.750E-01,9.900E-01,9.950E-01,1.000E+00,1.000E+00,
     *      1.000E+00/
       DATA XXY/3.300E-01,3.500E-01,3.650E-01,4.150E-01,4.500E-01, 
     *      4.950E-01,5.300E-01,5.700E-01,6.150E-01,6.550E-01,
     *      1.000E+00/

       IF(NTSTEP.EQ.1) THEN
          XXY0=XXY(1)
          DO N=1,11
             XXY(N)=XXY(N)-XXY0
          END DO
C          WRITE(999,*) 'ZONE T="ORGINAL"'
C          WRITE(999,'(2(1X,F14.6))') (XXV(N),XXY(N),N=1,11)
       
          CALL UGLYDK(11,1,2,XXY,XXV,0.0,0.0,XXCUB)

C          WRITE(999,*) 'ZONE T="SPLINE"'
       END IF

       DO KK=1,NBLADE
          DTBLA=DELK*FLOAT(KK-1)

          DO J=1,NPANEL
             RCP=SQRT(XCTP(J,2,1)**2+XCTP(J,3,1)**2)
             THP=ATAN2(XCTP(J,3,1),XCTP(J,2,1))+DTBLA
             T=THP-TSTEP
             CT1=COS(T)
             YG1=RCP*CT1
             YG1=YFS-YG1             

C...........Non-homogenous (due to the flat plate at F.S.) w/o 
C...........propeller
             IF(YG1.LE.ZERO) THEN
                VOX=ZERO
             ELSE IF(YG1.GT.XXY(11)) THEN
                VOX=ONE
             ELSE
                XXY1=YG1
                CALL EVALDKs(11,1,XXY,XXY1,XXV1,XXCUB)
                VOX=XXV1
             END IF             
             
C             WRITE(999,'(2(1X,F14.6))') VOX,YG1

C...........Add homogenous, rotational flow field
             COCP=XCTP(J,2,KK)/RCP
             SICP=XCTP(J,3,KK)/RCP
             WROVS=PI*RCP/ADVCO
             VOY=-WROVS*SICP
             VOZ=WROVS*COCP

             VELY=VEL(J,2)*COS(DTBLA)-VEL(J,3)*SIN(DTBLA) 
             VELZ=VEL(J,2)*SIN(DTBLA)+VEL(J,3)*COS(DTBLA) 
             BUG=VOX*VEL(J,1)+VOY*VELY+VOZ*VELZ

             IF(KK.EQ.1) THEN
                KK0=1
                DPDN(J)=-BUG
                VOX1(J)=VOX
                VOY1(J)=VOY
                VOZ1(J)=VOZ
             ELSE
                KK0=NBLADE+2-KK
             END IF

             STRGTH(J,KK0)=-BUG
             
          END DO

          IF(KK.EQ.1) THEN
C-----------------------------------------------------------------------
C     On-coming velocity in local coordinate system
C-----------------------------------------------------------------------
             DO N=1,NC
                DO M=1,MR
                   L=INDEXB(N,M)
                   VXIB(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
     *                  +VOZ1(L)*DIR(L,1,3)
                   VETAB(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
     *                  +VOZ1(L)*DIR(L,2,3)
                   VINFSB(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
                END DO
             END DO

             IF(IHUB.NE.0) THEN
                DO N=1,NHBX
                   DO M=1,MHBT
                      L=INDEXH(N,M)
                      VXIH(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
     *                     +VOZ1(L)*DIR(L,1,3)
                      VETAH(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
     *                     +VOZ1(L)*DIR(L,2,3)
                      VINFSH(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
                   END DO
                END DO
             END IF

          END IF

       END DO

       DO J=1,NPWAKS
          RCP=SQRT(XCPW(J,2,1)**2+XCPW(J,3,1)**2)
          THP=ATAN2(XCPW(J,3,1),XCPW(J,2,1))
          T=THP-TSTEP
          CT1=COS(T)
          YG1=RCP*CT1
          YG1=YFS-YG1             

C........Non-homogenous (due to the flat plate at F.S.) w/o 
C........propeller
          IF(YG1.LE.ZERO) THEN
             VOXW(J)=ZERO
          ELSE IF(YG1.GT.XXY(11)) THEN
             VOXW(J)=ONE
          ELSE
             XXY1=YG1
             CALL EVALDKs(11,1,XXY,XXY1,XXV1,XXCUB)
             VOXW(J)=XXV1
          END IF             
             
C          WRITE(999,'(2(1X,F14.6))') VOXW(J),YG1

C........Add homogenous, rotational flow field
          COCP=XCPW(J,2,1)/RCP
          SICP=XCPW(J,3,1)/RCP
          WROVS=PI*RCP/ADVCO
          VOYW(J)=-WROVS*SICP
          VOZW(J)=WROVS*COCP

       END DO

C       IF(NTSTEP.EQ.NTPREV) STOP
                
       RETURN
       END
