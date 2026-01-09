       SUBROUTINE SUBMPLT(NEL)

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'
       DIMENSION YB1(NBPZ,MBPZ),ZB1(NBPZ,MBPZ),YH1(NHPZ,MHPZ),
     *      ZH1(NHPZ,MHPZ),NEL(NSTEP)

C-----------------------------------------------------------------------
C     Plot submerged panels
C-----------------------------------------------------------------------
 1000  FORMAT(1x,'VARIABLES="X","Y","Z"')
 1005  FORMAT(1x,'ZONE T="Time=',F8.1,'", N=',I5,', E=',I5,
     *      ' F=FEPOINT, ET=QUADRILATERAL')
 1015  FORMAT(3(1X,F14.8))
 1020  FORMAT(4(1X,I6))
       
       OPEN(5,FILE='surface.plt',STATUS='UNKNOWN')
       WRITE(5,1000) 

       DO KK=1,NTPREV
          IT=(KK-1)*NDLTAT

          NEL1=0
          DO K1=1,NBLADE
             
             KK1=NTPOS1(K1,IT)
             NEL1=NEL1+NEL(KK1)
          END DO

          IF(NEL1.NE.0) THEN
             WRITE(5,1005) TT(KK),(NEL1+1)*4,NEL1+1
             
             WRITE(5,1015) -1.,YFS,-1.2
             WRITE(5,1015)  1.,YFS,-1.2
             WRITE(5,1015)  1.,YFS, 1.2             
             WRITE(5,1015) -1.,YFS, 1.2             
          END IF

          DO K1=1,NBLADE
              KK1=NTPOS1(K1,IT)
              DT1=-DELTAT*(KK1-1)

              DO M=1,MR+1
                 DO N=1,NC+1
                    YB1(N,M)=YB(N,M)*COS(DT1)-ZB(N,M)*SIN(DT1)
                    ZB1(N,M)=YB(N,M)*SIN(DT1)+ZB(N,M)*COS(DT1)
                 END DO

                 DO N=1,NWMIN+1
                    YW1(N,M)=YW(N,M)*COS(DT1)-ZW(N,M)*SIN(DT1)
                    ZW1(N,M)=YW(N,M)*SIN(DT1)+ZW(N,M)*COS(DT1)
                 END DO
              END DO

              DO M=1,MHBT+1
                 DO N=1,NHBX+1
                    YH1(N,M)=YH(N,M)*COS(DT1)-ZH(N,M)*SIN(DT1)
                    ZH1(N,M)=YH(N,M)*SIN(DT1)+ZH(N,M)*COS(DT1)
                 END DO
              END DO

C............Plot submerged panels on the blade.
              DO M=1,MR
                 DO IDR=1,2
                    IF(ICB(M,IDR,KK1).EQ.1) THEN
                       IF(IDR.EQ.1) THEN
                          N1=IC(1,M,KK1)
                          N2=IC(2,M,KK1)
                          KI=1
                       ELSE
                          N1=IW(2,M,KK1)
                          N2=IW(1,M,KK1)
                          KI=-1
                       END IF

                       DO N=N1,N2,KI
                         
                          IF(IDR.EQ.1) THEN
                             NN1=N
                             NN2=N+1
                          ELSE IF(IDR.EQ.2) THEN
                             NN2=N
                             NN1=N+1
                          END IF
                          
                          WRITE(5,1015) XB(NN1,M),YB1(NN1,M),
     *                         ZB1(NN1,M)
                          WRITE(5,1015) XB(NN1,M+1),YB1(NN1,M+1),
     *                         ZB1(NN1,M+1)
                          WRITE(5,1015) XB(NN2,M+1),YB1(NN2,M+1),
     *                         ZB1(NN2,M+1)
                          WRITE(5,1015) XB(NN2,M),YB1(NN2,M),
     *                         ZB1(NN2,M)
                          
                       END DO
                    END IF
                 END DO
              END DO

C............Plot submerged panels on the wake.
              DO M=1,MR
                 IF(ICW(M,KK1).EQ.1) THEN
                    N1=1
                    N2=MSW(M,KK1)           
                    DO N=N1,N2
                       NN1=N
                       NN2=N+1
                       
                       WRITE(5,1015) XW(NN1,M),YW1(NN1,M),
     *                      ZW1(NN1,M)
                       WRITE(5,1015) XW(NN1,M+1),YW1(NN1,M+1),
     *                      ZW1(NN1,M+1)
                       WRITE(5,1015) XW(NN2,M+1),YW1(NN2,M+1),
     *                      ZW1(NN2,M+1)
                       WRITE(5,1015) XW(NN2,M),YW1(NN2,M),
     *                      ZW1(NN2,M)
                    END DO
                 END IF
              END DO
              
C............Plot submerged panels on the hub.
              DO M=1,MHBT
                 DO N=1,NHBX
                    I=INDEXH(N,M)
                    IF(ISUBM(I,KK1).EQ.1) THEN
                       M1=M+1
                       N1=N+1
                       WRITE(5,1015) XH(N,M1),YH1(N,M1),ZH1(N,M1)
                       WRITE(5,1015) XH(N,M),YH1(N,M),ZH1(N,M)
                       WRITE(5,1015) XH(N1,M),YH1(N1,M),ZH1(N1,M)
                       WRITE(5,1015) XH(N1,M1),YH1(N1,M1),ZH1(N1,M1)
                    END IF
                 END DO
              END DO

           END DO
              
           IF(NEL1.GT.0) THEN
              DO II=1,NEL1+1
                 NI1=(II-1)*4
                 WRITE(5,1020) NI1+1,NI1+2,NI1+3,NI1+4
              END DO
           END IF
       
        END DO
        CLOSE(5)

        RETURN
        END 


