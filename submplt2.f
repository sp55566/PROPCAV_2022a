       SUBROUTINE SUBMPLT2

       use m_PSXYZ
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       DIMENSION XN(4,3)
!      COMMON/PSXYZ/XPS(NSTEP,NPBZ,5),YPS(NSTEP,NPBZ,5),
!    *     ZPS(NSTEP,NPBZ,5),IPSN(NPBZ,NSTEP),NEL(NSTEP)

       if (.NOT. allocated(xps)) then
         allocate(XPS(NSTEP,NPBZ,5),YPS(NSTEP,NPBZ,5))
         allocate(ZPS(NSTEP,NPBZ,5),IPSN(NPBZ,NSTEP),NEL(NSTEP))
       end if
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
             WRITE(5,1015)  XUWDK,YFS,-1.2
             WRITE(5,1015)  XUWDK,YFS, 1.2             
             WRITE(5,1015)  -1.,YFS, 1.2             
          END IF

          DO K1=1,NBLADE
              KK1=NTPOS1(K1,IT)
              DT1=-DELTAT*FLOAT(KK1-1)
              
              NPS=0

              DO M=MR,1,-1
                 
C...............Plot submerged panels on the blade.
                 DO N=1,NC
                    I=INDEXB(N,M)
                    
C..................Fully submerged panels
                    IF(ISUBM(I,KK1).EQ.1) THEN
                       CALL NODE(1,N,M,XN)
                       DO N1=1,4
                          CALL ROTATE(XN(N1,2),XN(N1,3),YY,ZZ,DT1)
                          WRITE(5,1015) XN(N1,1),YY,ZZ
                       END DO
                       
C..................Partially submerged panels
                    ELSE IF(ISUBM(I,KK1).EQ.2) THEN
                       NPS=NPS+1
                       CALL PLTPSUBM(NPS,KK1)
                    
                    END IF

                 END DO
              
                 DO N=1,NWMIN
                    I=IDXWAK(N,M)
                 
C..................Fully submerged panels
                    IF(ISUWM(I,KK1).EQ.1) THEN
                       CALL NODE(3,N,M,XN)
                       DO N1=1,4
                          CALL ROTATE(XN(N1,2),XN(N1,3),YY,ZZ,DT1)
                          WRITE(5,1015) XN(N1,1),YY,ZZ
                       END DO
                          
C..................Partially submerged panels
                    ELSE IF(ISUWM(I,KK1).EQ.2) THEN
                       NPS=NPS+1
                       CALL PLTPSUBM(NPS,KK1)

                    END IF
                    
                 END DO

              END DO

c              GO TO 3000

C............Plot submerged panels on the hub
              DO N=1,NHBX
                 DO M=1,MHBT
                    I=INDEXH(N,M)

C..................Fully submerged panels
                    IF(ISUBM(I,KK1).EQ.1) THEN
                       CALL NODE(2,N,M,XN)
                       DO N1=1,4
                          CALL ROTATE(XN(N1,2),XN(N1,3),YY,ZZ,DT1)
                          WRITE(5,1015) XN(N1,1),YY,ZZ
                       END DO
                       
C...................Partially submerged panels
                    ELSE IF(ISUBM(I,KK1).EQ.2) THEN
                       NPS=NPS+1
                       CALL PLTPSUBM(NPS,KK1)

                    END IF
                    
                 END DO
              END DO

 3000         CONTINUE

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


        SUBROUTINE PLTPSUBM(NPS,KK1)

        use m_PSXYZ
        INCLUDE 'PUFCAV.INC'
!       COMMON/PSXYZ/XPS(NSTEP,NPBZ,5),YPS(NSTEP,NPBZ,5),
!    *       ZPS(NSTEP,NPBZ,5),IPSN(NPBZ,NSTEP),NEL(NSTEP)
        if (.NOT. allocated(xps)) then
          allocate(XPS(NSTEP,NPBZ,5),YPS(NSTEP,NPBZ,5))
          allocate(ZPS(NSTEP,NPBZ,5),IPSN(NPBZ,NSTEP),NEL(NSTEP))
        end if
        
 1015   FORMAT(3(1X,F14.8))
        DO NN=1,3
           WRITE(5,1015) XPS(KK1,NPS,NN),
     *          YPS(KK1,NPS,NN),ZPS(KK1,NPS,NN)
        END DO
        IF(IPSN(NPS,KK1).EQ.3) THEN
           WRITE(5,1015) XPS(KK1,NPS,1),
     *          YPS(KK1,NPS,1),ZPS(KK1,NPS,1)
        ELSE IF(IPSN(NPS,KK1).GT.3) THEN
           WRITE(5,1015) XPS(KK1,NPS,4),
     *          YPS(KK1,NPS,4),ZPS(KK1,NPS,4) 
           IF(IPSN(NPS,KK1).GT.4) THEN
              WRITE(5,1015) XPS(KK1,NPS,4),
     *             YPS(KK1,NPS,4),ZPS(KK1,NPS,4)
              WRITE(5,1015) XPS(KK1,NPS,5),
     *             YPS(KK1,NPS,5),ZPS(KK1,NPS,5) 
              WRITE(5,1015) XPS(KK1,NPS,1),
     *             YPS(KK1,NPS,1),ZPS(KK1,NPS,1) 
              WRITE(5,1015) XPS(KK1,NPS,4),
     *             YPS(KK1,NPS,4),ZPS(KK1,NPS,4) 
           END IF
        END IF
        
        RETURN
        END
