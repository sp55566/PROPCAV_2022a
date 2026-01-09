C -----------------------------------------
      SUBROUTINE PROP2GBFLOW(AKT,AKQ)
C -----------------------------------------
C
C     This routine writes info. for GBFLOW input
C     Blade geometry and mean dipole and sources
C ------------------------------------------------

      USE TOGBFLOW
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
C      DIMENSION STEMP(NPANZ), PTMEAN(NPBZ),SRCMEAN(NPBZ)
C      DIMENSION WKTMP(NPAWZ), WPMEAN(NPAWZ)
C      COMMON /TOGBFLOW/ XCCC(NBHZP,MBPZ),YCCC(NBHZP,MBPZ),
C     %              ZCCC(NBHZP,MBPZ),CPBMEAN(NSTEP,NBHZP,MBPZ)
C      DIMENSION XCCP(NBHZP,MBPZ),YCCP(NBHZP,MBPZ),ZCCP(NBHZP,MBPZ)

      ALLOCATABLE :: STEMP(:), PTMEAN(:),SRCMEAN(:)
      ALLOCATABLE :: WKTMP(:), WPMEAN(:)
      ALLOCATABLE :: XCCP(:,:),YCCP(:,:),ZCCP(:,:)

CVV   
      ALLOCATE(STEMP(NPANZ), PTMEAN(NPBZ),SRCMEAN(NPBZ))
      ALLOCATE(WKTMP(NPAWZ), WPMEAN(NPAWZ))
      ALLOCATE(XCCP(NBHZP,MBPZ),YCCP(NBHZP,MBPZ),ZCCP(NBHZP,MBPZ))
CVV

      OPEN(565,FILE='prop2gb.dat',STATUS='UNKNOWN')

C --Write Blade Coordinate

      NREAD = NWMIN * MR
 
      WRITE(565,*) ' Blade '
      WRITE(565,*) NCP, MRP
      WRITE(565,*) ' N, M , X , Y, Z '
      DO N = 1, NCP
         DO M = 1 , MRP
            WRITE(565,*) N, M, XB(N,M),YB(N,M),ZB(N,M)
         ENDDO
      ENDDO

      WRITE(565,*) ' WAKE '
      WRITE(565,*) NWMIN+1, MRP
      WRITE(565,*) ' N, M , X , Y, Z '
      DO N = 1 , NWMIN+1
         DO M = 1 , MRP
            WRITE(565,*) N , M , XW(N,M),YW(N,M),ZW(N,M)
         ENDDO
      ENDDO

C -- Write Mean strength of Source and Dipole

      DO I = 1, NTPREV 
         PTMEAN(I) = 0.0
         SRCMEAN(I) = 0.0
      ENDDO

      DO I = 1 , NTPREV
         CALL READ2(45,I,STEMP,NPANEL)  
         DO J = 1 , NPANB
            PTMEAN(J) = PTMEAN(J) + STEMP(J) / NTPREV
         ENDDO
      ENDDO

      DO I = 1 , NTPREV
         CALL READ2(47,I,STEMP,NPANEL)  
         DO J = 1 , NPANB
            SRCMEAN(J) = SRCMEAN(J) + STEMP(J) /NTPREV
         ENDDO
      ENDDO

      WRITE(565,*) 'Dipole and Source strength from Blade'
      WRITE(565,*) NPANB
      DO I = 1 , NPANB
         WRITE(565,*) I , PTMEAN(I) , SRCMEAN(I)
      ENDDO


      DO I = 1, NTPREV
         CALL READ2(46,I,WKTMP,NREAD)  
         DO J = 1 , NPANB
            WPMEAN(J) = WPMEAN(J) + WKTMP(J) / NTPREV
         ENDDO
      ENDDO

      WRITE(565,*) 'Dipole strength from Wake'
      WRITE(565,*) NREAD
      DO I = 1 , NREAD
         WRITE(565,*) I , WPMEAN(I)
      ENDDO

      close(565)

C ---- Input to GBFLOW (Circulation distribution) ---
C
      OPEN(565,FILE='prop2gb.gam',STATUS='UNKNOWN')

      WRITE(565,*)
      WRITE(565,*) ADVCO,AKT,AKQ
      WRITE(565,*) MR, NBLADE
      DO I = 1 , MR
         WRITE(565,*) I , I
      ENDDO
      WRITE(565,*) 1
      WRITE(565,*) 
      
      DO I = 1 , NTPREV
         akt = xktv(I,1)*nblade
         akq = xktv(I,4)*nblade
         WRITE(565,*) 1, -TT(I), akt, akq
         WRITE(565,*) 1.0
      ENDDO

      CLOSE(565)

      OPEN(565,FILE='prop2gb.wbf',STATUS='UNKNOWN')

      WRITE(565,*) NH, MR
      WRITE(565,*) ((XCCC(N,M) , N=1, NH), M=1, MR)
      WRITE(565,*) ((YCCC(N,M) , N=1, NH), M=1, MR)
      WRITE(565,*) ((ZCCC(N,M) , N=1, NH), M=1, MR)
      DO I = 1 , NTPREV
         THET = -TT(I) * PI / 180.0
         DO N = 1 , NH
            DO M = 1 , MR
               YYY = YCON(N,M) * COS(THET) + ZCON(N,M) * SIN(THET)
               ZZZ = -YCON(N,M) * SIN(THET) + ZCON(N,M) * COS(THET)

               XCCP(N,M) = 0.5*CPBMEAN(I,N,M) * XCON(N,M)
               YCCP(N,M) = 0.5*CPBMEAN(I,N,M) * YYY
               ZCCP(N,M) = 0.5*CPBMEAN(I,N,M) * ZZZ
            ENDDO
         ENDDO
         WRITE(565,*) I
         WRITE(565,*) ((XCCP(N,M), N=1,NH),M=1,MR)
         WRITE(565,*) ((YCCP(N,M), N=1,NH),M=1,MR)
         WRITE(565,*) ((ZCCP(N,M), N=1,NH),M=1,MR)
      ENDDO
         
      CLOSE(565)


CVV   
      DEALLOCATE(STEMP,PTMEAN,SRCMEAN)
      DEALLOCATE(WKTMP,WPMEAN)
      DEALLOCATE(XCCP,YCCP,ZCCP)
CVV

      RETURN
      END
