      MODULE m_rpan_infcof
        type t_rpan_infcof
          real::XVrpan(4)
          real::YVrpan(4)
          real::SQrpan(15)
          real::SIDErpan(4)
        end type t_rpan_infcof
      END MODULE m_rpan_infcof

      MODULE m_hypot_infcof
        type t_hypot_infcof
          real:: XM1hypot(3)
          real:: XM2hypot(3)
          real:: XM3hypot(3)
          real:: XM4hypot(3)
          real:: XMChypot(3)
        end type t_hypot_infcof
      END MODULE m_hypot_infcof

      SUBROUTINE INFCOF
******************************************************************************
C 14-08-2017  S.N.KIM  ----- OpenMP parallelized to enhance computing
C                            efficiency. It is expected to reduct a huge 
C                            amound of computing time when intense geometrical 
C                            operation is included, such as repaneling.
******************************************************************************

      use omp_lib
      use m_rpan_infcof
      use m_hypot_infcof
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
c XM YU 10/2011
      INCLUDE 'PUFCAVC.INC'
c XM YU 10/2011

      TYPE(t_rpan_infcof)::rpan_geo_data
      TYPE(t_hypot_infcof)::hypot_geo_data
      DIMENSION A_TMP_INF(NPANZ,NC,MR,KZ)
      DIMENSION TEMP1_TMP_INF(NPANZ,NC,MR,KZ)
      DIMENSION A_TMP_INF2(NPANZ,MHBT,NHBX,KZ)
      DIMENSION TEMP1_TMP_INF2(NPANZ,MHBT,NHBX,KZ)
      DIMENSION A_TMP_INF3(NPANZ,NDUCT,MDUCT,KZ)
      DIMENSION TEMP1_TMP_INF3(NPANZ,NDUCT,MDUCT,KZ)
      DIMENSION A_TMP_INF4(NPANZ,NAXT,MTUNEL,KZ)
      DIMENSION TEMP1_TMP_INF4(NPANZ,NAXT,MTUNEL,KZ)

C-----------------------------------------------------------------------
C     Initialization
C-----------------------------------------------------------------------

C.....FILE 41 (IUMA)  -- dipole inf. functions for each blade
C.....FILE 42 (IUMB)  -- source inf. functions for each blade
c      write(*,*) 'beta infcof'

      IUMA=41
      IUMB=42

      IDDUCT2 = NPANB + NPANH
      IDDUCT3 = NPANB + NPANH + NPAND

      IAAA1 = NPANB + NPANH + NPAND
      IAAA2 = NPANB + NPANH + NPAND + NPANTN

      sum1 = 0.0
      sum2 = 0.0

      CALL CLEAR(A,NPANZ)
      CALL CLEAR(TEMP1,NPANZ)

      IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2

!$omp parallel
!$omp&private(M,N,J,K,rpan_geo_data,hypot_geo_data,IMR0,ER1,ER2,ER3,ER4,ERRMIN,
!$omp&        KK,I,CTLT,XLOC,YLOC,ZLOC,IMR,IXYZ,FD,FS,FSX,FSY,FDX,FDY,FDZ,
!$omp&        XV12,XV13,KT,II,JJ)
!$omp&shared(MR,NC,XVP,YVP,SID,SS,NBLADE,NPANEL,ZERO,XCTP,XCT,DIR,TWOPI,PI,
!$omp&       ICON,WINGIMAG,ITUN,IAAA1,IAAA2,IDUCT,IDOPT,IDDUCT2,IDDUCT3,
!$omp&       IVISC,W,A_TMP_INF,TEMP1_TMP_INF,XB,YB,ZB,A,
!$omp&       IHUB,NHBX,MHBT,XH,YH,ZH,MDUCT,NDUCT,XD,YD,ZD,
!$omp&       NAXT,MTUNEL,XTUN,YTUN,ZTUN,AVL
!$omp&       A_TMP_INF2,A_TMP_INF3,A_TMP_INF4,
!$omp&       TEMP1_TMP_INF2,TEMP1_TMP_INF3,TEMP1_TMP_INF4)
!$omp do
      DO M=MR,1,-1
         DO N=1,NC
            J=INDEXB(N,M)
            DO K=1,4
               rpan_geo_data%XVrpan(K)=XVP(J,K)
               rpan_geo_data%YVrpan(K)=YVP(J,K)
               rpan_geo_data%SIDErpan(K)=SID(J,K)
            END DO
            DO K=1,15
               rpan_geo_data%SQrpan(K)=SS(J,K)
            END DO
            CTLT = CHRLEPS(J) 

            hypot_geo_data%XM1hypot(1)=XB(N,M)
            hypot_geo_data%XM1hypot(2)=YB(N,M)
            hypot_geo_data%XM1hypot(3)=ZB(N,M)
            hypot_geo_data%XM2hypot(1)=XB(N,M+1)
            hypot_geo_data%XM2hypot(2)=YB(N,M+1)
            hypot_geo_data%XM2hypot(3)=ZB(N,M+1)
            hypot_geo_data%XM3hypot(1)=XB(N+1,M+1)
            hypot_geo_data%XM3hypot(2)=YB(N+1,M+1)
            hypot_geo_data%XM3hypot(3)=ZB(N+1,M+1)
            hypot_geo_data%XM4hypot(1)=XB(N+1,M)
            hypot_geo_data%XM4hypot(2)=YB(N+1,M)
            hypot_geo_data%XM4hypot(3)=ZB(N+1,M)
            IMR0=0

          ER1=ABS(hypot_geo_data%XM2hypot(1)-hypot_geo_data%XM1hypot(1))
     &       +ABS(hypot_geo_data%XM2hypot(2)-hypot_geo_data%XM1hypot(2))
     &       +ABS(hypot_geo_data%XM2hypot(3)-hypot_geo_data%XM1hypot(3))
          ER2=ABS(hypot_geo_data%XM3hypot(1)-hypot_geo_data%XM2hypot(1))
     &       +ABS(hypot_geo_data%XM3hypot(2)-hypot_geo_data%XM2hypot(2))
     &       +ABS(hypot_geo_data%XM3hypot(3)-hypot_geo_data%XM2hypot(3))
          ER3=ABS(hypot_geo_data%XM4hypot(1)-hypot_geo_data%XM3hypot(1))
     &       +ABS(hypot_geo_data%XM4hypot(2)-hypot_geo_data%XM3hypot(2))
     &       +ABS(hypot_geo_data%XM4hypot(3)-hypot_geo_data%XM3hypot(3))
          ER4=ABS(hypot_geo_data%XM1hypot(1)-hypot_geo_data%XM4hypot(1))
     &       +ABS(hypot_geo_data%XM1hypot(2)-hypot_geo_data%XM4hypot(2))
     &       +ABS(hypot_geo_data%XM1hypot(3)-hypot_geo_data%XM4hypot(3))

            ERRMIN=AMIN1(ER1,ER2,ER3,ER4)

            IF(ERRMIN.LE.1.0E-6) THEN
               IMR0=1
            END IF

            DO KK=1,NBLADE
               DO I=1,NPANEL
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO K=1,3
                     XLOC=XLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                  END DO

                  IMR=IMR0
                  CALL RPAN_PAR3(XLOC,YLOC,ZLOC,CTLT,
     *                 FS,FD,FSX,FSY,FDX,FDY,FDZ,0,IMR,rpan_geo_data)
                  
                  IF(IMR.EQ.2) THEN
                     DO IXYZ=1,3
                        hypot_geo_data%XMChypot(IXYZ)=XCTP(I,IXYZ,KK)
                     END DO
                     CALL HYPOT_PAR2(FD,FS,hypot_geo_data)
                  END IF

                  IF(I.EQ.J .AND. KK.EQ.1) THEN
                     IF(IMR.EQ.1) THEN
                        FD=TWOPI
                     END IF
                     IF(FD.LT.0.) THEN
                        FD=4.0*PI+FD
                     END IF

                     IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN
                        FD=2.0*PI
                     END IF

                  ELSE
                     IF(ABS(FD).GT.6.28) THEN
                        FD=0.0
                     END IF

                  END IF

c IMAGE MODEL
                  if ((icon.eq.5).and.(wingimag.ne.0)) then
                     call rudim_par(j,i,kk,imr,xv12,xv13,hypot_geo_data)
                     fd=fd+xv12
                     fs=fs+xv13
                  endif
c IMAGE MODEL

                  IF(ITUN .NE. 0) THEN
                     IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                      FD = -FD
                      FS = -FS
                     ENDIF
                  ENDIF

                  IF(IDUCT .EQ. 1) THEN
                     IF(IDOPT .EQ. 1) THEN
                        IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                           FD = 0.0
                           FS = 0.0
                        ENDIF
                     ENDIF
                  ENDIF

c                  A(I)=FD
c                  TEMP1(I)=FS
                  A_TMP_INF(I,N,M,KK) = FD
                  TEMP1_TMP_INF(I,N,M,KK) = FS

                  if (IVISC.eq.1) then
c XM YU 10/2011 ------------------------------------------------------------
c for viscous calculation
                    if(kk.eq.1) then
                      kt=nc*(mr-m)+1
                      if(i.ge.kt.and.i.lt.kt+nc)then
                        if(j.ge.kt.and.j.lt.kt+nc)then
                          ii=i+1-kt
                          jj=j+1-kt
                          AVL(ii,jj,m)=A_TMP_INF(I,N,M,KK)
ccc..... add wake influence
                          if(n.eq.1)then
                            AVL(ii,jj,m)=avl(ii,jj,m)-w(i,m)
                          else if(n.eq.nc)then
                            AVL(ii,jj,m)=avl(ii,jj,m)+w(i,m)
                          endif
                        endif
                      endif
                    endif
                  end if
c XM YU 10/2011-----------------------------------------------------------

               END DO
c               CALL WRITE1(IUMA,A,NPANEL)
c               CALL WRITE1(IUMB,TEMP1,NPANEL)
            END DO
         END DO
      END DO
!$omp end do

!$omp master
      DO M = MR, 1, -1
       DO N = 1, NC
         DO KK = 1, NBLADE
            CALL WRITE1(IUMA,A_TMP_INF(1,N,M,KK),NPANEL)
            CALL WRITE1(IUMB,TEMP1_TMP_INF(1,N,M,KK),NPANEL)
          ENDDO
        ENDDO
      ENDDO
!$omp end master

C-----------------------------------------------------------------------
C     Compute influence coefficients due to the hub
C-----------------------------------------------------------------------

      IF(IHUB.NE.0) THEN

!$omp do
        DO N=1,NHBX
          DO M=1,MHBT
            J=INDEXH(N,M)

            DO K=1,4
              rpan_geo_data%XVrpan(K)=XVP(J,K)
              rpan_geo_data%YVrpan(K)=YVP(J,K)
              rpan_geo_data%SIDErpan(K)=SID(J,K)
            END DO
            DO K=1,15
              rpan_geo_data%SQrpan(K)=SS(J,K)
            END DO
            CTLT = CHRLEPS(J)

            hypot_geo_data%XM1hypot(1)=XH(N,M+1)
            hypot_geo_data%XM1hypot(2)=YH(N,M+1)
            hypot_geo_data%XM1hypot(3)=ZH(N,M+1)
            hypot_geo_data%XM2hypot(1)=XH(N,M)
            hypot_geo_data%XM2hypot(2)=YH(N,M)
            hypot_geo_data%XM2hypot(3)=ZH(N,M)
            hypot_geo_data%XM3hypot(1)=XH(N+1,M)
            hypot_geo_data%XM3hypot(2)=YH(N+1,M)
            hypot_geo_data%XM3hypot(3)=ZH(N+1,M)
            hypot_geo_data%XM4hypot(1)=XH(N+1,M+1)
            hypot_geo_data%XM4hypot(2)=YH(N+1,M+1)
            hypot_geo_data%XM4hypot(3)=ZH(N+1,M+1)
            IMR0=0

          ER1=ABS(hypot_geo_data%XM2hypot(1)-hypot_geo_data%XM1hypot(1))
     &       +ABS(hypot_geo_data%XM2hypot(2)-hypot_geo_data%XM1hypot(2))
     &       +ABS(hypot_geo_data%XM2hypot(3)-hypot_geo_data%XM1hypot(3))
          ER2=ABS(hypot_geo_data%XM3hypot(1)-hypot_geo_data%XM2hypot(1))
     &       +ABS(hypot_geo_data%XM3hypot(2)-hypot_geo_data%XM2hypot(2))
     &       +ABS(hypot_geo_data%XM3hypot(3)-hypot_geo_data%XM2hypot(3))
          ER3=ABS(hypot_geo_data%XM4hypot(1)-hypot_geo_data%XM3hypot(1))
     &       +ABS(hypot_geo_data%XM4hypot(2)-hypot_geo_data%XM3hypot(2))
     &       +ABS(hypot_geo_data%XM4hypot(3)-hypot_geo_data%XM3hypot(3))
          ER4=ABS(hypot_geo_data%XM1hypot(1)-hypot_geo_data%XM4hypot(1))
     &       +ABS(hypot_geo_data%XM1hypot(2)-hypot_geo_data%XM4hypot(2))
     &       +ABS(hypot_geo_data%XM1hypot(3)-hypot_geo_data%XM4hypot(3))

            ERRMIN=AMIN1(ER1,ER2,ER3,ER4)

            IF(ERRMIN.LE.1.0E-6) THEN
               WRITE(*,2000) M,N
 2000          FORMAT(' ...... triangular panels (m,n) -- > ',
     *              2I4,' (hub) ')
               IMR0=1
            END IF

            DO KK=1,NBLADE
               DO I=1,NPANEL
                XLOC=0.
                YLOC=0.
                ZLOC=0.
                DO K=1,3
                  XLOC=XLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                  YLOC=YLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                  ZLOC=ZLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                END DO

                IMR=IMR0

                CALL RPAN_PAR3(XLOC,YLOC,ZLOC,CTLT,FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR,rpan_geo_data)

                IF(IMR.EQ.2) THEN
                   DO IXYZ=1,3
                      hypot_geo_data%XMChypot(IXYZ)=XCTP(I,IXYZ,KK)
                   END DO
                   CALL HYPOT_PAR2(FD,FS,hypot_geo_data)
                END IF

                IF(I.EQ.J .AND. KK.EQ.1) THEN
                   IF(IMR.EQ.1) THEN
                      FD=TWOPI
                   END IF
                   IF(FD.LT.0.) THEN
                      FD=4.0*PI+FD
                   END IF

                   IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.5.0) THEN
                      FD=2.0*PI
                   END IF
                ELSE
                   IF(ABS(FD).GT.6.28) THEN
                      FD=0.0
                   END IF

                END IF

                IF(ITUN .NE. 0) THEN
                   IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                      FD = -FD
                      FS = -FS
                   ENDIF
                ENDIF

                IF(IDUCT .EQ. 1) THEN
                   IF(IDOPT .EQ. 1) THEN
                      IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                         FD = 0.0
                         FS = 0.0
                      ENDIF

                   ENDIF
                ENDIF

                A_TMP_INF2(I,M,N,KK) = FD
                TEMP1_TMP_INF2(I,M,N,KK) = FS

              END DO
c              CALL WRITE1(IUMA,A,NPANEL)
c              CALL WRITE1(IUMB,TEMP1,NPANEL)
            END DO
          END DO
        END DO
!$omp end do 

!$omp master
      DO N = 1, NHBX
        DO M = 1, MHBT
          DO KK = 1, NBLADE  
            CALL WRITE1(IUMA,A_TMP_INF2(1,M,N,KK),NPANEL)
            CALL WRITE1(IUMB,TEMP1_TMP_INF2(1,M,N,KK),NPANEL)
          ENDDO
        ENDDO
      ENDDO
!$omp end master 

      END IF

C-----------------------------------------------------------------------
C     Compute influence coefficients due to DUCT
C-----------------------------------------------------------------------

      IF(IDUCT .NE. 0) THEN
!$omp do
         DO M=1,MDUCT
            DO N=1,NDUCT
               J=INDEXD(N,M)

               DO K=1,4
                  rpan_geo_data%XVrpan(K)=XVP(J,K)
                  rpan_geo_data%YVrpan(K)=YVP(J,K)
                  rpan_geo_data%SIDErpan(K)=SID(J,K)
               ENDDO

               DO K=1,15
                  rpan_geo_data%SQrpan(K)=SS(J,K)
               ENDDO
               CTLT = CHRLEPS(J)

               hypot_geo_data%XM1hypot(1)=XD(N,M+1)
               hypot_geo_data%XM1hypot(2)=YD(N,M+1)
               hypot_geo_data%XM1hypot(3)=ZD(N,M+1)
               hypot_geo_data%XM2hypot(1)=XD(N,M)
               hypot_geo_data%XM2hypot(2)=YD(N,M)
               hypot_geo_data%XM2hypot(3)=ZD(N,M)
               hypot_geo_data%XM3hypot(1)=XD(N+1,M)
               hypot_geo_data%XM3hypot(2)=YD(N+1,M)
               hypot_geo_data%XM3hypot(3)=ZD(N+1,M)
               hypot_geo_data%XM4hypot(1)=XD(N+1,M+1)
               hypot_geo_data%XM4hypot(2)=YD(N+1,M+1)
               hypot_geo_data%XM4hypot(3)=ZD(N+1,M+1)

               IMR0=0

          ER1=ABS(hypot_geo_data%XM2hypot(1)-hypot_geo_data%XM1hypot(1))
     &       +ABS(hypot_geo_data%XM2hypot(2)-hypot_geo_data%XM1hypot(2))
     &       +ABS(hypot_geo_data%XM2hypot(3)-hypot_geo_data%XM1hypot(3))
          ER2=ABS(hypot_geo_data%XM3hypot(1)-hypot_geo_data%XM2hypot(1))
     &       +ABS(hypot_geo_data%XM3hypot(2)-hypot_geo_data%XM2hypot(2))
     &       +ABS(hypot_geo_data%XM3hypot(3)-hypot_geo_data%XM2hypot(3))
          ER3=ABS(hypot_geo_data%XM4hypot(1)-hypot_geo_data%XM3hypot(1))
     &       +ABS(hypot_geo_data%XM4hypot(2)-hypot_geo_data%XM3hypot(2))
     &       +ABS(hypot_geo_data%XM4hypot(3)-hypot_geo_data%XM3hypot(3))
          ER4=ABS(hypot_geo_data%XM1hypot(1)-hypot_geo_data%XM4hypot(1))
     &       +ABS(hypot_geo_data%XM1hypot(2)-hypot_geo_data%XM4hypot(2))
     &       +ABS(hypot_geo_data%XM1hypot(3)-hypot_geo_data%XM4hypot(3))

               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)

               IF(ERRMIN.LE.1.0E-6) THEN
                  WRITE(*,2100) M,N
 2100             FORMAT(' ...... triangular panels (m,n) -- > ',
     *                 2I4,' (DUCT) ')
                  IMR0=1
               END IF

               DO KK=1,NBLADE
                  DO I=1,NPANEL
                     XLOC=0.
                     YLOC=0.
                     ZLOC=0.
                     DO K=1,3
                        XLOC=XLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                        YLOC=YLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                        ZLOC=ZLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                     ENDDO
                     IMR=IMR0

                     CALL RPAN_PAR3(XLOC,YLOC,ZLOC,CTLT,FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR,rpan_geo_data)

                     IF(IMR.EQ.2) THEN
                        DO IXYZ=1,3
                           hypot_geo_data%XMChypot(IXYZ)=XCTP(I,IXYZ,KK)
                        ENDDO
                        CALL HYPOT_PAR2(FD,FS,hypot_geo_data)
                     END IF

                     IF(I.EQ.J .AND. KK.EQ.1) THEN
                        IF(IMR.EQ.1) THEN
                           FD=TWOPI
                        END IF
                        IF(FD.LT.0.) THEN
                           FD=4.0*PI+FD
                        END IF
                        IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.5.0) THEN
                           FD=2.0*PI
                        END IF
                     ELSE

                        IF(ABS(FD).GT.6.28) THEN
                           FD=0.0
                        END IF

                     END IF

                     IF(ITUN .NE. 0) THEN
                        IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                          FD = -FD
                          FS = -FS
                        ENDIF
                    ENDIF

                    IF(IDOPT .EQ. 1) THEN
                       IF(I .LE. IDDUCT2 .OR. I .GT. IDDUCT3) THEN
                          FD = 0.0
                          FS = 0.0
                       ENDIF

                    ENDIF

                    A_TMP_INF3(I,N,M,KK)=FD
                    TEMP1_TMP_INF3(I,N,M,KK) = FS

                  ENDDO
c                  CALL WRITE1(IUMA,A,NPANEL)
c                  CALL WRITE1(IUMB,TEMP1,NPANEL)
               ENDDO

            ENDDO
         ENDDO
!$omp end do

!$omp master
         DO M = 1, MDUCT
           DO N = 1, NDUCT
             DO KK = 1, NBLADE
               CALL WRITE1(IUMA,A_TMP_INF3(1,N,M,KK),NPANEL)
               CALL WRITE1(IUMB,TEMP1_TMP_INF3(1,N,M,KK),NPANEL)
             ENDDO
           ENDDO
         ENDDO 
!$omp end master 

      END IF

C-----------------------------------------------------------------------
C     Compute influence coefficients due to the tunnel walls
C-----------------------------------------------------------------------

      IF(ITUN.NE.0) THEN
!$omp do
         DO N=1,NAXT
            DO M=1,MTUNEL
               J=INDEXTN(N,M)

               DO K=1,4
                  rpan_geo_data%XVrpan(K)=XVP(J,K)
                  rpan_geo_data%YVrpan(K)=YVP(J,K)
                  rpan_geo_data%SIDErpan(K)=SID(J,K)
               ENDDO

               DO K=1,15
                  rpan_geo_data%SQrpan(K)=SS(J,K)
               ENDDO
               CTLT = CHRLEPS(J)

               hypot_geo_data%XM1hypot(1)=XTUN(N,M)
               hypot_geo_data%XM1hypot(2)=YTUN(N,M)
               hypot_geo_data%XM1hypot(3)=ZTUN(N,M)
               hypot_geo_data%XM2hypot(1)=XTUN(N,M+1)
               hypot_geo_data%XM2hypot(2)=YTUN(N,M+1)
               hypot_geo_data%XM2hypot(3)=ZTUN(N,M+1)
               hypot_geo_data%XM3hypot(1)=XTUN(N+1,M+1)
               hypot_geo_data%XM3hypot(2)=YTUN(N+1,M+1)
               hypot_geo_data%XM3hypot(3)=ZTUN(N+1,M+1)
               hypot_geo_data%XM4hypot(1)=XTUN(N+1,M)
               hypot_geo_data%XM4hypot(2)=YTUN(N+1,M)
               hypot_geo_data%XM4hypot(3)=ZTUN(N+1,M)

               IMR0=0

          ER1=ABS(hypot_geo_data%XM2hypot(1)-hypot_geo_data%XM1hypot(1))
     &       +ABS(hypot_geo_data%XM2hypot(2)-hypot_geo_data%XM1hypot(2))
     &       +ABS(hypot_geo_data%XM2hypot(3)-hypot_geo_data%XM1hypot(3))
          ER2=ABS(hypot_geo_data%XM3hypot(1)-hypot_geo_data%XM2hypot(1))
     &       +ABS(hypot_geo_data%XM3hypot(2)-hypot_geo_data%XM2hypot(2))
     &       +ABS(hypot_geo_data%XM3hypot(3)-hypot_geo_data%XM2hypot(3))
          ER3=ABS(hypot_geo_data%XM4hypot(1)-hypot_geo_data%XM3hypot(1))
     &       +ABS(hypot_geo_data%XM4hypot(2)-hypot_geo_data%XM3hypot(2))
     &       +ABS(hypot_geo_data%XM4hypot(3)-hypot_geo_data%XM3hypot(3))
          ER4=ABS(hypot_geo_data%XM1hypot(1)-hypot_geo_data%XM4hypot(1))
     &       +ABS(hypot_geo_data%XM1hypot(2)-hypot_geo_data%XM4hypot(2))
     &       +ABS(hypot_geo_data%XM1hypot(3)-hypot_geo_data%XM4hypot(3))

               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)

               IF(ERRMIN.LE.1.0E-6) THEN
                  WRITE(*,2200) M,N
 2200             FORMAT(' ...... triangular panels (m,n) -- > ',
     *                 2I4,' (tunnel) ')
                  IMR0=1
               END IF

               DO KK=1,NBLADE
                  DO I=1,NPANEL
                     XLOC=0.
                     YLOC=0.
                     ZLOC=0.
                     DO K=1,3
                        XLOC=XLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                        YLOC=YLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                        ZLOC=ZLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                     ENDDO
                     IMR=IMR0

                     CALL RPAN_PAR3(XLOC,YLOC,ZLOC,CTLT,FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR,rpan_geo_data)

                     IF(IMR.EQ.2) THEN
                        DO IXYZ=1,3
                           hypot_geo_data%XMChypot(IXYZ)=XCTP(I,IXYZ,KK)
                        ENDDO
                        CALL HYPOT_PAR2(FD,FS,hypot_geo_data)
                     END IF

                     IF(I.EQ.J .AND. KK.EQ.1) THEN
                        IF(IMR.EQ.1) THEN
                           FD=TWOPI
                        END IF
                        IF(FD.LT.0.) THEN
                           FD=4.0*PI+FD
                        END IF
                        IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.5.0) THEN
                           FD=2.0*PI
                        END IF
                     ELSE
                        IF(ABS(FD).GT.6.28) THEN
                           FD=0.0
                        END IF
                        FD = -FD
                     END IF

                     IF(IDUCT .NE. 0.AND. IDOPT .EQ. 1) THEN
                        IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                           FD = 0.0
                           FS = 0.0
                        ENDIF
                     ENDIF

                     A_TMP_INF4(I,N,M,KK)=FD
                     TEMP1_TMP_INF4(I,N,M,KK) = FS

                  ENDDO
c                  CALL WRITE1(IUMA,A,NPANEL)
c                  CALL WRITE1(IUMB,TEMP1,NPANEL)
               ENDDO

            ENDDO
         ENDDO
!$omp end do
 
!$omp master
         DO N = 1, NAXT
           DO M = 1, MTUNEL
             DO KK = 1, NBLADE
               CALL WRITE1(IUMA,A_TMP_INF4(1,N,M,KK),NPANEL)
               CALL WRITE1(IUMB,TEMP1_TMP_INF4(1,N,M,KK),NPANEL)
             ENDDO
           ENDDO
         ENDDO 
!$omp end master 

      END IF
!$omp end parallel

cC-----------------------------------------------------------------------
cC     Compute influence coefficients due to the tip vortex cavity
cC     Do not include effects of the tip vortex and tip bulb
cC     corresponding to the other blades.
cC     (Consider only KEY Blade's effect)
cC-----------------------------------------------------------------------
c
c      IF(IAN.EQ.2) THEN
c
c         DO N=1,nthx
c            DO M=1,mcvt
c
c               J=INDEXT(N,M)
c
c               IMR0=0
c
c               DO K=1,4
c                  XV(K)=XVP(J,K)
c                  YV(K)=YVP(J,K)
c                  SIDE(K)=SID(J,K)
c               ENDDO
c               DO K=1,15
c                  S(K)=SS(J,K)
c               ENDDO
c
c               XM1(1)=XCH(N+1,M+1)
c               XM1(2)=YCH(N+1,M+1)
c               XM1(3)=ZCH(N+1,M+1)
c               XM2(1)=XCH(N+1,M)
c               XM2(2)=YCH(N+1,M)
c               XM2(3)=ZCH(N+1,M)
c               XM3(1)=XCH(N,M)
c               XM3(2)=YCH(N,M)
c               XM3(3)=ZCH(N,M)
c               XM4(1)=XCH(N,M+1)
c               XM4(2)=YCH(N,M+1)
c               XM4(3)=ZCH(N,M+1)
c
c               CRLT = chrleps(j)
c
c               ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
c     %              +ABS(XM2(3)-XM1(3))
c               ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
c     %              +ABS(XM3(3)-XM2(3))
c               ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
c     %              +ABS(XM4(3)-XM3(3))
c               ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
c     %              +ABS(XM1(3)-XM4(3))
c               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)
c               IF(ERRMIN.LE.1.0E-6) THEN
c                  IMR0=1
c               END IF
c
c               DO I=1,NPANEL
c                  XLOC=0.
c                  YLOC=0.
c                  ZLOC=0.
c
c                  DO K=1,3
c                     XLOC=XLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,1,K)
c                     YLOC=YLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,2,K)
c                     ZLOC=ZLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,3,K)
c                  ENDDO
c
c                  IMR=IMR0
c
c                  CALL RPAN(XLOC,YLOC,ZLOC,CRLT,FS,FD,FSX,FSY,
c     *                 FDX,FDY,FDZ,0,IMR)
c
c                  IF(IMR.EQ.2) THEN
c                     DO IXYZ=1,3
c                        XMC(IXYZ)=XCTP(I,IXYZ,1)
c                     ENDDO
c                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
c                  END IF
c
c                  IF(I.EQ.J) THEN
c                     IF(IMR.EQ.1) THEN
c                        FD=TWOPI
c                     END IF
c                     IF(FD.LT.0.) THEN
c                        FD=4.0*PI+FD
c                     END IF
c
c                     IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN
c                        FD=2.0*PI
c                     END IF
c                  ELSE
c                     IF(ABS(FD).GT.6.28) THEN
c                        FD=0.0
c                     END IF
c                  END IF
c
c                  IF(IDUCT .NE. 0 .AND. IDOPT .EQ. 1) THEN
c                     IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
c                        FD = 0.0
c                        FS = 0.0
c                     ENDIF
c                  ENDIF
c
c                  A(I)=FD
c                  TEMP1(I)=FS
c
c               ENDDO
c               CALL WRITE1(IUMA,A,NPANEL)
c               CALL WRITE1(IUMB,TEMP1,NPANEL)
c            ENDDO
c         ENDDO
c
c         DO N=1,ncvx
c            DO M=1,mcvt
c
c               J=INDEXC(N,M)
c
c               IMR0 = 0
c
c               DO K=1,4
c                  XV(K)=XVP(J,K)
c                  YV(K)=YVP(J,K)
c                  SIDE(K)=SID(J,K)
c               ENDDO
c               DO K=1,15
c                  S(K)=SS(J,K)
c               ENDDO
c
c               XM1(1)=XVC(N+1,M+1)
c               XM1(2)=YVC(N+1,M+1)
c               XM1(3)=ZVC(N+1,M+1)
c               XM2(1)=XVC(N+1,M)
c               XM2(2)=YVC(N+1,M)
c               XM2(3)=ZVC(N+1,M)
c               XM3(1)=XVC(N,M)
c               XM3(2)=YVC(N,M)
c               XM3(3)=ZVC(N,M)
c               XM4(1)=XVC(N,M+1)
c               XM4(2)=YVC(N,M+1)
c               XM4(3)=ZVC(N,M+1)
c
c               ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
c     %              +ABS(XM2(3)-XM1(3))
c               ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
c     %              +ABS(XM3(3)-XM2(3))
c               ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
c     %              +ABS(XM4(3)-XM3(3))
c               ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
c     %              +ABS(XM1(3)-XM4(3))
c               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)
c               IF(ERRMIN.LE.1.0E-6) THEN
c                  IMR0=1
c               END IF
c
c               CRLT = chrleps(j)
c
c               DO I=1,NPANEL
c                  XLOC=0.
c                  YLOC=0.
c                  ZLOC=0.
c                  DO K=1,3
c                     XLOC=XLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,1,K)
c                     YLOC=YLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,2,K)
c                     ZLOC=ZLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,3,K)
c                  ENDDO
c
c                  IMR=IMR0
c
c                  CALL RPAN(XLOC,YLOC,ZLOC,CRLT,FS,FD,FSX,FSY,
c     *                 FDX,FDY,FDZ,0,IMR)
c
c                  IF(IMR.EQ.2) THEN
c                     DO IXYZ=1,3
c                        XMC(IXYZ)=XCTP(I,IXYZ,1)
c                     ENDDO
c                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
c                  END IF
c                  IF(I.EQ.J) THEN
c                     IF(IMR.EQ.1) THEN
c                        FD=TWOPI
c                     END IF
c                     IF(FD.LT.0.) THEN
c                        FD=4.0*PI+FD
c                     END IF
c
c                     IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN
c                        FD=2.0*PI
c                     END IF
c                  ELSE
c                     IF(ABS(FD).GT.6.28) THEN
c                        FD=0.0
c                     END IF
c                  END IF
c                  A(I)=FD
c                  TEMP1(I)=FS
c               ENDDO
c
c               IF(IDUCT .NE. 0 .AND. IDOPT .EQ. 1) THEN
c                  IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
c                     FD = 0.0
c                     FS = 0.0
c                  ENDIF
c               ENDIF
c
c               CALL WRITE1(IUMA,A,NPANEL)
c               CALL WRITE1(IUMB,TEMP1,NPANEL)
c
c            ENDDO
c         ENDDO
c      ENDIF

CSH--REPLACE, NBLADE=1-------------------------
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------

      RETURN
      END SUBROUTINE

      SUBROUTINE INFCOF2
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
c XM YU 10/2011
      INCLUDE 'PUFCAVC.INC'
c XM YU 10/2011

      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)

C-----------------------------------------------------------------------
C     Initialization
C-----------------------------------------------------------------------

C.....FILE 41 (IUMA)  -- dipole inf. functions for each blade
C.....FILE 42 (IUMB)  -- source inf. functions for each blade

      IUMA=41
      IUMB=42

      IDDUCT2 = NPANB + NPANH
      IDDUCT3 = NPANB + NPANH + NPAND

      IAAA1 = NPANB + NPANH + NPAND
      IAAA2 = NPANB + NPANH + NPAND + NPANTN

      sum1 = 0.0
      sum2 = 0.0

      CALL CLEAR(A,NPANZ)
      CALL CLEAR(TEMP1,NPANZ)

      IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2

      DO M=MR,1,-1
         DO N=1,NC
            J=INDEXB(N,M)
            DO K=1,4
               XV(K)=XVP(J,K)
               YV(K)=YVP(J,K)
               SIDE(K)=SID(J,K)
            END DO
            DO K=1,15
               S(K)=SS(J,K)
            END DO

            XM1(1)=XB(N,M)
            XM1(2)=YB(N,M)
            XM1(3)=ZB(N,M)
            XM2(1)=XB(N,M+1)
            XM2(2)=YB(N,M+1)
            XM2(3)=ZB(N,M+1)
            XM3(1)=XB(N+1,M+1)
            XM3(2)=YB(N+1,M+1)
            XM3(3)=ZB(N+1,M+1)
            XM4(1)=XB(N+1,M)
            XM4(2)=YB(N+1,M)
            XM4(3)=ZB(N+1,M)

            DO KK=1,NBLADE

!$OMP parallel private(I,K,XMC,FD,FS,XV12,XV13)
!$OMP&shared(KK,NPANEL,XCTP,XM1,XM2,XM3,XM4,J,PI,icon,wingimag,
!$OMP&       imr,ITUN,IAAA1,IAAA2,IDUCT,IDOPT,IDDUCT2,IDDUCT3,A,TEMP1)
!$OMP DO SCHEDULE(DYNAMIC,50)
               DO I=1,NPANEL
                  DO K=1,3
                     XMC(K)=XCTP(I,K,KK)
                  END DO
                  CALL ANRPAN(XM1,XM2,XM3,XM4,XMC,FD,FS)

                  IF(I.EQ.J .AND. KK.EQ.1) THEN
                     IF(FD.LT.0.) THEN
                        FD=4.0*PI+FD
                     END IF
                  END IF
c IMAGE MODEL
                  if ((icon.eq.5).and.(wingimag.ne.0)) then
                     call rudim(j,i,kk,imr,xm1,xm2,xm3,xm4,xv12,xv13)
                     fd=fd+xv12
                     fs=fs+xv13
                  endif
c IMAGE MODEL

                  IF(ITUN .NE. 0) THEN
                     IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                      FD = -FD
                      FS = -FS
                     ENDIF
                  ENDIF

                  IF(IDUCT .EQ. 1) THEN
                     IF(IDOPT .EQ. 1) THEN
                        IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                           FD = 0.0
                           FS = 0.0
                        ENDIF
                     ENDIF
                  ENDIF
                  A(I)=FD
                  TEMP1(I)=FS
               END DO
!$OMP END DO
!$OMP end parallel

               DO I=1,NPANEL
                 if (IVISC.eq.1) then
c XM YU 10/2011 ------------------------------------------------------------
c for viscous calculation
                   if(kk.eq.1) then
                     kt=nc*(mr-m)+1
                     if(i.ge.kt.and.i.lt.kt+nc)then
                       if(j.ge.kt.and.j.lt.kt+nc)then
                         ii=i+1-kt
                         jj=j+1-kt
                         AVL(ii,jj,m)=A(I)
ccc..... add wake influence
                         if(n.eq.1)then
                           AVL(ii,jj,m)=avl(ii,jj,m)-w(i,m)
                         else if(n.eq.nc)then
                           AVL(ii,jj,m)=avl(ii,jj,m)+w(i,m)
                         endif
                      endif
                    endif
                  endif
                end if
c XM YU 10/2011-----------------------------------------------------------
              END DO


               CALL WRITE1(IUMA,A,NPANEL)
               CALL WRITE1(IUMB,TEMP1,NPANEL)
            END DO
         END DO
      END DO

      WRITE(*,*) "Finished calculating blade inf-cof with NRPAN"

C-----------------------------------------------------------------------
C     Compute influence coefficients due to the hub
C-----------------------------------------------------------------------

      IF(IHUB.NE.0) THEN

        DO N=1,NHBX
          DO M=1,MHBT
            J=INDEXH(N,M)

            DO K=1,4
              XV(K)=XVP(J,K)
              YV(K)=YVP(J,K)
              SIDE(K)=SID(J,K)
            END DO
            DO K=1,15
              S(K)=SS(J,K)
            END DO

            XM1(1)=XH(N,M+1)
            XM1(2)=YH(N,M+1)
            XM1(3)=ZH(N,M+1)
            XM2(1)=XH(N,M)
            XM2(2)=YH(N,M)
            XM2(3)=ZH(N,M)
            XM3(1)=XH(N+1,M)
            XM3(2)=YH(N+1,M)
            XM3(3)=ZH(N+1,M)
            XM4(1)=XH(N+1,M+1)
            XM4(2)=YH(N+1,M+1)
            XM4(3)=ZH(N+1,M+1)

            DO KK=1,NBLADE

!$OMP parallel private(I,K,XMC,FD,FS)
!$OMP&shared(KK,NPANEL,XCTP,XM1,XM2,XM3,XM4,J,PI,
!$OMP&       ITUN,IAAA1,IAAA2,IDUCT,IDOPT,IDDUCT2,IDDUCT3,A,TEMP1)
!$OMP DO SCHEDULE(DYNAMIC,50)
               DO I=1,NPANEL
                  DO K=1,3
                     XMC(K)=XCTP(I,K,KK)
                  END DO
                  CALL ANRPAN(XM1,XM2,XM3,XM4,XMC,FD,FS)

                IF(I.EQ.J .AND. KK.EQ.1) THEN
                   IF(FD.LT.0.) THEN
                      FD=4.0*PI+FD
                   END IF
                END IF

                IF(ITUN .NE. 0) THEN
                   IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                      FD = -FD
                      FS = -FS
                   ENDIF
                ENDIF

                IF(IDUCT .EQ. 1) THEN
                   IF(IDOPT .EQ. 1) THEN
                      IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                         FD = 0.0
                         FS = 0.0
                      ENDIF

                   ENDIF
                ENDIF

                A(I)=FD
                TEMP1(I)=FS

              END DO
!$OMP END DO
!$OMP end parallel


              CALL WRITE1(IUMA,A,NPANEL)
              CALL WRITE1(IUMB,TEMP1,NPANEL)
            END DO
          END DO
        END DO
      WRITE(*,*) "Finished calculating hub inf-cof with NRPAN"
      END IF


C-----------------------------------------------------------------------
C     Compute influence coefficients due to DUCT
C-----------------------------------------------------------------------

      IF(IDUCT .NE. 0) THEN

         DO M=1,MDUCT
            DO N=1,NDUCT
               J=INDEXD(N,M)

               DO K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
               ENDDO

               DO K=1,15
                  S(K)=SS(J,K)
               ENDDO

               XM1(1)=XD(N,M+1)
               XM1(2)=YD(N,M+1)
               XM1(3)=ZD(N,M+1)
               XM2(1)=XD(N,M)
               XM2(2)=YD(N,M)
               XM2(3)=ZD(N,M)
               XM3(1)=XD(N+1,M)
               XM3(2)=YD(N+1,M)
               XM3(3)=ZD(N+1,M)
               XM4(1)=XD(N+1,M+1)
               XM4(2)=YD(N+1,M+1)
               XM4(3)=ZD(N+1,M+1)

               DO KK=1,NBLADE
                  DO I=1,NPANEL
                    DO K=1,3
                       XMC(K)=XCTP(I,K,KK)
                    END DO
                    CALL ANRPAN(XM1,XM2,XM3,XM4,XMC,FD,FS)

                    IF(I.EQ.J .AND. KK.EQ.1) THEN
                        IF(FD.LT.0.) THEN
                           FD=4.0*PI+FD
                        END IF
                     END IF

                     IF(ITUN .NE. 0) THEN
                        IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                          FD = -FD
                          FS = -FS
                        ENDIF
                    ENDIF

                    IF(IDOPT .EQ. 1) THEN
                       IF(I .LE. IDDUCT2 .OR. I .GT. IDDUCT3) THEN
                          FD = 0.0
                          FS = 0.0
                       ENDIF
                    ENDIF

                    A(I)=FD
                    TEMP1(I) = FS

                  ENDDO
                  CALL WRITE1(IUMA,A,NPANEL)
                  CALL WRITE1(IUMB,TEMP1,NPANEL)
               ENDDO

            ENDDO
         ENDDO
         WRITE(*,*) "Finished calculating duct inf-cof with NRPAN"
      END IF


C-----------------------------------------------------------------------
C     Compute influence coefficients due to the tunnel walls
C-----------------------------------------------------------------------

      IF(ITUN.NE.0) THEN

         DO N=1,NAXT
            DO M=1,MTUNEL
               J=INDEXTN(N,M)

               DO K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
               ENDDO

               DO K=1,15
                  S(K)=SS(J,K)
               ENDDO

               XM1(1)=XTUN(N,M)
               XM1(2)=YTUN(N,M)
               XM1(3)=ZTUN(N,M)
               XM2(1)=XTUN(N,M+1)
               XM2(2)=YTUN(N,M+1)
               XM2(3)=ZTUN(N,M+1)
               XM3(1)=XTUN(N+1,M+1)
               XM3(2)=YTUN(N+1,M+1)
               XM3(3)=ZTUN(N+1,M+1)
               XM4(1)=XTUN(N+1,M)
               XM4(2)=YTUN(N+1,M)
               XM4(3)=ZTUN(N+1,M)


               DO KK=1,NBLADE
                  DO I=1,NPANEL
                    DO K=1,3
                       XMC(K)=XCTP(I,K,KK)
                    ENDDO
                    CALL ANRPAN(XM1,XM2,XM3,XM4,XMC,FD,FS)

                     IF(I.EQ.J .AND. KK.EQ.1) THEN
                        IF(FD.LT.0.) THEN
                           FD=4.0*PI+FD
                        END IF
                     ELSE
                        FD = -FD
                     END IF

                     IF(IDUCT .NE. 0.AND. IDOPT .EQ. 1) THEN
                        IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                           FD = 0.0
                           FS = 0.0
                        ENDIF
                     ENDIF

                     A(I)=FD
                     TEMP1(I) = FS

                  ENDDO
                  CALL WRITE1(IUMA,A,NPANEL)
                  CALL WRITE1(IUMB,TEMP1,NPANEL)
               ENDDO

            ENDDO
         ENDDO
         WRITE(*,*) "Finished calculating tunnel inf-cof with NRPAN"
      END IF

C************************************************************************
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
C************************************************************************
cC-----------------------------------------------------------------------
cC     Compute influence coefficients due to the tip vortex cavity
cC     Do not include effects of the tip vortex and tip bulb
cC     corresponding to the other blades.
cC     (Consider only KEY Blade's effect)
cC-----------------------------------------------------------------------
c
c      IF(IAN.EQ.2) THEN
c
c         DO N=1,nthx
c            DO M=1,mcvt
c
c               J=INDEXT(N,M)
c
c               IMR0=0
c
c               DO K=1,4
c                  XV(K)=XVP(J,K)
c                  YV(K)=YVP(J,K)
c                  SIDE(K)=SID(J,K)
c               ENDDO
c               DO K=1,15
c                  S(K)=SS(J,K)
c               ENDDO
c
c               XM1(1)=XCH(N+1,M+1)
c               XM1(2)=YCH(N+1,M+1)
c               XM1(3)=ZCH(N+1,M+1)
c               XM2(1)=XCH(N+1,M)
c               XM2(2)=YCH(N+1,M)
c               XM2(3)=ZCH(N+1,M)
c               XM3(1)=XCH(N,M)
c               XM3(2)=YCH(N,M)
c               XM3(3)=ZCH(N,M)
c               XM4(1)=XCH(N,M+1)
c               XM4(2)=YCH(N,M+1)
c               XM4(3)=ZCH(N,M+1)
c
c               CRLT = chrleps(j)
c
c               ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
c     %              +ABS(XM2(3)-XM1(3))
c               ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
c     %              +ABS(XM3(3)-XM2(3))
c               ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
c     %              +ABS(XM4(3)-XM3(3))
c               ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
c     %              +ABS(XM1(3)-XM4(3))
c               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)
c               IF(ERRMIN.LE.1.0E-6) THEN
c                  IMR0=1
c               END IF
c
c               DO I=1,NPANEL
c                  XLOC=0.
c                  YLOC=0.
c                  ZLOC=0.
c
c                  DO K=1,3
c                     XLOC=XLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,1,K)
c                     YLOC=YLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,2,K)
c                     ZLOC=ZLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,3,K)
c                  ENDDO
c
c                  IMR=IMR0
c
c                  CALL RPAN(XLOC,YLOC,ZLOC,CRLT,FS,FD,FSX,FSY,
c     *                 FDX,FDY,FDZ,0,IMR)
c
c                  IF(IMR.EQ.2) THEN
c                     DO IXYZ=1,3
c                        XMC(IXYZ)=XCTP(I,IXYZ,1)
c                     ENDDO
c                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
c                  END IF
c
c                  IF(I.EQ.J) THEN
c                     IF(IMR.EQ.1) THEN
c                        FD=TWOPI
c                     END IF
c                     IF(FD.LT.0.) THEN
c                        FD=4.0*PI+FD
c                     END IF
c
c                     IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN
c                        FD=2.0*PI
c                     END IF
c                  ELSE
c                     IF(ABS(FD).GT.6.28) THEN
c                        FD=0.0
c                     END IF
c                  END IF
c
c                  IF(IDUCT .NE. 0 .AND. IDOPT .EQ. 1) THEN
c                     IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
c                        FD = 0.0
c                        FS = 0.0
c                     ENDIF
c                  ENDIF
c
c                  A(I)=FD
c                  TEMP1(I)=FS
c
c               ENDDO
c               CALL WRITE1(IUMA,A,NPANEL)
c               CALL WRITE1(IUMB,TEMP1,NPANEL)
c            ENDDO
c         ENDDO
c
c         DO N=1,ncvx
c            DO M=1,mcvt
c
c               J=INDEXC(N,M)
c
c               IMR0 = 0
c
c               DO K=1,4
c                  XV(K)=XVP(J,K)
c                  YV(K)=YVP(J,K)
c                  SIDE(K)=SID(J,K)
c               ENDDO
c               DO K=1,15
c                  S(K)=SS(J,K)
c               ENDDO
c
c               XM1(1)=XVC(N+1,M+1)
c               XM1(2)=YVC(N+1,M+1)
c               XM1(3)=ZVC(N+1,M+1)
c               XM2(1)=XVC(N+1,M)
c               XM2(2)=YVC(N+1,M)
c               XM2(3)=ZVC(N+1,M)
c               XM3(1)=XVC(N,M)
c               XM3(2)=YVC(N,M)
c               XM3(3)=ZVC(N,M)
c               XM4(1)=XVC(N,M+1)
c               XM4(2)=YVC(N,M+1)
c               XM4(3)=ZVC(N,M+1)
c
c               ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
c     %              +ABS(XM2(3)-XM1(3))
c               ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
c     %              +ABS(XM3(3)-XM2(3))
c               ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
c     %              +ABS(XM4(3)-XM3(3))
c               ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
c     %              +ABS(XM1(3)-XM4(3))
c               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)
c               IF(ERRMIN.LE.1.0E-6) THEN
c                  IMR0=1
c               END IF
c
c               CRLT = chrleps(j)
c
c               DO I=1,NPANEL
c                  XLOC=0.
c                  YLOC=0.
c                  ZLOC=0.
c                  DO K=1,3
c                     XLOC=XLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,1,K)
c                     YLOC=YLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,2,K)
c                     ZLOC=ZLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,3,K)
c                  ENDDO
c
c                  IMR=IMR0
c
c                  CALL RPAN(XLOC,YLOC,ZLOC,CRLT,FS,FD,FSX,FSY,
c     *                 FDX,FDY,FDZ,0,IMR)
c
c                  IF(IMR.EQ.2) THEN
c                     DO IXYZ=1,3
c                        XMC(IXYZ)=XCTP(I,IXYZ,1)
c                     ENDDO
c                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
c                  END IF
c                  IF(I.EQ.J) THEN
c                     IF(IMR.EQ.1) THEN
c                        FD=TWOPI
c                     END IF
c                     IF(FD.LT.0.) THEN
c                        FD=4.0*PI+FD
c                     END IF
c
c                     IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN
c                        FD=2.0*PI
c                     END IF
c                  ELSE
c                     IF(ABS(FD).GT.6.28) THEN
c                        FD=0.0
c                     END IF
c                  END IF
c                  A(I)=FD
c                  TEMP1(I)=FS
c               ENDDO
c
c               IF(IDUCT .NE. 0 .AND. IDOPT .EQ. 1) THEN
c                  IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
c                     FD = 0.0
c                     FS = 0.0
c                  ENDIF
c               ENDIF
c
c               CALL WRITE1(IUMA,A,NPANEL)
c               CALL WRITE1(IUMB,TEMP1,NPANEL)
c
c            ENDDO
c         ENDDO
c      ENDIF
C************************************************************************
C/e S.N.KIM | Aug. 2018.
C************************************************************************

CSH--REPLACE, NBLADE=1-------------------------
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------

      RETURN
      END SUBROUTINE

  	SUBROUTINE RUDIM_PAR(iv1,iv2,iv3,iv4,xv12,xv13,hypot_geo_data)
csingh-------Subroutine added on 08/26/2008 (Sowmitra Singh)------------
csingh-------iv1 = J (panel index) -------------------------------------
csingh-------iv2 = I (panel at which the influence is considered)-------
csingh-------iv3 = KK (Blade index)-------------------------------------
csingh-------iv4 = IMR (determines if Rpan needs to be called or Hypot)-
csingh-------xv5 = XM1 (stores coordinates of panel N,M)----------------
csingh-------xv6 = XM2 (Stores coordinates of panel N,M+1)--------------
csingh-------xv7 = XM3 (Stores coordinates of panel N+1,M+1)------------
csingh-------xv8 = XM4 (Stores coordinates of panel N+1,M)--------------
       
        USE m_hypot_infcof
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC'
        TYPE(t_hypot_infcof)::hypot_geo_data 
        DIMENSION xv5(3),xv6(3),xv7(3),xv8(3),XMC(3)

        DO I = 1, 3
          xv5(I) = DBLE(hypot_geo_data%XM1hypot(I))
          xv6(I) = DBLE(hypot_geo_data%XM2hypot(I))
          xv7(I) = DBLE(hypot_geo_data%XM3hypot(I))
          xv8(I) = DBLE(hypot_geo_data%XM4hypot(I))
        END DO   

        xv12 = ZERO
        xv13 = ZERO
        xv9 = 1
        xv14 = YB(1,1)
csingh------xv14 (lower wall), xv9 (upper wall) ------------------------
        xv18 = xv9 - xv14 
        iv15 = wingimag
csingh-----sets of images that are considered on each wall-------------- 
        Do 100 iv16=1,iv15
csingh--------Considering top wall ----- -------------------------------
            XLOC=ZERO
            YLOC=ZERO 
            ZLOC=ZERO
            Do 10 K=1,3
               IF (K.EQ.2) THEN
                  xv17 = xv9 - XCTP(iv2,K,iv3)
                  IF (mod(iv16,2).GT.0) THEN 
                     xv10 = xv9 + (iv16-1)*xv18 + xv17 
                  ELSE
                     xv10 = xv9 + iv16*xv18 - xv17 
                  ENDIF
csingh-------xv10 = y coordinate of the image --------------------------

                  XLOC = XLOC + (xv10-XCT(iv1,K))*DIR(iv1,1,K)
                  YLOC = YLOC + (xv10-XCT(iv1,K))*DIR(iv1,2,K)
                  ZLOC = ZLOC + (xv10-XCT(iv1,K))*DIR(iv1,3,K)
               ELSE
                  XLOC = XLOC + (XCTP(iv2,K,iv3)-XCT(iv1,K))*
     *                   DIR(iv1,1,K)
                  YLOC = YLOC + (XCTP(iv2,K,iv3)-XCT(iv1,K))*
     *                   DIR(iv1,2,K)
                  ZLOC = ZLOC + (XCTP(iv2,K,iv3)-XCT(iv1,K))*
     *                   DIR(iv1,3,K)
               ENDIF
10          CONTINUE

            CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                 FSY,FDX,FDY,FDZ,0,iv4)

            IF (iv4.EQ.2) THEN
csingh-----HYPOT being called ------------------------------------------
               DO 20 iv11=1,3
                  IF (iv11.EQ.2) THEN
                     xv17 = xv9 - XCTP(iv2,iv11,iv3)
                     IF (mod(iv16,2).GT.0) THEN 
                        XMC(iv11) = xv9 + (iv16-1)*xv18 + xv17 
                     ELSE
                        XMC(iv11) = xv9 + iv16*xv18 - xv17 
                     ENDIF
csingh--------y co-ordinate of the image--------------------------
                  ELSE
                     XMC(iv11)=XCTP(iv2,iv11,iv3)
                  ENDIF
20             CONTINUE 
               CALL HYPOT(xv5,xv6,xv7,xv8,XMC,FD,FS)
            ENDIF

            xv12 =  xv12 + FD
            xv13 =  xv13 + FS
csingh--------Considering bottom wall ----- ----------------------------
            XLOC=ZERO
            YLOC=ZERO 
            ZLOC=ZERO
            Do 30 iv20=1,3
               IF (iv20.EQ.2) THEN
                  xv19 = XCTP(iv2,iv20,iv3) - xv14 
                  IF (mod(iv16,2).GT.0) THEN 
                     xv10 = xv14 - (iv16-1)*xv18 - xv19 
                  ELSE
                     xv10 = xv14 - iv16*xv18 + xv19 
                  ENDIF
csingh-------xv10 = y coordinate of the image --------------------------

                  XLOC = XLOC + (xv10-XCT(iv1,iv20))*DIR(iv1,1,iv20)
                  YLOC = YLOC + (xv10-XCT(iv1,iv20))*DIR(iv1,2,iv20)
                  ZLOC = ZLOC + (xv10-XCT(iv1,iv20))*DIR(iv1,3,iv20)
               ELSE
                  XLOC = XLOC + (XCTP(iv2,iv20,iv3)-XCT(iv1,iv20))*
     *                   DIR(iv1,1,iv20)
                  YLOC = YLOC + (XCTP(iv2,iv20,iv3)-XCT(iv1,iv20))*
     *                   DIR(iv1,2,iv20)
                  ZLOC = ZLOC + (XCTP(iv2,iv20,iv3)-XCT(iv1,iv20))*
     *                   DIR(iv1,3,iv20)
               ENDIF
30         CONTINUE

           CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                FSY,FDX,FDY,FDZ,0,iv4)

           IF (iv4.EQ.2) THEN
csingh-----HYPOT being called ------------------------------------------
              DO 40 iv21=1,3
                 IF (iv21.EQ.2) THEN
                    xv19 = XCTP(iv2,iv21,iv3) - xv14 
                    IF (mod(iv16,2).GT.0) THEN 
                       XMC(iv21) = xv14 - (iv16-1)*xv18 - xv19
                    ELSE
                       XMC(iv21) = xv14 - iv16*xv18 + xv19
                    ENDIF
csingh--------y co-ordinate of the image-------------------------------
                 ELSE
                    XMC(iv21)=XCTP(iv2,iv21,iv3)
                 ENDIF
40            CONTINUE 

              CALL HYPOT(xv5,xv6,xv7,xv8,XMC,FD,FS)
           ENDIF

csingh------------------------------------------------------------------

           xv12 =  xv12 + FD
           xv13 =  xv13 + FS

100     CONTINUE

        RETURN
        END


      SUBROUTINE RPAN_PAR3(Xi,Yi,Zi,CHRLENSi, FSs,FDd,FSXi,FSYi,
     %                  FDXi,FDYi,FDZi,ID,IFLAG,rpan_geo_data)
      use m_rpan_infcof
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      type(t_rpan_infcof)::rpan_geo_data
      REAL XI,YI,ZI,CHRLENSI,FSS,FDD,FSXI,FSYI,FDXI,FDYI,FDZI
      DIMENSION XV(4),YV(4),S(15),SIDE(4)
      DIMENSION R(4),RR(4),RI(4),XRI(4),YRI(4),FE(4),
     *          B(5),XMXV(4),YMYV(4),N1(4)

C----------------------------------------------------------------------
C
C  DATA ARRAY B FOR 6D RATIONAL APPROXIMATION OF ARCTANGENT
C  SEE HART, ET AL, COMPUTER APPROXIMATIONS, WILEY, 1968, TABLE 5090
C
C----------------------------------------------------------------------
      DATA B/ 2.4091197D-01, 3.7851122D+00, 5.6770721D+00,
     *        5.6772854D+00, 5.6770747D+00/
C      DATA PI/ 3.1415927D+00 /, PI2/ 1.5707963D+00 /
C      DATA TWOPI/ 6.2831853D+00 /
      DATA ONE/ 1.D+00/, A3/ 3.D+00/, A5/ 5.D+00/
      DATA A7/ 7.D+00/, A9/ 9.D+00/, A11/ 11.D+00/
      DATA A14/ 14.D+00/, A35/ 35.D+00/
      DATA A49/ 49.D+00/, A63/ 63.D+00/, A99/ 99.D+00/
C      DATA ONE10/.1D+00/, ONE6/.1666667D+00/, ONE3/ .3333333D+00/
C      DATA FIVE3/ 1.666667D+00/, SEVEN3/2.333333D+00/, ZERO/ 0.0D0/
      DATA ZERO/ 0.0D0/
      DATA FT3/ 4.666667D+00/, TOL/ 1D-08 /
      DATA N1/ 2, 3, 4, 1 /

      pi = dacos(-1.d0)
      pi2 = 0.5d0*pi
      twopi = 2.d0*pi
      one10 = 0.1d0
      one6 = 1.d0/6.d0
      one3 = 1.d0/3.d0
      five3 = 5.d0/3.d0
      seven3 = 7.d0/3.d0


C --- Change Input variable to double precision

      x = dble(xi)
      y = dble(yi)
      z = dble(zi)
      chrlens = dble(chrlensi)
      do i = 1 , 4
      xv(i) = dble(rpan_geo_data%XVrpan(i))
      yv(i) = dble(rpan_geo_data%YVrpan(i))
      side(i) = dble(rpan_geo_data%SIDErpan(i))
      enddo
      do i = 1 , 15
        s(i) = dble(rpan_geo_data%SQrpan(i))
      enddo

      XMXC=X-S(6)
      YMYC=Y-S(2)
      XX=XMXC*XMXC
      YY=YMYC*YMYC
      ZZ=Z*Z
      RRC=XX+YY+ZZ
      IF (RRC.LT.100.d0*CHRLENS) THEN
         IF(IFLAG.EQ.1) THEN
            GO TO 11
         ELSE IF (IFLAG.EQ.0) THEN
            IFLAG=2
            RETURN
         END IF
      END IF
C----------------------------------------------------------------------
C
C   TWO-TERM MULTIPOLE EXPANSION INCLUDING SECOND MOMENTS
C
C----------------------------------------------------------------------
      R2=ONE/RRC
      R1=DSQRT(R2)
      R3=R1*R2
      R5=R3*R2
      ZR2=Z*R2
      XY=XMXC*YMYC
      SS1=S(1)*R1
      SS3=-(S(3)+S(10))*R3
      SS5=(XX*S(10)+XY*S(7)+YY*S(3))*R5
      FS=SS1+ONE3*SS3+SS5
      FDSUM=SS1+SS3+A5*SS5
      FD=ZR2*FDSUM
C----------------------------------------------------------------------
C    DELETE NEXT  14  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 8
      RSS3=R2*SS1
      SSX3=-XMXC*RSS3
      SSY3=-YMYC*RSS3
      SSX5=(XMXC*(S(3)+A3*S(10))+YMYC*S(7))*R5
      SSY5=(YMYC*(S(10)+A3*S(3))+XMXC*S(7))*R5
      A5R2=A5*R2
      RSS7=-A5R2*SS5
      SSX7=XMXC*RSS7
      SSY7=YMYC*RSS7
      FSX=SSX3+SSX5+SSX7
      FSY=SSY3+SSY5+SSY7
      FDX=ZR2*(A3*SSX3+A5*SSX5+A7*SSX7)
      FDY=ZR2*(A3*SSY3+A5*SSY5+A7*SSY7)
      ZZR4=ZR2*ZR2
      FDZ=R2*FDSUM-ZZR4*(A3*SS1+A5*SS3+A35*SS5)
  8   IF (RRC.GT.15.d0*CHRLENS) GO TO 99
C----------------------------------------------------------------------
C
C    THIRD AND FOURTH MOMENTS ADDED FOR RRC/AREA BETWEEN 40 AND 150
C
C----------------------------------------------------------------------
      S914=S(9)+S(14)
      S813=S(8)+S(13)
      S411=S(4)+S(11)
      S512=S(5)+S(12)
      S1215=S(12)+S(15)
      R7=R5*R2
      R9=R7*R2
      SS5=(-XMXC*S813-YMYC*S411+ONE10*(S512+S1215))*R5
      SS7=(FIVE3*((XMXC*XX*S(13)+YMYC*YY*S(4))+A3*XY*(XMXC*S(11)
     *+YMYC*S(8)))-XX*S1215-YY*S512-XY*S914)*R7
      SS9=(A7*(ONE6*(XX*XX*S(15)+YY*YY*S(5))+XX*YY*S(12))
     *+SEVEN3*XY*(XX*S(14)+YY*S(9)))*R9
      FS=FS+SS5+SS7+SS9
      FDSUM=A5*SS5+A7*SS7+A9*SS9
      FD=FD+ZR2*FDSUM
C----------------------------------------------------------------------
C    DELETE NEXT  20  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 99
      TXY=XY+XY
      SSX5=-S813*R5
      SSY5=-S411*R5
      RSS7=A5R2*SS5
      SSX7=(A5*(XX*S(13)+TXY*S(11)+YY*S(8))-S1215*(XMXC+XMXC)
     *-YMYC*S914)*R7-XMXC*RSS7
      SSY7=(A5*(YY*S(4)+XX*S(11)+TXY*S(8))-S512*(YMYC+YMYC)
     *-XMXC*S914)*R7-YMYC*RSS7
      RSS9=A7*SS7*R2
      SSX9=(FT3*XMXC*XX*S(15)+A14*XMXC*YY*S(12)+A49*YMYC*(XX*S(14)
     *+ONE3*YY*S(9)))*R9-XMXC*RSS9
      SSY9=(FT3*YMYC*YY*S(5)+A14*YMYC*XX*S(12)+A49*XMXC*(YY*S(9)
     *+ONE3*XX*S(14)))*R9-YMYC*RSS9
      RSS11=A9*SS9*R2
      SSX11=-XMXC*RSS11
      SSY11=-YMYC*RSS11
      FSX=FSX+SSX5+SSX7+SSX9+SSX11
      FSY=FSY+SSY5+SSY7+SSY9+SSY11
      FDX=FDX+ZR2*(A5*SSX5+A7*SSX7+A9*SSX9+A11*SSX11)
      FDY=FDY+ZR2*(A5*SSY5+A7*SSY7+A9*SSY9+A11*SSY11)
      FDZ=FDZ+R2*FDSUM-ZZR4*(A35*SS5+A63*SS7+A99*SS9)
      GO TO 99
C----------------------------------------------------------------------
C
C    NEAR-FIELD SECTION USES EXACT FORMULATION
C      SET Z=TOL IF Z.LT.TOL TO AVOID INDETERMINACY ON PANEL
C      ZVTX IS USED TO DETERMINE PROXIMITY TO VERTEX NORMALS
C      MFLAG=1 IF NEAR VERTEX NORMALS
C
C----------------------------------------------------------------------
  11  FD=ZERO
      FS=ZERO
      ABZ=DABS(Z)
      IF(ABZ.GT.TOL) GO TO 12
        Z=TOL
        ZZ=Z*Z
        ABZ=TOL
  12  ZVTX=1.005D0*ABZ
      MFLAG=0
C----------------------------------------------------------------------
C    DELETE NEXT 5 LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      FSX=ZERO
      FSY=ZERO
      FDX=ZERO
      FDY=ZERO
      FDZ=ZERO
C----------------------------------------------------------------------
C
C    LOOP FOR CORNER FUNCTIONS
C
C----------------------------------------------------------------------
      DO 13 N=1,4
        XMXV(N)=X-XV(N)
        YMYV(N)=Y-YV(N)
        XX=XMXV(N)*XMXV(N)
        YY=YMYV(N)*YMYV(N)
        FE(N)=ZZ+XX
        RR(N)=FE(N)+YY
        R(N)=DSQRT(RR(N))
        IF (R(N).LT.ZVTX) MFLAG=1
        RI(N)=ONE/R(N)
C----------------------------------------------------------------------
C    DELETE NEXT   3   LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 13
        XRI(N)=RI(N)*XMXV(N)
        YRI(N)=RI(N)*YMYV(N)
   13 CONTINUE
C----------------------------------------------------------------------
C
C    LOOP FOR SIDE FUNCTIONS AND SUMS OVER FOUR SIDES
C
C----------------------------------------------------------------------
      DO 33 N=1,4
        IF (SIDE(N).LT.TOL) GO TO 33
        SIDI=ONE/SIDE(N)
        CT=(XV(N1(N))-XV(N))*SIDI
        ST=(YV(N1(N))-YV(N))*SIDI
        V=XMXV(N)*ST-YMYV(N)*CT
        VV=V*V
        RADS=VV+ZZ
        U1=XMXV(N)*CT+YMYV(N)*ST
        U2=XMXV(N1(N))*CT+YMYV(N1(N))*ST
        RSUM=R(N)+R(N1(N))
        FLAG=RI(N)*RI(N1(N))*U1*U2
C----------------------------------------------------------------------
C
C       FLAG=1 ON EXTENSIONS, -1 ON SIDES
C         IN FOLLOWING SUBSECTIONS FS,FSX,FSY,FDZ ARE EVALUATED FROM
C           LAST FORM OF (3.9) IN NORMAL CASE
C         ELSE
C           FIRST FORM OF (3.9) IF NEAR SIDE OF PANEL
C
C----------------------------------------------------------------------
        IF (FLAG.GT.-.99D0) THEN
          RSP=RSUM+SIDE(N)
          RSM=RSUM-SIDE(N)
          FLN=DLOG(RSP/RSM)
        ELSE
          RU1=R(N)+U1
          RU2=R(N1(N))-U2
          RADI=ONE/RADS
          FLN=DLOG(RU1*RU2*RADI)
        ENDIF
        FS=FS+V*FLN
C----------------------------------------------------------------------
C    DELETE FOLLOWING LINES TO NEXT ENDIF TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 14
        IF (FLAG.GT.-.99D0) THEN
          FAC=V*(SIDE(N)+SIDE(N))/(RSP*RSM)
          FSX=FSX+FLN*ST-FAC*(XRI(N)+XRI(N1(N)))
          FSY=FSY-FLN*CT-FAC*(YRI(N)+YRI(N1(N)))
          FDZ=FDZ-FAC*(RI(N)+RI(N1(N)))
        ELSE
          RU1I=ONE/RU1
          RU2I=ONE/RU2
          FA=RU1I-RU2I
          FB=-(V+V)*RADI
          FSX=FSX+FLN*ST+V*(FA*CT+FB*ST+RU1I*XRI(N)+RU2I*XRI(N1(N)))
          FSY=FSY-FLN*CT+V*(FA*ST-FB*CT+RU1I*YRI(N)+RU2I*YRI(N1(N)))
          FDZ=FDZ+FB+V*(RU1I*RI(N)+RU2I*RI(N1(N)))
        ENDIF
C----------------------------------------------------------------------
C
C        IN FOLLOWING SUBSECTIONS FACTORS IN (2.15) ARE EVALUATED FROM
C          (2.7) IN NORMAL CASE
C        ELSE
C          (2.14) IF NEAR NORMAL TO A VERTEX
C
C----------------------------------------------------------------------
   14   IF (MFLAG.EQ.0) THEN
          S1=V*R(N)
          C1=ABZ*U1
          S2=V*R(N1(N))
          C2=ABZ*U2
        ELSE
          FH1=XMXV(N)*YMYV(N)
          FH2=XMXV(N1(N))*YMYV(N1(N))
          S1=FE(N)*ST-FH1*CT
          C1=ABZ*R(N)*CT
          S2=FE(N1(N))*ST-FH2*CT
          C2=ABZ*R(N1(N))*CT
        ENDIF
        S12=S1*C2-S2*C1
        C12=C1*C2+S1*S2
C----------------------------------------------------------------------
C
C    EVALUATE THIRD ARCTANGENT IN (2.15)
C         ANGLE (MODULO PI) BETWEEN -PI/4 AND PI/4
C       ELSE
C         USE INVERSE COTANGENT AND ADD/SUBTRACT PI/2
C
C----------------------------------------------------------------------
        IF (DABS(S12).LE.DABS(C12)) THEN
          U=S12/C12
          IF (C12.LT.ZERO) FD=FD+DSIGN(PI,S12)
        ELSE
          U=-C12/S12
          FD=FD+DSIGN(PI2,S12)
        ENDIF
        UU=U*U
        FD=FD+U*((B(1)*UU+B(2))*UU+B(3))/((UU+B(4))*UU+B(5))
C----------------------------------------------------------------------
C
C    FOLLOWING THREE SECTIONS EVALUATE FDX,FDY FOR
C       FIELD POINT NEAR NORMAL TO VERTEX (MFLAG=1)
C       NORMAL CASE (FLAG.LT.0.99)
C       NEAR EXTENSIONS OF SIDES (FLAG.GE.0.99)
C
C    DELETE ALL FOLLOWING LINES ABOVE LABEL 33 TO OMIT DERIVATIVES
C----------------------------------------------------------------------
        IF (ID.EQ.0) GO TO 33
        IF (MFLAG.EQ.0) GO TO 20
          FAC=C1/((C1*C1+S1*S1)*RR(N))
        IF(Z.LT.ZERO) FAC=-FAC   ! Fix sign of derivatives, Newman
          FDX=FDX+(RR(N)*V+FH1*U1)*FAC
          FDY=FDY-FE(N)*U1*FAC
          FAC=C2/((C2*C2+S2*S2)*RR(N1(N)))
        IF(Z.LT.ZERO) FAC=-FAC   ! Fix sign of derivatives, Newman
          FDX=FDX-(RR(N1(N))*V+FH2*U2)*FAC
          FDY=FDY+FE(N1(N))*U2*FAC
          GO TO 33
  20    IF (FLAG.LT.0.99D0) THEN
          U1V=U1*V
          FAC=Z/(C1*C1+S1*S1)
          FDX=FDX+(U1V*XRI(N)+R(N)*YMYV(N))*FAC
          FDY=FDY+(U1V*YRI(N)-R(N)*XMXV(N))*FAC
          U2V=U2*V
          FAC=Z/(C2*C2+S2*S2)
          FDX=FDX-(U2V*XRI(N1(N))+R(N1(N))*YMYV(N1(N)))*FAC
          FDY=FDY-(U2V*YRI(N1(N))-R(N1(N))*XMXV(N1(N)))*FAC
        ELSE
          ZS=Z*SIDE(N)
          USUM=U1+U2
          VRADS=V*RADS
          SFAC=VRADS*USUM
          SFS=-SFAC*ZS
          SFA=SFAC*C12
          CFAC=U2*R(N)+U1*R(N1(N))
          SFB=SFAC*CFAC
          CCF=C12*CFAC
          PA=(CCF+CCF)*VRADS-SFA*RSUM-SFB*ZZ*USUM
          PB=CCF*USUM*(VV+VV+RADS)-SFB*(S1+S1)*R(N1(N))
          PC=-SFA*U2-SFB*VV*R(N1(N))
          PD=-SFA*U1-SFB*VV*R(N)
          FAC=ZS/(CCF*CCF+SFS*SFS)
          FDX=FDX-(PA*CT+PB*ST+PC*XRI(N)+PD*XRI(N1(N)))*FAC
          FDY=FDY-(PA*ST-PB*CT+PC*YRI(N)+PD*YRI(N1(N)))*FAC
        ENDIF
  33  CONTINUE
      IF (FD.LT.ZERO) FD=FD+TWOPI
      IF (Z.LT.ZERO) FD=-FD
      FS=FS-Z*FD
C----------------------------------------------------------------------
C    DELETE NEXT  3  LINES TO OMIT DERIVATIVES
C----------------------------------------------------------------------
      IF (ID.EQ.0) GO TO 99
      FSX=FSX-Z*FDX
      FSY=FSY-Z*FDY

99    CONTINUE

      fdd = sngl(fd)
      fss = sngl(fs)
      if(id .ne. 0) then
        fsxi = sngl(fsx)
        fsyi = sngl(fsy)
        fdxi = sngl(fdx)
        fdyi = sngl(fdy)
        fdzi = sngl(fdz)
      endif

      RETURN
      END


C**************************** MIT PSF10 ********************************
C
C                PROPELLER STEADY FLOW ANALYSIS PROGRAM
C                        PANEL METHOD SOLUTION         
C
C                             Version 1.0
C
C       Copyright (c) Massachusetts Institute of Technology 1997
C
C                     Release Date: 1 January 1997
C
C***********************************************************************
      SUBROUTINE HYPOT_PAR2(FDR,FSR,hypot_geo_data)
C**********************************************************************
C     DOUBLE PRECISION
C     Compute the potential due to sources and dipoles based on 
C      Morino's formula
C     -- x1, x2, x3 and x4 should be input in the clock wise direction
C     -- This subroutine CAN NOT calculate the influence functions of 
C        a TRIANGULAR panel
C     
C     10-04-88 C.Y.HSIN @MHL
C     14-08-17 S.N.KIM  @MHL (UT-AUSTIN)
C
      USE m_hypot_infcof
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TYPE(t_hypot_infcof)::hypot_geo_data
      REAL FDR,FSR 
      DIMENSION X1(3),X2(3),X3(3),X4(3),P(3),PC(3),P1(3),P2(3),P3(3),
     *          XCC(3),RC(3),XI0(4),ETA0(4),FDP(4),FSP(4),  
     *          R(3),A1(3),A2(3),U(3),XNC(3),XN(3),RXA1(3),RXA2(3),
     *          NNEG(4)
      DATA XI0 /-1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA ETA0/-1.0D0,1.0D0,1.0D0,-1.0D0/
      DATA ZERO,QUAD,HALF,TWO,FOUR/0.0D0,0.25D0,0.50D0,
     *                                 2.0D0,4.0D0/
C
C.....Transfer to double precision
      DO 10 I=1,3
         X1(I)=DBLE(hypot_geo_data%XM1hypot(I))         
         X2(I)=DBLE(hypot_geo_data%XM2hypot(I))         
         X3(I)=DBLE(hypot_geo_data%XM3hypot(I))         
         X4(I)=DBLE(hypot_geo_data%XM4hypot(I))         
         XCC(I)=DBLE(hypot_geo_data%XMChypot(I))         
 10   CONTINUE
C      PI=3.14159265D0
      PI = DACOS(-1.D0)
      TWOPI=TWO*PI 
      FOURPI=FOUR*PI 
      PI2=HALF*PI
      DO 20 I=1,3
         PC(I)=QUAD*( X1(I)+X2(I)+X3(I)+X4(I) )
         P1(I)=QUAD*( X3(I)+X4(I)-X2(I)-X1(I) )           
         P2(I)=QUAD*( X3(I)-X4(I)+X2(I)-X1(I) )           
         P3(I)=QUAD*( X3(I)-X4(I)-X2(I)+X1(I) )           
 20   CONTINUE  
C
C.....normal vector at the center point
C
      XI=ZERO
      ETA=ZERO
      DO 30 I=1,3
         P(I)=PC(I)+XI*P1(I)+ETA*P2(I)+XI*ETA*P3(I)
         RC(I)=P(I)-XCC(I)            
         A1(I)=P1(I)+ETA*P3(I)
         A2(I)=P2(I)+XI*P3(I)
 30   CONTINUE
      CALL EXPROD(A1,A2,U)
      UL=DSQRT( U(1)*U(1)+U(2)*U(2)+U(3)*U(3) )
      RCL=DSQRT( RC(1)*RC(1)+RC(2)*RC(2)+RC(3)*RC(3) )
      DO 40 I=1,3 
         XNC(I)=U(I)/UL
 40   CONTINUE        
C
C.....Compute Is, Id at four corner points       
      DO 50 J=1,4
         XI=XI0(J)           
         ETA=ETA0(J)
         DO 60 I=1,3
            P(I)=PC(I)+XI*P1(I)+ETA*P2(I)+XI*ETA*P3(I)
            R(I)=P(I)-XCC(I)
            A1(I)=P1(I)+ETA*P3(I)
            A2(I)=P2(I)+XI*P3(I)
 60      CONTINUE
         CALL EXPROD(A1,A2,U)
         UL=DSQRT( U(1)*U(1)+U(2)*U(2)+U(3)*U(3) )
         DO 70 I=1,3 
            XN(I)=U(I)/UL
 70      CONTINUE
         CALL EXPROD(R,A1,RXA1)
         CALL EXPROD(R,A2,RXA2)
         CALL ENPROD(RXA1,RXA2,RA1RA2)
         CALL ENPROD(R,U,RU)
         CALL ENPROD(R,A1,RA1)
         CALL ENPROD(R,A2,RA2)
         CALL ENPROD(RXA1,XNC,RXA1NC)
         CALL ENPROD(RXA2,XNC,RXA2NC)
         CALL ENPROD(RC,XNC,RNC)
         RL=DSQRT( R(1)*R(1)+R(2)*R(2)+R(3)*R(3) )
         RXA1L=DSQRT( RXA1(1)*RXA1(1)+RXA1(2)*RXA1(2)+RXA1(3)*RXA1(3) )        
         RXA2L=DSQRT( RXA2(1)*RXA2(1)+RXA2(2)*RXA2(2)+RXA2(3)*RXA2(3) )        
         A1L=DSQRT( A1(1)*A1(1)+A1(2)*A1(2)+A1(3)*A1(3) )        
         A2L=DSQRT( A2(1)*A2(1)+A2(2)*A2(2)+A2(3)*A2(3) )        
         DEN=RL*RU
         FDP(J)=ATAN3( DEN,RA1RA2 ) 
C........Determine the quardrants
         IF(FDP(J).GE.ZERO) THEN
            NNEG(J)=0
         ELSE
            NNEG(J)=1
         END IF
         FS1=RNC*FDP(J)
         FS2=RXA1NC*DASINH(RA1/RXA1L)/A1L
         FS3=RXA2NC*DASINH(RA2/RXA2L)/A2L
         FSP(J)=-FS2+FS3
 50   CONTINUE
C
C     Compute potential of source and dipole
      FD=FDP(1)-FDP(2)+FDP(3)-FDP(4)
      MNEG=NNEG(1)+NNEG(2)+NNEG(3)+NNEG(4)
      IF(MNEG.EQ.1) THEN
         IF(NNEG(1).EQ.1.OR.NNEG(3).EQ.1) THEN
            FD=FD+TWOPI 
         END IF
      ELSE IF(MNEG.EQ.3) THEN
         IF(NNEG(1).EQ.0.OR.NNEG(3).EQ.0) THEN
            FD=FD-TWOPI 
         END IF
      END IF
      IF(FD.LT.-6.28318531) THEN
         FD=FOURPI+FD
      END IF
      FS1=FSP(1)-FSP(2)+FSP(3)-FSP(4)
      FS=-FD*RNC+FSP(1)-FSP(2)+FSP(3)-FSP(4)
      FDR=SNGL(FD)
      FSR=SNGL(FS)
      RETURN
C))))))))))))))))))))))) End of subroutine HYPOT  (((((((((((((((((((((
      END

      SUBROUTINE INFCOF_DEAD
************************************************************************
C/s S.N.KIM | This subroutine is the old version of 'INFCOF' used before
C           | PROPCAV V3.3, which was not OpenMP parallelized.
C/e S.N.KIM | Aug. 2018.
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
c XM YU 10/2011
      INCLUDE 'PUFCAVC.INC'
c XM YU 10/2011

      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)

C-----------------------------------------------------------------------
C     Initialization
C-----------------------------------------------------------------------

C.....FILE 41 (IUMA)  -- dipole inf. functions for each blade
C.....FILE 42 (IUMB)  -- source inf. functions for each blade
c      write(*,*) 'beta infcof'

      IUMA=41
      IUMB=42

      IDDUCT2 = NPANB + NPANH
      IDDUCT3 = NPANB + NPANH + NPAND

      IAAA1 = NPANB + NPANH + NPAND
      IAAA2 = NPANB + NPANH + NPAND + NPANTN

      sum1 = 0.0
      sum2 = 0.0

C      write(*,*) 'npanb = ',npanb
C      write(*,*) 'npanh = ',npanh
C      write(*,*) 'npand = ',npand
C      write(*,*) 'npantn = ',npantn

      CALL CLEAR(A,NPANZ)
      CALL CLEAR(TEMP1,NPANZ)

      IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2

      DO M=MR,1,-1
         DO N=1,NC
            J=INDEXB(N,M)
            DO K=1,4
               XV(K)=XVP(J,K)
               YV(K)=YVP(J,K)
               SIDE(K)=SID(J,K)
            END DO
            DO K=1,15
               S(K)=SS(J,K)
            END DO

            XM1(1)=XB(N,M)
            XM1(2)=YB(N,M)
            XM1(3)=ZB(N,M)
            XM2(1)=XB(N,M+1)
            XM2(2)=YB(N,M+1)
            XM2(3)=ZB(N,M+1)
            XM3(1)=XB(N+1,M+1)
            XM3(2)=YB(N+1,M+1)
            XM3(3)=ZB(N+1,M+1)
            XM4(1)=XB(N+1,M)
            XM4(2)=YB(N+1,M)
            XM4(3)=ZB(N+1,M)
            IMR0=0

            ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))+ABS(XM2(3)-XM1(3))
            ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))+ABS(XM3(3)-XM2(3))
            ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))+ABS(XM4(3)-XM3(3))
            ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))+ABS(XM1(3)-XM4(3))
            ERRMIN=AMIN1(ER1,ER2,ER3,ER4)

            IF(ERRMIN.LE.1.0E-6) THEN
               IMR0=1
            END IF

            DO KK=1,NBLADE
               DO I=1,NPANEL
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO K=1,3
                     XLOC=XLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                  END DO

                  IMR=IMR0
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),
     *                 FS,FD,FSX,FSY,FDX,FDY,FDZ,0,IMR)

                  IF(IMR.EQ.2) THEN
                     DO IXYZ=1,3
                        XMC(IXYZ)=XCTP(I,IXYZ,KK)
                     END DO
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF

                  IF(I.EQ.J .AND. KK.EQ.1) THEN
                     IF(IMR.EQ.1) THEN
                        FD=TWOPI
                     END IF
                     IF(FD.LT.0.) THEN
                        FD=4.0*PI+FD
                     END IF

                     IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN
C                        WRITE(*,*) 'INFCOF-B',NTSTEP,KK,I,M,N,IMR,FD
                        FD=2.0*PI
                     END IF

                  ELSE
                     IF(ABS(FD).GT.6.28) THEN
C                        WRITE(*,*) 'INFCOF-D',NTSTEP,KK,I,M,N,IMR,FD
                        FD=0.0
                     END IF

                  END IF
c IMAGE MODEL
                  if ((icon.eq.5).and.(wingimag.ne.0)) then

                     call rudim(j,i,kk,imr,xm1,xm2,xm3,xm4,xv12,xv13)
                     fd=fd+xv12
                     fs=fs+xv13
                  endif
c IMAGE MODEL

                  IF(ITUN .NE. 0) THEN
                     IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                      FD = -FD
                      FS = -FS
                     ENDIF
                  ENDIF

                  IF(IDUCT .EQ. 1) THEN
                     IF(IDOPT .EQ. 1) THEN
                        IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                           FD = 0.0
                           FS = 0.0
                        ENDIF
                     ENDIF
                  ENDIF

                  A(I)=FD
                  TEMP1(I)=FS

                if (IVISC.eq.1) then
c XM YU 10/2011 ------------------------------------------------------------
c for viscous calculation
                if(kk.eq.1) then
                   kt=nc*(mr-m)+1
                   if(i.ge.kt.and.i.lt.kt+nc)then
                      if(j.ge.kt.and.j.lt.kt+nc)then
                         ii=i+1-kt
                         jj=j+1-kt
                         AVL(ii,jj,m)=A(I)
ccc..... add wake influence
                         if(n.eq.1)then
                              AVL(ii,jj,m)=avl(ii,jj,m)-w(i,m)
                         else if(n.eq.nc)then
                              AVL(ii,jj,m)=avl(ii,jj,m)+w(i,m)
                         endif
                      endif
                   endif
                endif
                end if
c XM YU 10/2011-----------------------------------------------------------


               END DO
               CALL WRITE1(IUMA,A,NPANEL)
               CALL WRITE1(IUMB,TEMP1,NPANEL)
            END DO
         END DO
      END DO
C-----------------------------------------------------------------------
C     Compute influence coefficients due to the hub
C-----------------------------------------------------------------------

      IF(IHUB.NE.0) THEN

        DO N=1,NHBX
          DO M=1,MHBT
            J=INDEXH(N,M)

            DO K=1,4
              XV(K)=XVP(J,K)
              YV(K)=YVP(J,K)
              SIDE(K)=SID(J,K)
            END DO
            DO K=1,15
              S(K)=SS(J,K)
            END DO

            XM1(1)=XH(N,M+1)
            XM1(2)=YH(N,M+1)
            XM1(3)=ZH(N,M+1)
            XM2(1)=XH(N,M)
            XM2(2)=YH(N,M)
            XM2(3)=ZH(N,M)
            XM3(1)=XH(N+1,M)
            XM3(2)=YH(N+1,M)
            XM3(3)=ZH(N+1,M)
            XM4(1)=XH(N+1,M+1)
            XM4(2)=YH(N+1,M+1)
            XM4(3)=ZH(N+1,M+1)
            IMR0=0
            ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
     *           +ABS(XM2(3)-XM1(3))
            ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
     *           +ABS(XM3(3)-XM2(3))
            ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
     *           +ABS(XM4(3)-XM3(3))
            ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
     *           +ABS(XM1(3)-XM4(3))
            ERRMIN=AMIN1(ER1,ER2,ER3,ER4)

            IF(ERRMIN.LE.1.0E-6) THEN
               WRITE(*,2000) M,N
 2000          FORMAT(' ...... triangular panels (m,n) -- > ',
     *              2I4,' (hub) ')
               IMR0=1
            END IF

            DO KK=1,NBLADE
               DO I=1,NPANEL
                XLOC=0.
                YLOC=0.
                ZLOC=0.
                DO K=1,3
                  XLOC=XLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                  YLOC=YLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                  ZLOC=ZLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                END DO

                IMR=IMR0

                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR)

                IF(IMR.EQ.2) THEN
                   DO IXYZ=1,3
                      XMC(IXYZ)=XCTP(I,IXYZ,KK)
                   END DO
                   CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                END IF

                IF(I.EQ.J .AND. KK.EQ.1) THEN
                   IF(IMR.EQ.1) THEN
                      FD=TWOPI
                   END IF
                   IF(FD.LT.0.) THEN
                      FD=4.0*PI+FD
                   END IF

                   IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.5.0) THEN
C                      WRITE(*,*) 'INFCOF-E',NTSTEP,KK,I,M,N,IMR,FD
                      FD=2.0*PI
                   END IF
                ELSE
                   IF(ABS(FD).GT.6.28) THEN
C                      WRITE(*,*) 'INFCOF-F',NTSTEP,KK,I,M,N,IMR,FD
                      FD=0.0
                   END IF

                END IF

                IF(ITUN .NE. 0) THEN
                   IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                      FD = -FD
                      FS = -FS
                   ENDIF
                ENDIF

                IF(IDUCT .EQ. 1) THEN
                   IF(IDOPT .EQ. 1) THEN
                      IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                         FD = 0.0
                         FS = 0.0
                      ENDIF

                   ENDIF
                ENDIF

                A(I)=FD
                TEMP1(I)=FS

              END DO
              CALL WRITE1(IUMA,A,NPANEL)
              CALL WRITE1(IUMB,TEMP1,NPANEL)
            END DO
          END DO
        END DO

      END IF



C-----------------------------------------------------------------------
C     Compute influence coefficients due to DUCT
C-----------------------------------------------------------------------

      IF(IDUCT .NE. 0) THEN

         DO M=1,MDUCT
            DO N=1,NDUCT
               J=INDEXD(N,M)

               DO K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
               ENDDO

               DO K=1,15
                  S(K)=SS(J,K)
               ENDDO

               XM1(1)=XD(N,M+1)
               XM1(2)=YD(N,M+1)
               XM1(3)=ZD(N,M+1)
               XM2(1)=XD(N,M)
               XM2(2)=YD(N,M)
               XM2(3)=ZD(N,M)
               XM3(1)=XD(N+1,M)
               XM3(2)=YD(N+1,M)
               XM3(3)=ZD(N+1,M)
               XM4(1)=XD(N+1,M+1)
               XM4(2)=YD(N+1,M+1)
               XM4(3)=ZD(N+1,M+1)

               IMR0=0
               ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
     *              +ABS(XM2(3)-XM1(3))
               ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
     *              +ABS(XM3(3)-XM2(3))
               ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
     *              +ABS(XM4(3)-XM3(3))
               ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
     *              +ABS(XM1(3)-XM4(3))
               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)

               IF(ERRMIN.LE.1.0E-6) THEN
                  WRITE(*,2100) M,N
 2100             FORMAT(' ...... triangular panels (m,n) -- > ',
     *                 2I4,' (DUCT) ')
                  IMR0=1
               END IF

               DO KK=1,NBLADE
                  DO I=1,NPANEL
                     XLOC=0.
                     YLOC=0.
                     ZLOC=0.
                     DO K=1,3
                        XLOC=XLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                        YLOC=YLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                        ZLOC=ZLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                     ENDDO
                     IMR=IMR0

                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR)

                     IF(IMR.EQ.2) THEN
                        DO IXYZ=1,3
                           XMC(IXYZ)=XCTP(I,IXYZ,KK)
                        ENDDO
                        CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                     END IF

                     IF(I.EQ.J .AND. KK.EQ.1) THEN
                        IF(IMR.EQ.1) THEN
                           FD=TWOPI
                        END IF
                        IF(FD.LT.0.) THEN
                           FD=4.0*PI+FD
                        END IF
                        IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.5.0) THEN
C                           WRITE(*,*) 'INFCOF-E',NTSTEP,KK,I,M,N,IMR,FD
                           FD=2.0*PI
                        END IF
                     ELSE

                        IF(ABS(FD).GT.6.28) THEN
C                           WRITE(*,*) 'INFCOF-F',NTSTEP,KK,I,M,N,IMR,FD
                           FD=0.0
                        END IF

                     END IF

                     IF(ITUN .NE. 0) THEN
                        IF(I .GT. IAAA1 .AND. I .LE. IAAA2) THEN
                          FD = -FD
                          FS = -FS
                        ENDIF
                    ENDIF

                    IF(IDOPT .EQ. 1) THEN
                       IF(I .LE. IDDUCT2 .OR. I .GT. IDDUCT3) THEN
                          FD = 0.0
                          FS = 0.0
                       ENDIF

                    ENDIF

                    A(I)=FD
                    TEMP1(I) = FS

                  ENDDO
                  CALL WRITE1(IUMA,A,NPANEL)
                  CALL WRITE1(IUMB,TEMP1,NPANEL)
               ENDDO

            ENDDO
         ENDDO
      END IF

C-----------------------------------------------------------------------
C     Compute influence coefficients due to the tunnel walls
C-----------------------------------------------------------------------

      IF(ITUN.NE.0) THEN

         DO N=1,NAXT
            DO M=1,MTUNEL
               J=INDEXTN(N,M)

               DO K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
               ENDDO

               DO K=1,15
                  S(K)=SS(J,K)
               ENDDO

               XM1(1)=XTUN(N,M)
               XM1(2)=YTUN(N,M)
               XM1(3)=ZTUN(N,M)
               XM2(1)=XTUN(N,M+1)
               XM2(2)=YTUN(N,M+1)
               XM2(3)=ZTUN(N,M+1)
               XM3(1)=XTUN(N+1,M+1)
               XM3(2)=YTUN(N+1,M+1)
               XM3(3)=ZTUN(N+1,M+1)
               XM4(1)=XTUN(N+1,M)
               XM4(2)=YTUN(N+1,M)
               XM4(3)=ZTUN(N+1,M)

               IMR0=0
               ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
     *              +ABS(XM2(3)-XM1(3))
               ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
     *              +ABS(XM3(3)-XM2(3))
               ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
     *              +ABS(XM4(3)-XM3(3))
               ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
     *              +ABS(XM1(3)-XM4(3))
               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)

               IF(ERRMIN.LE.1.0E-6) THEN
                  WRITE(*,2200) M,N
 2200             FORMAT(' ...... triangular panels (m,n) -- > ',
     *                 2I4,' (tunnel) ')
                  IMR0=1
               END IF

               DO KK=1,NBLADE
                  DO I=1,NPANEL
                     XLOC=0.
                     YLOC=0.
                     ZLOC=0.
                     DO K=1,3
                        XLOC=XLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,1,K)
                        YLOC=YLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,2,K)
                        ZLOC=ZLOC+(XCTP(I,K,KK)-XCT(J,K))*DIR(J,3,K)
                     ENDDO
                     IMR=IMR0

                     CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                    FDX,FDY,FDZ,0,IMR)

                     IF(IMR.EQ.2) THEN
                        DO IXYZ=1,3
                           XMC(IXYZ)=XCTP(I,IXYZ,KK)
                        ENDDO
                        CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                     END IF

                     IF(I.EQ.J .AND. KK.EQ.1) THEN
                        IF(IMR.EQ.1) THEN
                           FD=TWOPI
                        END IF
                        IF(FD.LT.0.) THEN
                           FD=4.0*PI+FD
                        END IF
                        IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.5.0) THEN
C                           WRITE(*,*) 'INFCOF-E',NTSTEP,KK,I,M,N,IMR,FD
                           FD=2.0*PI
                        END IF
                     ELSE
                        IF(ABS(FD).GT.6.28) THEN
C                           WRITE(*,*) 'INFCOF-F',NTSTEP,KK,I,M,N,IMR,FD
                           FD=0.0
                        END IF
                        FD = -FD
                     END IF

                     IF(IDUCT .NE. 0.AND. IDOPT .EQ. 1) THEN
                        IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                           FD = 0.0
                           FS = 0.0
                        ENDIF
                     ENDIF

                     A(I)=FD
                     TEMP1(I) = FS

                  ENDDO
                  CALL WRITE1(IUMA,A,NPANEL)
                  CALL WRITE1(IUMB,TEMP1,NPANEL)
               ENDDO

            ENDDO
         ENDDO
      END IF

C-----------------------------------------------------------------------
C     Compute influence coefficients due to the tip vortex cavity
C     Do not include effects of the tip vortex and tip bulb
C     corresponding to the other blades.
C     (Consider only KEY Blade's effect)
C-----------------------------------------------------------------------

      IF(IAN.EQ.2) THEN

         DO N=1,nthx
            DO M=1,mcvt

               J=INDEXT(N,M)

               IMR0=0

               DO K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
               ENDDO
               DO K=1,15
                  S(K)=SS(J,K)
               ENDDO

               XM1(1)=XCH(N+1,M+1)
               XM1(2)=YCH(N+1,M+1)
               XM1(3)=ZCH(N+1,M+1)
               XM2(1)=XCH(N+1,M)
               XM2(2)=YCH(N+1,M)
               XM2(3)=ZCH(N+1,M)
               XM3(1)=XCH(N,M)
               XM3(2)=YCH(N,M)
               XM3(3)=ZCH(N,M)
               XM4(1)=XCH(N,M+1)
               XM4(2)=YCH(N,M+1)
               XM4(3)=ZCH(N,M+1)

               CRLT = chrleps(j)

               ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
     %              +ABS(XM2(3)-XM1(3))
               ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
     %              +ABS(XM3(3)-XM2(3))
               ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
     %              +ABS(XM4(3)-XM3(3))
               ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
     %              +ABS(XM1(3)-XM4(3))
               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)
               IF(ERRMIN.LE.1.0E-6) THEN
                  IMR0=1
               END IF

               DO I=1,NPANEL
                  XLOC=0.
                  YLOC=0.
                  ZLOC=0.

                  DO K=1,3
                     XLOC=XLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,3,K)
                  ENDDO

                  IMR=IMR0

                  CALL RPAN(XLOC,YLOC,ZLOC,CRLT,FS,FD,FSX,FSY,
     *                 FDX,FDY,FDZ,0,IMR)

                  IF(IMR.EQ.2) THEN
                     DO IXYZ=1,3
                        XMC(IXYZ)=XCTP(I,IXYZ,1)
                     ENDDO
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF

                  IF(I.EQ.J) THEN
                     IF(IMR.EQ.1) THEN
                        FD=TWOPI
                     END IF
                     IF(FD.LT.0.) THEN
                        FD=4.0*PI+FD
                     END IF

                     IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN
                        FD=2.0*PI
                     END IF
                  ELSE
                     IF(ABS(FD).GT.6.28) THEN
                        FD=0.0
                     END IF
                  END IF

                  IF(IDUCT .NE. 0 .AND. IDOPT .EQ. 1) THEN
                     IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                        FD = 0.0
                        FS = 0.0
                     ENDIF
                  ENDIF

                  A(I)=FD
                  TEMP1(I)=FS

               ENDDO
               CALL WRITE1(IUMA,A,NPANEL)
               CALL WRITE1(IUMB,TEMP1,NPANEL)
            ENDDO
         ENDDO

         DO N=1,ncvx
            DO M=1,mcvt

               J=INDEXC(N,M)

               IMR0 = 0

               DO K=1,4
                  XV(K)=XVP(J,K)
                  YV(K)=YVP(J,K)
                  SIDE(K)=SID(J,K)
               ENDDO
               DO K=1,15
                  S(K)=SS(J,K)
               ENDDO

               XM1(1)=XVC(N+1,M+1)
               XM1(2)=YVC(N+1,M+1)
               XM1(3)=ZVC(N+1,M+1)
               XM2(1)=XVC(N+1,M)
               XM2(2)=YVC(N+1,M)
               XM2(3)=ZVC(N+1,M)
               XM3(1)=XVC(N,M)
               XM3(2)=YVC(N,M)
               XM3(3)=ZVC(N,M)
               XM4(1)=XVC(N,M+1)
               XM4(2)=YVC(N,M+1)
               XM4(3)=ZVC(N,M+1)

               ER1=ABS(XM2(1)-XM1(1))+ABS(XM2(2)-XM1(2))
     %              +ABS(XM2(3)-XM1(3))
               ER2=ABS(XM3(1)-XM2(1))+ABS(XM3(2)-XM2(2))
     %              +ABS(XM3(3)-XM2(3))
               ER3=ABS(XM4(1)-XM3(1))+ABS(XM4(2)-XM3(2))
     %              +ABS(XM4(3)-XM3(3))
               ER4=ABS(XM1(1)-XM4(1))+ABS(XM1(2)-XM4(2))
     %              +ABS(XM1(3)-XM4(3))
               ERRMIN=AMIN1(ER1,ER2,ER3,ER4)
               IF(ERRMIN.LE.1.0E-6) THEN
                  IMR0=1
               END IF

               CRLT = chrleps(j)

               DO I=1,NPANEL
                  XLOC=0.
                  YLOC=0.
                  ZLOC=0.
                  DO K=1,3
                     XLOC=XLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,1)-XCT(J,K))*DIR(J,3,K)
                  ENDDO

                  IMR=IMR0

                  CALL RPAN(XLOC,YLOC,ZLOC,CRLT,FS,FD,FSX,FSY,
     *                 FDX,FDY,FDZ,0,IMR)

                  IF(IMR.EQ.2) THEN
                     DO IXYZ=1,3
                        XMC(IXYZ)=XCTP(I,IXYZ,1)
                     ENDDO
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF
                  IF(I.EQ.J) THEN
                     IF(IMR.EQ.1) THEN
                        FD=TWOPI
                     END IF
                     IF(FD.LT.0.) THEN
                        FD=4.0*PI+FD
                     END IF

                     IF(ABS(FD).GT.7.0.OR.ABS(FD).LT.6.0) THEN
                        FD=2.0*PI
                     END IF
                  ELSE
                     IF(ABS(FD).GT.6.28) THEN
                        FD=0.0
                     END IF
                  END IF
                  A(I)=FD
                  TEMP1(I)=FS
               ENDDO

               IF(IDUCT .NE. 0 .AND. IDOPT .EQ. 1) THEN
                  IF(I .GT. IDDUCT2 .AND. I .LE. IDDUCT3) THEN
                     FD = 0.0
                     FS = 0.0
                  ENDIF
               ENDIF

               CALL WRITE1(IUMA,A,NPANEL)
               CALL WRITE1(IUMB,TEMP1,NPANEL)

            ENDDO
         ENDDO
      ENDIF

CSH--REPLACE, NBLADE=1-------------------------
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------

      RETURN
      END SUBROUTINE
