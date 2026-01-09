
C ===============================================================
        SUBROUTINE INTERBLADE( XXX,ZZZ,NBLADE,NH,MHBT,PITL,PITT)
C ===============================================================
C
C       Transverse grid line and circumferential grid points are
C       generated : Hub & DUCT
C
C       Programmed by H.S. Lee
C-------------------------------------------------------------
        use m_NACAIO
        INCLUDE 'PARAM.INC'
!       COMMON /NACAIO/ XKEY(NBPZ),YKEY(NBPZ),ZKEY(NBPZ)
!    #                  ,XNEAR(NBPZ),ZNEAR(NBPZ),CROSS(3,NBPZ)
        DIMENSION COFT(4,2),ABUG(4,4),CBUG(4)
        DIMENSION COFM(4,MHPZ)
        DIMENSION XXX(NBPZ,MHPZ+1),ZZZ(NBPZ,MHPZ+1)
        dimension sdum(301), xxdum(301), zzdum(301),
     %            xxdumco(1200), zzdumco(1200)
C
C       Sectons will be virtually corrected as following scheme.
C
C          1. Its shape is originally loacted on the conical cylinder.
C             It is projected on a circular cylinder for convenience's sake
C             During the step, axial coordinates (X) is not changed but
C             vertical (Y) and horizontal (Z) ones is corrected.
C          2. The surface on the circular cylinder is unrolled to planar
C             plane so that the 3 dimensional coordinate is translated
C             to 2 dimensinal one.
C          3. The coordinates in global coordinate system is transformed
C             to local one whose origin is on the mid-chord position of
C             key section and which is rotated by tip pitch angle
C          4. Transverse grid lines connecting nodes on back and face of
C             key and near section are generated using cubic function.
C          5. During all the steps, the checks are performed transverse grid
C             line crosses section and near grid line or not. If crossing
C             occurs, slope of the cubic function is adjusted and re-build
C             the function
C          6. After that, circumferential grid points are calculated and
C             reverse transform and projection is fulfilled.
C
C       OUTPUT : XXX(N,M) , ZZZ(N,M)  , N=1,NN , M=1...MHBT1
C
C                XXX : Axial Position
C                ZZZ : Angle
C
C
C       NOTE : when pitch angle is greater than zero , i.e.,
C              right handed vane case, transverse grid line is connected
C              from back of key section to face of near section.
C
C                                              L.E.
C                               /  F        B   /
C                             /<-A--------A-- /
C                           /  C        C   /
C                         /  E        K   /
C                                      T.E.
C
C              when pitch angle is less than zero , i.e.,
C              left handed vane case, transverse grid line is connected
C              from face of key section to back of near section.
C
C                                     L.E.
C                       \  B         F  \
C                         \<-A-------- A  \
C                           \  C         C  \
C                             \  K         E  \
C                                            T.E.
C
C------ Data preparations for Key and near section
C

        ZERO = 0.0
        HALF = 0.5
        ONE = 1.0
        TWO = 2.0
        PI = ACOS(-1.)
        TWOPI = TWO*PI

        NN = NH

        write(*,*) ' ... Inter-Blade Shaping '

        DTHET2 = TWOPI / FLOAT(NBLADE)

        DO I=1,2*NN+1
           XKEY(I) = CROSS(1,I)
           YKEY(I) = CROSS(2,I)
           ZKEY(I) = CROSS(3,I)
           XNEAR(I) = CROSS(1,I)
        ENDDO

C
C------ Find Mean Radius : Mean of Radii of L.E. and T.E.
C -- For Constant radius cyinder : Rmean = RHUB
C
        RMEAN = HALF *
     #        ( SQRT( YKEY(NN+1)**2 + ZKEY(NN+1)**2 )
     #        + SQRT( YKEY(1)**2    + ZKEY(1)**2    ) )

C
C------ Project sections on the circular cylinder
C------ Expand the circular cylinder into 2-D Plane: 3D (X,Y,Z) ---> 2D (X,Z)
C
C  X -> X, Y & Z -> R * theta
C
        DO I=1,2*NN+1
           TBUG = ATAN2( ZKEY(I) , YKEY(I) )
           ZKEY(I) = RMEAN * TBUG
           TBUG = TBUG + DTHET2
           ZNEAR(I) = RMEAN * TBUG
        ENDDO
C
C------ Find Tip Pitch Angle
C
        IF( ZKEY(NN+1).EQ.ZKEY(1) ) THEN
           TPHI = HALF * PI
        ELSE
           TPHI = ATAN( ( XKEY(NN+1) - XKEY(1) )
     #                 / ( ZKEY(NN+1) - ZKEY(1) ) )
        ENDIF

C
C------ Coordinate Transform : Global ---> Local
C
        XORG = XKEY(1)
        ZORG = ZKEY(1)

        DO I=1,2*NN+1
           XTEMP = XKEY(I)
           ZTEMP = ZKEY(I)
           IF(TPHI.GT.ZERO) THEN
              XKEY(I) =  ( XTEMP-XORG ) * SIN( TPHI )
     #                  +( ZTEMP-ZORG ) * COS( TPHI )
              ZKEY(I) = -( XTEMP-XORG ) * COS( TPHI )
     #                  +( ZTEMP-ZORG ) * SIN( TPHI )
           ELSE
              XKEY(I) =  ( XTEMP-XORG ) * SIN( TPHI )
     #                  -( ZTEMP-ZORG ) * COS( TPHI )
              ZKEY(I) =  ( XTEMP-XORG ) * COS( TPHI )
     #                  +( ZTEMP-ZORG ) * SIN( TPHI )
           ENDIF
           XTEMP = XNEAR(I)
           ZTEMP = ZNEAR(I)
           IF(TPHI.GT.ZERO) THEN
              XNEAR(I) =  ( XTEMP-XORG ) * SIN( TPHI )
     #                  +( ZTEMP-ZORG ) * COS( TPHI )
              ZNEAR(I) = -( XTEMP-XORG ) * COS( TPHI )
     #                  +( ZTEMP-ZORG ) * SIN( TPHI )
           ELSE
              XNEAR(I) =  ( XTEMP-XORG ) * SIN( TPHI )
     #                  -( ZTEMP-ZORG ) * COS( TPHI )
              ZNEAR(I) =  ( XTEMP-XORG ) * COS( TPHI )
     #                  +( ZTEMP-ZORG ) * SIN( TPHI )
           ENDIF
        ENDDO

C
C------ Initialize 4 X 4 matrix for cubic function
C
C  F =  A*x**3 + B*x**2 + C*X + D
C
C    ==> [ABUG][COEFT] = [CBUG]
C
C
        DO I=1,4
           DO J=1,4
              ABUG(I,J) = ZERO
           ENDDO
           CBUG(I) = ZERO
           DO J=1,2
              COFT(I,J) = ZERO
           ENDDO
        ENDDO

C
C       When TPHI > 0
C               N       : Back of key section
C               L       : Face of near section
        NLE = NN + 1
        NTE = 2 * NN + 1
        ID = 1
C
C       When TPHI < 0
C               N       : Face of key section
C               L       : Back of near section
C
        IF( TPHI.LT.ZERO ) NTE = 1
        IF( TPHI.LT.ZERO ) ID = -1
        LLE = 2 * NN + 2 - NLE
        LTE = 2 * NN + 2 - NTE
C
C------ Find Index at which ZKEY is max. and ZNEAR is min.
C
        BIG = ZKEY(1)
        IKEY = 1
        DO I=2,2*NN+1
           IF( ZKEY(I).GT.BIG ) THEN
              BIG = ZKEY(I)
              IKEY = I
           ENDIF
        ENDDO

        SMALL = ZNEAR(1)
        INEAR = 1
        DO I=2,2*NN+1
           IF( ZNEAR(I).LT.SMALL ) THEN
              SMALL = ZNEAR(I)
              INEAR = I
           ENDIF
        ENDDO

C
C------ Transverse Grid Line and Grid Points at Edges
C
        write(*,*) ' ... Transverse Grid at Edges'

        DO 1700 K=1,2
C
C--------- K=1 : Leading Edge , K=2 : Trailing Edge
C
           I = NLE
           I2 = LLE
           IF( K.EQ.2 ) THEN
              I = NTE
              I2 = LTE
           ENDIF
           ABUG(1,1) = ZKEY(I)**3
           ABUG(1,2) = ZKEY(I)**2
           ABUG(1,3) = ZKEY(I)
           ABUG(1,4) = ONE
           ABUG(2,1) = ZNEAR(I2)**3
           ABUG(2,2) = ZNEAR(I2)**2
           ABUG(2,3) = ZNEAR(I2)
           ABUG(2,4) = ONE
           ABUG(3,1) = 3.0 * ZKEY(I)**2
           ABUG(3,2) = TWO * ZKEY(I)
           ABUG(3,3) = ONE
           ABUG(3,4) = ZERO
           ABUG(4,1) = 3.0 * ZNEAR(I2)**2
           ABUG(4,2) = TWO * ZNEAR(I2)
           ABUG(4,3) = ONE
           ABUG(4,4) = ZERO
           CBUG(1) = XKEY(I)
           CBUG(2) = XNEAR(I2)
C
C--------- Solve matrix to get coefficents of cubic function
C
           SLP1 = ZERO
           SLP2 = ZERO

           IF( K.EQ.1 ) THEN                   ! L.E.

              I5 = I + 1
              I6 = I2 - 1
              IF( TPHI.LT.ZERO ) THEN
                 I5 = I - 1
                 I6 = I2 + 1
              ENDIF

              SLP1 = ( XKEY(I5) - XKEY(I) )
     #             / ( ZKEY(I5) - ZKEY(I) )

              SLP2 = ( XKEY(I) - XKEY(I6) )
     #             / ( ZKEY(I) - ZKEY(I6) )

              TBUG = - PITL + ATAN( SLP1 )
              CBUG(3) = TAN( TBUG )

              TBUG =   PITL + ATAN( SLP2 )
              CBUG(4) = TAN( TBUG )

            ELSE

              I5 = I - 1
              I6 = I2 + 1
              IF( TPHI.LT.ZERO ) THEN
                 I5 = I + 1
                 I6 = I2 - 1
              ENDIF

              SLP1 = ( XKEY(I5) - XKEY(I) )
     #             / ( ZKEY(I5) - ZKEY(I) )

              SLP2 = ( XKEY(I) - XKEY(I6) )
     #             / ( ZKEY(I) - ZKEY(I6) )

              TBUG =   PITT + ATAN( SLP1 )
              CBUG(3) = TAN( TBUG )

              TBUG = - PITT + ATAN( SLP2 )
              CBUG(4) = TAN( TBUG )
           ENDIF
C
C--------- Solve matrix to get coefficents of cubic function
C
C          When number of iteration for modifying cubic functions
C          is greater than 30 , it means that thickness is too thick
C          or pitch is too small. So internal modifications of them
C          is performed to fit them.
C
           LCOUNT = 0
 1800      IF(LCOUNT.GE.10000) CALL ERRMSG( 1 )
           CALL JORDAN( ABUG , CBUG , IERR , COFT(1,K) )

C
C--------- Check whether transverse grid line crosses key section or not
C          and adjust slope at key section
C
           IF( K.EQ.1 ) THEN
              IDD = ID
           ELSE
              IDD = -ID
           ENDIF
           IERR = 0
           DO 1900 J=I,IKEY,IDD
              IF( J.EQ.I ) GOTO 1900
              ZBUG1 = COFT(1,K) * ZKEY(J)**3
     #              + COFT(2,K) * ZKEY(J)**2
     #              + COFT(3,K) * ZKEY(J)
     #              + COFT(4,K)

              IF( K.EQ.1 ) THEN

                 IF( ZBUG1.GE.XKEY(J) ) THEN            ! Crossing Occurs
                     SBUG = ( SLP1 - CBUG(3) ) / 100.0
                     CBUG(3) = CBUG(3) - SBUG
                     IERR = 1
                 ENDIF

              ELSE

                 IF( ZBUG1.LE.XKEY(J) ) THEN            ! Crossing Occurs
                     SBUG = ( SLP1 - CBUG(3) ) / 100.0
                     CBUG(3) = CBUG(3) + SBUG
                     IERR = 1
                 ENDIF

              ENDIF

              IF(IERR.EQ.1) THEN
                 IF(LCOUNT.EQ.49)
     #            write(*,*) ' ........ Key section crossing'
                 LCOUNT = LCOUNT + 1
                 GOTO 1800
              ENDIF

 1900      CONTINUE
C
C--------- Check whether transverse grid line crosses near section or not
C          and adjust slope at near section
C
           IERR = 0
           DO 1950 J=I2,INEAR,-IDD
              IF( J.EQ.I2 ) GOTO 1950
              ZBUG1 = COFT(1,K) * ZNEAR(J)**3
     #              + COFT(2,K) * ZNEAR(J)**2
     #              + COFT(3,K) * ZNEAR(J)
     #              + COFT(4,K)

              IF( K.EQ.1 ) THEN

                 IF( ZBUG1.GE.XNEAR(J) ) THEN           ! Crossing Occurs
                     SBUG = ( SLP2 - CBUG(4) ) / 100.0
                     CBUG(4) = CBUG(4) + SBUG
                     IERR = 1
                 ENDIF

              ELSE

                 IF( ZBUG1.LE.XNEAR(J) ) THEN           ! Crossing Occurs
                     SBUG = ( SLP2 - CBUG(4) ) / 100.0
                     CBUG(4) = CBUG(4) - SBUG
                    IERR = 1
                 ENDIF

              ENDIF
              IF(IERR.EQ.1) THEN
                 IF(LCOUNT.EQ.49)
     #            write(*,*) ' ........ Near section crossing'
                 LCOUNT = LCOUNT + 1
                 GOTO 1800
              ENDIF
 1950      CONTINUE

C
C--------- Calculates Transverse Grid Points
C
 1910      write(*,*) ' ... Spacing Transeverse Grid at Edges'

           IF( K.EQ.1 ) N = 1                   ! L.E.
           IF( K.EQ.2 ) N = NN + 1              ! T.E
           ZDEL = ( ZNEAR(I2) - ZKEY(I) ) / FLOAT( MHBT )

           DO 2200 M=1,MHBT+1
              IF(M.EQ.1) THEN
                 ZZZ(N,M) = ZKEY(I)
                 XXX(N,M) = XKEY(I)
              ELSE
                IF(M.EQ.MHBT+1) THEN
                 ZZZ(N,M) = ZNEAR(I2)
                 XXX(N,M) = XNEAR(I2)
                ELSE
                 ZZZ(N,M) = ZZZ(N,M-1) + ZDEL
                 XXX(N,M) = COFT(1,K) * ZZZ(N,M)**3
     #                 + COFT(2,K) * ZZZ(N,M)**2
     #                 + COFT(3,K) * ZZZ(N,M)
     #                 + COFT(4,K)
                ENDIF
             ENDIF
 2200      CONTINUE

 1700   CONTINUE

C
C------ Build Grid points
C
 2290   N = 1
        DO 2300 I=NLE,NTE,ID
           IF( I.EQ.NLE .OR. I.EQ.NTE ) GOTO 2300
           N = N + 1
           I2 = 2 * NN + 2 - I
           DO 2400 M=1,MHBT+1,MHBT

              IF(M.EQ.1) THEN
                 ZZZ(N,M) = ZKEY(I)
                 XXX(N,M) = XKEY(I)
              ELSE
                 ZZZ(N,M) = ZNEAR(I2)
                 XXX(N,M) = XNEAR(I2)
              ENDIF

 2400      CONTINUE

 2300   CONTINUE
C
C------ Find Index at Mid-chord
C
        IKEY = ( NN + 2 ) / 2
        INEAR = IKEY
C
C------ Meridional Grid Line
C
        write(*,*) ' ... Spacing Meridional Grid'

        ZDEL = ( ZZZ(INEAR,MHBT+1) - ZZZ(IKEY,1) )
     #       / FLOAT(MHBT )

        DO 4100 M=2,MHBT
           ABUG(1,1) = XXX(NN+1,M)**3
           ABUG(1,2) = XXX(NN+1,M)**2
           ABUG(1,3) = XXX(NN+1,M)
           ABUG(1,4) = ONE
           ABUG(2,1) = XXX(1,M)**3
           ABUG(2,2) = XXX(1,M)**2
           ABUG(2,3) = XXX(1,M)
           ABUG(2,4) = ONE
           ZBUG = ZZZ(IKEY,1) + ZDEL * ( M - 1 )

           XBUG = HALF * ( XXX(1,M) + XXX(NN+1,M) )
           ABUG(3,1) = XBUG**3
           ABUG(3,2) = XBUG**2
           ABUG(3,3) = XBUG
           ABUG(3,4) = ONE
           ABUG(4,1) = 3.0 * XBUG**2
           ABUG(4,2) = TWO * XBUG
           ABUG(4,3) = ONE
           ABUG(4,4) = ZERO
           CBUG(1) = ZZZ(NN+1,M)
           CBUG(2) = ZZZ(1,M)
           CBUG(3) = ZBUG
           CBUG(4) = ZERO

           LCOUNT = 0
 4150      IF(LCOUNT.GE.1000) CALL ERRMSG( 1 )
           CALL JORDAN( ABUG , CBUG , IERR , COFM(1,M) )
C
C--------- Check whether Meridional Grid Line crosses section or not
C
           M1 = M - 1

 4170      DO 4200 N=1,NN+1
              ZBUG1 = COFM(1,M) * XXX(N,M1)**3
     #              + COFM(2,M) * XXX(N,M1)**2
     #              + COFM(3,M) * XXX(N,M1)
     #              + COFM(4,M)
              IERR = 0
              IF( M1.EQ.M-1 ) THEN
                 IF( ZBUG1.LE.ZZZ(N,M1) ) THEN
                     IF(N.GE.IKEY) THEN
                         CBUG(4) = CBUG(4) + 0.02
                     ELSE
                         CBUG(4) = CBUG(4) - 0.02
                     ENDIF
                     IERR = 1
                 ENDIF
              ELSE
                 IF( ZBUG1.GE.ZZZ(N,M1) ) THEN
                     IF(N.GE.INEAR) THEN
                         CBUG(4) = CBUG(4) - 0.02
                     ELSE
                         CBUG(4) = CBUG(4) + 0.02
                     ENDIF
                     IERR = 1
                 ENDIF
              ENDIF
              IF(IERR.EQ.1) THEN
                 IF(LCOUNT.EQ.49) write(*,'(a,i2,a,a,i2,a)')
     #                       ' ........ ',M,'-th Meridonal'
     #                      ,' Grid Crosses',M1,'-th One'
                 LCOUNT = LCOUNT + 1
                 GOTO 4150
              ENDIF

 4200      CONTINUE

           IF(M.EQ.MHBT.AND.M1.EQ.MHBT-1) THEN
              M1 = MHBT+1
              GOTO 4170
           ENDIF

C -----------------------------------------------------
 4400      xdel = (xxx(nn+1,m) - xxx(1,m) ) / real(200)
           sdum(1) = 0.0
           xxdum(1) = xxx(1,m)
           zzdum(1) = zzz(1,m)
           do m1 = 1 , 200
              xxdum(m1+1) = xxdum(m1) + xdel
              zzdum(m1+1) = cofm(1,m) * xxdum(m1+1)**3
     %                    + cofm(2,m) * xxdum(m1+1)**2
     %                    + cofm(3,m) * xxdum(m1+1) + cofm(4,m)
           enddo

           do m1 =1 , 200
              xlen = sqrt ( (xxdum(m1+1) - xxdum(m1))**2 +
     %                      (zzdum(m1+1) - zzdum(m1))**2 )
              sdum(m1+1) = sdum(m1) + xlen
           enddo

C ---- YiranSu ---- improve the inter-blade interpolation scheme to
C                   avoid extream panels
C           call spacenew(1, nn+1, xkey, zero, sdum(201))
      NHP=NH+1
      XKEY(1)=0.0E0
      DO II=2,NHP
      XKEY(II)=XKEY(II-1)+SQRT((CROSS(1,NHP+1-II)-CROSS(1,NHP-II+2))**2
     *                        +(CROSS(2,NHP+1-II)-CROSS(2,NHP-II+2))**2
     *                        +(CROSS(3,NHP+1-II)-CROSS(3,NHP-II+2))**2)
      END DO
      TTLEN=XKEY(NHP)
      DO II=2,NHP
         XKEY(II)=XKEY(II)/TTLEN*SDUM(201)
      END DO
C Done!


           CALL UGLYDK(201,1,1,sdum,xxdum,0.0,0.0,xxdumco)
           CALL UGLYDK(201,1,1,sdum,zzdum,0.0,0.0,zzdumco)

           CALL EVALDK(201,nn+1,sdum,xkey,xxdum,xxdumco)
           CALL EVALDK(201,nn+1,sdum,xkey,zzdum,zzdumco)

           do n = 1 , nn+1
              xxx(n,m) = xxdum(n)
              zzz(n,m) = zzdum(n)
           enddo

 4100   CONTINUE

C
C------ Coordinate Re-Transform : Local ---> Global
C
        DO 2500 N=1,NN+1
           DO 2600 M=1,MHBT+1
              XTEMP = XXX(N,M)
              ZTEMP = ZZZ(N,M)
              IF(TPHI.GT.ZERO) THEN
                 XXX(N,M) = XTEMP * SIN( TPHI )
     #                    - ZTEMP * COS( TPHI ) + XORG
                 ZZZ(N,M) = XTEMP * COS( TPHI )
     #                    + ZTEMP * SIN( TPHI ) + ZORG
              ELSE
                 XXX(N,M) = XTEMP * SIN( TPHI )
     #                    + ZTEMP * COS( TPHI ) + XORG
                 ZZZ(N,M) =-XTEMP * COS( TPHI )
     #                    + ZTEMP * SIN( TPHI ) + ZORG
              ENDIF

 2600      CONTINUE
           DO 2700 M=1,MHBT+1
              ZZZ(N,M) = ZZZ(N,M) / RMEAN
 2700      CONTINUE
 2500   CONTINUE

        RETURN
        END


C ===================================
        SUBROUTINE ERRMSG( I )
C ===================================
        IF( I.EQ.1 ) THEN
           WRITE(*,*) ' .... Error ... Too Thick Section '
     #               ,' or Too Small Pitch '
        ENDIF
        STOP
        END


C====================================================================
        SUBROUTINE JORDAN(AA,BB,IERR,X)
C====================================================================
C       IERR =   0 : NO ERROR
C       IERR = 999 : A MATRIX IS SINGULAR
C====================================================================
C
        DIMENSION A(4,4) , B(4) , X(4)
        DIMENSION AA(4,4) , BB(4)
        DIMENSION IPIV(4),INDXR(4),INDXC(4)

        N = 4
        IERR = 0
        DO 11 J=1,N
           IPIV(J) = 0
 11     CONTINUE
        DO 177 I=1,N
           DO 188 J=1,N
              A(I,J) = AA(I,J)
 188       CONTINUE
           B(I) = BB(I)
 177    CONTINUE
        DO 22 I=1,N
           BIG = 0.0
           DO 13 J=1,N
              IF(IPIV(J).NE.1) THEN
                 DO 12 K=1,N
                    IF(IPIV(K).EQ.0) THEN
                       IF(ABS(A(J,K)).GE.BIG) THEN
                          BIG = ABS(A(J,K))
                          IROW = J
                          ICOL = K
                       ENDIF
                    ELSEIF(IPIV(K).GT.1) THEN
                       IERR = 999
                       RETURN
                    ENDIF
 12              CONTINUE
              ENDIF
 13        CONTINUE
           IPIV(ICOL) = IPIV(ICOL) + 1
           IF(IROW.NE.ICOL) THEN
              DO 14 M=1,N
                 DUM = A(IROW,M)
                 A(IROW,M) = A(ICOL,M)
                 A(ICOL,M) = DUM
 14           CONTINUE
              DUM = B(IROW)
              B(IROW) = B(ICOL)
              B(ICOL) = DUM
           ENDIF
           INDXR(I) = IROW
           INDXC(I) = ICOL
           IF(A(ICOL,ICOL).EQ.0.0) THEN
              IERR = 999
              RETURN
           ENDIF

           PIVINV = 1.0 / A(ICOL,ICOL)
           A(ICOL,ICOL) = 1.0

           DO 16 M=1,N
              A(ICOL,M) = A(ICOL,M) * PIVINV
 16        CONTINUE
           B(ICOL) = B(ICOL) * PIVINV

           DO 21 MM=1,N
              IF(MM.NE.ICOL) THEN
                 DUM = A(MM,ICOL)
                 A(MM,ICOL) = 0.0
                 DO 18 M=1,N
                    A(MM,M) = A(MM,M) - A(ICOL,M) * DUM
 18              CONTINUE
                 B(MM) = B(MM) - B(ICOL) * DUM
              ENDIF
 21        CONTINUE
 22     CONTINUE

        DO 24 M=N,1,-1
           IF(INDXR(M).NE.INDXC(M)) THEN
              DO 23 K=1,N
                 DUM = A(K,INDXR(M))
                 A(K,INDXR(M)) = A(K,INDXC(M))
                 A(K,INDXC(M)) = DUM
 23           CONTINUE
           ENDIF
 24     CONTINUE

        DO 91 I=1,N
           X(I) = B(I)
 91     CONTINUE

        RETURN
        END


C ==========================================================
      SUBROUTINE SPACENEW( IK , N , X , XI , XE )
C ==========================================================
C
C      SPACING
C
C      IK = 1       : COSINE SPACING
C      IK = 2       : HALF-COSINE SPACING FROM XI TO XE
C      IK = -2 : HALF-COSINE SPACING FROM XE TO XI
C      N       : NUMBER OF POINTS TO BE GOTTEN ( N-1 PANELS )
C      X(N)       : SPACED VALUE
C      XI         : INITIAL VALUE
C      XE         : FINAL VALUE    ( XI < XE )
C
C-------------------------------------------------------------
      DIMENSION X(N)
      ONE = 1.0
      HALF = 0.5
      PI = ACOS( -ONE )

      IF( XI.GE.XE ) THEN
         WRITE(*,*) ' ... ERROR ... INVALID XI,XE IN SPACE '
         WRITE(*,*) '               XI = ',XI
         WRITE(*,*) '               XE = ',XE
         STOP
      ENDIF
      IF( N.EQ.1 ) THEN
         WRITE(*,*) ' ... ERROR ... INVALID N IN SPACE '
         WRITE(*,*) '               N  = ',N
         STOP
      ENDIF
      BUG = XE - XI
      IF( IK.EQ.1 ) THEN
         DO 1100 I=1,N-1
              FRACT = ( I-1 ) / FLOAT( N-1 )
            X(I) = XI + ( HALF - HALF * COS( PI*FRACT ) ) * BUG
 1100         CONTINUE
         X(N) = XE
      ENDIF
      IF( IK.EQ.2 ) THEN
         DO 1200 I=1,N-1
              FRACT = ( I-1 ) / FLOAT( N-1 )
            X(I) = XI + ( ONE - COS( HALF*PI*FRACT ) ) * BUG
 1200         CONTINUE
         X(N) = XE
      ENDIF
      IF( IK.EQ.-2 ) THEN
         DO 1300 I=1,N-1
              FRACT = ONE + ( I-1 ) / FLOAT( N-1 )
            X(I) = XI - COS( HALF*PI*FRACT ) * BUG
 1300         CONTINUE
         X(N) = XE
      ENDIF

      RETURN
      END

