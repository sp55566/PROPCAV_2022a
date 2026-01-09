C ===================================================
      SUBROUTINE DUCTINPUT
C ===================================================
      INCLUDE 'PUFCAV.INC'
      COMMON /DIMAGE/ DGAPSIZE,ENDFOIL
      DIMENSION SSD(101),DYTC(101),DTCC(101)
      DIMENSION SSDC(201),DYCCC(201),DCUBIC(1000),D2YDX(1001),DYCC(1001)
      DIMENSION X1(101),Y1(101)
      DIMENSION xbb1(100),ybb1(100),cubb(400)

C-- Read duct data file
      
      OPEN(9,FILE='duct.dat',STATUS='OLD')

C =====================================================================================================
C     XLED, YLED : Coordinates of duct Leading and Trailing edges in propeller coordinates (x/R, r/R)
C     DCHORD : Duct Chord Length (normalized by Propeller Radius)
C     DUCTANG : Duct Angle of Attack (Degrees)
C     IDTHK, IDCAM : Duct Thick and Camber distribution
C        if IDTHK and IDCAM = -1, read section data
C     DTHICK, DCAMBER : tmax/C, fmax/C
C     NDFWD, NDAFT : Number of panels on duct before and after propeller blade.
C     NDBACK : No. of panels on duct backside
C     IDGEO  =  0 Generate duct panel using propeller Pitch
C               1 Read Pitch angle for duct panel
C     IDOPT = 0 Solve duct+Propeller
C           = 1 Solver Duct and propeller seperately
C     DUCTPT : Duct pitch angle
C     IDALIGN : 0 : Blade wake not align
C               1 : Blade wake align
C =====================================================================================================
      
      READ(9,*) XLED, YLED
      READ(9,*) DCHORD
      READ(9,*) DUCTANG, DWAKEL
      READ(9,*) IDTHK, IDCAM, IDGEO, IDOPT
      READ(9,*) DTHICK, DCAMBER, DUCTGAP
      READ(9,*) NDUCT, MDUCT 
      READ(9,*) DUCTPT, IDALIGN, IREPANEL, NWAKEP

C/s S.N.KIM | Repaneling Option
      IF(IREPANEL.EQ.1.AND.IDUCT.EQ.0) THEN
        WRITE(*,*) 'Repaneling opion is valid only for ducted propeller
     & cases (IDUCT=1).'
        STOP
      ENDIF
      IF(IREPANEL.EQ.1.AND.IDGEO.EQ.1) THEN
        WRITE(*,*) 'IDGEO=1, therefore repaneling option is ignored.'
        IREPANEL = 0
      ENDIF
C/e S.N.KIM | Repaneling Option

      DGAPSIZE = DUCTGAP

      IF(IDGEO .EQ. 1) THEN
         DUCTPT = DUCTPT * PI / 180.
      ENDIF

C ==================================================================
C -- IF IDTHK and IDCAM .eq. -1 : Read duct geometry from data (x,r)
C    Need (x,r) coordinates from trailing edge of face side 
C    to trailing edge of back side. 
C    The values of (x and r) should be normalized by propeller radius.
C ----------------------------------------------------------

      IF(IDTHK .EQ. -1 .AND. IDCAM .EQ. -1) THEN
         READ(9,*) NDDAT
         DO N = 1 , NDDAT
            READ(9,*) XID(N), ETAD(N)
         ENDDO
         CLOSE(9)

C -- Correct geometry wrt gap size

         NN1 = NDDAT/2+1

         DO N = 1, NN1
            N1 = NN1 - N + 1
            XBB1(N) = XID(N1)
            YBB1(N) = ETAD(N1)
         ENDDO

         CALL UGLYDK(NN1,1,1,XBB1,YBB1,ZERO,ZERO,CUBB)
         CALL EVALDKs(NN1,1,XBB1,0.0,RRR,CUBB)
        

         RRGAP = RRR - 1.0

         DO N = 1, NDDAT
            ETAD(N) = ETAD(N) - RRGAP + DUCTGAP
         ENDDO
        
         GO TO 1000
      ENDIF
C =================================================================

C -- Set unit duct geometry using cosine spacing

      PI = ACOS(-1.0)
      NNDD = 50
c      NNDD = 35
      NNDDP = NNDD + 1
      DELPP = PI / REAL(NNDD)

      COST = COS(DUCTANG * PI / 180.)
      SINT = SIN(DUCTANG * PI / 180.)

      DO N = 1 , NNDDP
         SSD(N) = HALF*(1. - COS((N-1)*DELPP))
      ENDDO

      IF(DTHICK .EQ. 0.0) THEN
         DO N = 1 , NNDDP
            DYTC(N) = 0.0
         ENDDO
      ELSE
         IF(IDTHK .EQ. 1) THEN
            CALL NACA66(NNDDP,DTHICK,DRLE,SSD,DYTC)
         ELSEIF(IDTHK .EQ. 2) THEN
            CALL RAE(NNDDP,DTHICK,SSD,DYTC)
         ELSE IF(IDTHK.EQ.3) THEN
            CALL NACA65A(NNDDP,DTHICK,DRLE,SSD,DYTC)
         ELSE IF(IDTHK.EQ.4) THEN 
            CALL NACA64A(NNDDP,DTHICK,DRLE,SSD,DYTC)
         ELSE IF(IDTHK.EQ.5) THEN
            CALL NACA00(NNDDP,DTHICK,DRLE,SSD,DYTC)
         ELSE IF(IDTHK.EQ.6) THEN
            CALL ELIPSE(NNDDP,DTHICK,DRLE,SSD,DYTC)
         ELSE IF(IDTHK.EQ.9) THEN
            CALL NACA16(NNDDP,DTHICK,DRLE,SSD,DYTC)
         ELSE IF(IDTHK.EQ.10) THEN
            CALL TKPB(NNDDP,DTHICK,DRLE,SSD,DYTC)
         ELSE IF(IDTHK .EQ. 99) THEN
            READ(9,*) NNDDP
            NNDD = NNDDP - 1
            DO N = 1 , NNDDP
               READ(9,*) SSD(N) , DYTC(N)
            ENDDO
         END IF            
      ENDIF
       
      IF(IDCAM .EQ. 0) THEN
         CALL A8ML(NNDDP,DCAMBER,SSD,DYCC,DTCC)
      ELSE IF(IDCAM.EQ.1) THEN
         CALL CBPB(NNDDP,DCAMBER,SSD,DYCC,DTCC)
      ELSE IF(IDCAM.EQ.99) THEN
         READ(9,*) NDC
         DO N = 1 , NDC
            READ(9,*) SSDC(N), DYCCC(N)
         ENDDO
C -- Interpolate to the points SSD(N)

         CALL UGLYDK(NDC,1,1,SSDC,DYCCC,ZERO,ZERO,DCUBIC)
         CALL EVALDK(NDC,NNDDP,SSDC,SSD,DYCC,DCUBIC)
         CALL DRIVDK(NDC,NNDDP,SSDC,SSD,DTCC,D2YDX,DCUBIC)
         DO N=1,NNDDP
            DTCC(N)=ATAN(DTCC(N))
         ENDDO
      END IF         

      CLOSE(9)

      DO N = 1 , NNDDP
      
         NBOT = NNDDP - N + 1
         NTOP = NNDDP + N - 1
            WRITE(*,*)N,NBOT,NTOP
         XXX = SSD(N) * DCHORD
         YYY = DYCC(N) * DCHORD 
      
         DXT = DYTC(N) * DCHORD * SIN(DTCC(N)) 
         DYT = DYTC(N) * DCHORD * COS(DTCC(N))

         XU = XXX - DXT
         XL = XXX + DXT

         YU = YYY + DYT
         YL = YYY - DYT
      
         XID(NTOP) = XU * COST + YU * SINT + XLED
         ETAD(NTOP) = -XU * SINT + YU * COST + YLED

         XID(NBOT) = XL * COST + YL * SINT + XLED
         ETAD(NBOT) = -XL * SINT + YL * COST + YLED

         XXX1 = XXX
         YYY1 = DYCC(N) * DCHORD

         X1(N) = XXX1 * COST + YYY1 * SINT + XLED
         Y1(N) = -XXX1 * SINT + YYY1 * COST + YLED

      ENDDO

      NDDAT = 2 * NNDD + 1
      

 1000 CONTINUE

      OPEN(9,FILE='ductsect.plt',STATUS='UNKNOWN')

      IF(IDUCT .NE. 0 .AND. IDTHK .NE. -1 .AND. IDCAM .NE. -1) THEN 
         write(9,*) 'zone T="Duct Camber Surface"'
         
         DO N = 1 , NNDDP
            WRITE(9,*) X1(N) , Y1(N)
         ENDDO
      ENDIF

      write(9,*) ' zone T="Duct Section"'

      DO N = 1 , NDDAT
         WRITE(9,*) XID(N) , ETAD(N)
      ENDDO

      ENDFOIL = XID(1)

      write(*,*) 'End of Foil =', endfoil

      RETURN
      END


C ==============================================================
      SUBROUTINE DUCTGEO
C ==============================================================
      INCLUDE 'PUFCAV.INC'
 
      IF(IDGEO. EQ. 0) THEN
        IF(IAN.EQ.6) THEN  ! S.N.KIM made DUCTGEO1 to adapt the duct panels to the blade wake. This subroutine is called with the Full Wake Alignment.
          CALL DUCTGEO1
        ELSE
          CALL DUCTGEO11
        ENDIF
      ELSEIF(IDGEO .EQ. 1 .AND. DUCTPT .EQ. 0.0) THEN
        CALL DUCTGEO3(1)
      ELSEIF(IDGEO .EQ. 1 .AND. DUCTPT .NE. 0.0) THEN
        CALL DUCTGEO2
      ENDIF

      write(9,*) ' zone T="Duct Section2"'

      DO N = 1 , NDUCTP
         RDD = SQRT( YD(N,1)**2 + ZD(N,1)**2)
         WRITE(9,*) XD(N,1) , RDD
      ENDDO

      CLOSE(9)

      RETURN
      END
      
