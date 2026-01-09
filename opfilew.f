C     
C     Subroutine to open file for the fully wetted results
C     
C     -----------------------------------------
      SUBROUTINE OPFILEW
C     -----------------------------------------
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      CHARACTER*29 FNCIRW,FNWKTKQ,FNWPRS,FNWETP,FNDET
C      CHARACTER*29 FNCPH
      CHARACTER*29 FNWKTKQ_BL
      CHARACTER*29 FNCPBUSTD,FILENAME3,FNCPBUSTD_BL

      CALL CHRLEN(FN,LENCH)

      IF(ISCAV .EQ. 0) THEN
        FNCIRW=FN(1:LENCH)//'.cirw'
        FNWKTKQ=FN(1:LENCH)//'.wktkq'
        FNWETP=FN(1:LENCH)//'.wetp'
        FNWKTKQ_BL=FN(1:LENCH)//'.wktkq_bl' !Hong

        IF(ICON .NE. 5 .AND. ISTEADY .NE. 0) THEN
          OPEN(15,FILE=FNCIRW,STATUS='UNKNOWN',FORM='FORMATTED')
          WRITE(15,1100)
          WRITE(15,6100)
        ENDIF

        OPEN(73,FILE=FNWKTKQ,STATUS='UNKNOWN')

        WRITE(73,1200)
        WRITE(73,6200)

C..........Added viscous forces output  by Hong Sun, 2006            
        IF(IVISC.EQ.1) THEN
          OPEN(74,FILE=FNWKTKQ_BL,STATUS='UNKNOWN')
          WRITE(74,1200)
          WRITE(74,6200)
        ENDIF

        OPEN(620, FILE=FNWETP, STATUS='UNKNOWN')
        WRITE(620,1300)
        WRITE(620,6301)

        IF(IHUB .NE. 0 .AND. IPHUB .EQ. 1) THEN
          OPEN(17,FILE='hub-solw-3D.plt',STATUS='UNKNOWN')
          WRITE(17,2100)
          WRITE(17,6300)

          OPEN(311,FILE='hub-wprs.plt',STATUS='UNKNOWN')
          WRITE(311,2200)
          WRITE(311,6400)

C ---S.H.CHANG 03/10/2010

          OPEN(317, FILE='blade-velw1-3D.plt', STATUS='UNKNOWN')
          WRITE(317,4260)
          WRITE(317,6350)
          WRITE(317,6800) 0.0,NCP,MRP

          OPEN(318, FILE='blade-velw2-3D.plt', STATUS='UNKNOWN')
          WRITE(318,4260)
          WRITE(318,*) 'VARIABLES="X","Y","Z","U","V","W" '
          WRITE(318,*) 'ZONE T ="VEC" '

          OPEN(319, FILE='hub-velw1-3D.plt', STATUS='UNKNOWN')
          WRITE(319,4240)
          WRITE(319,6350)
          WRITE(319,6800) 0.0,NHBX+1,MHBT+1

          OPEN(320, FILE='hub-velw2-3D.plt', STATUS='UNKNOWN')
          WRITE(320,4240)
          WRITE(320,*) 'VARIABLES="X","Y","Z","U","V","W" '
          WRITE(320,*) 'ZONE T ="VEC" '
C ---S.H.CHANG 03/10/2010

        ENDIF

        IF(ITUN .eq. 1) THEN
          OPEN(621, FILE='tun-solw-3D.plt', STATUS='UNKNOWN')
          WRITE(621,3100)
          WRITE(621,6300)

C ---S.H.CHANG 03/10/2010
          OPEN(321, FILE='tun-velw1-3D.plt', STATUS='UNKNOWN')
          WRITE(321,4270)
          WRITE(321,6350)
          WRITE(321,6800) 0.0,NAXT+1,MTUNEL+1

          OPEN(322, FILE='tun-velw2-3D.plt', STATUS='UNKNOWN')
          WRITE(322,4270)
          WRITE(322,*) 'VARIABLES="X","Y","Z","U","V","W" '
          WRITE(322,*) 'ZONE T ="VEC" '
C ---S.H.CHANG 03/10/2010

        ENDIF

        IF(IDUCT .NE. 0) THEN

          OPEN(313,FILE='duct-wprs.plt', STATUS='UNKNOWN')            
          WRITE(313,4100)
          WRITE(313,6400)

          OPEN(312, FILE='duct-solw-3D.plt', STATUS='UNKNOWN')
          WRITE(312,4200)
          WRITE(312,6300)

          OPEN(314,FILE='duct-mean-wprs.plt',STATUS='UNKNOWN')
          WRITE(314,4300)
          WRITE(314,6500)

C ---S.H.CHANG 03/10/2010
          OPEN(315, FILE='duct-velw1-3D.plt', STATUS='UNKNOWN')
          WRITE(315,4250)
C          WRITE(315,6350)
          WRITE(315,6351)
C          WRITE(315,6800) 0.0,NDUCTP,MDUCTP
          WRITE(315,6801) 0.0,NDUCTP,MDUCTP

          OPEN(316, FILE='duct-velw2-3D.plt', STATUS='UNKNOWN')
          WRITE(316,4250)
          WRITE(316,*) 'VARIABLES="X","Y","Z","U","V","W" '
          WRITE(316,*) 'ZONE T ="VEC" '
C ---S.H.CHANG 03/10/2010

          OPEN(622,FILE='duct-wktkq.plt',STATUS='UNKNOWN')
          WRITE(622,4400)
          WRITE(622,6200)

C..........Added viscous forces output  by Hong Sun, 2006            
          IF(IVISC.EQ.1.AND.IDOPT.EQ.1) THEN
            OPEN(623,FILE='duct-wktkq_bl.plt',STATUS='UNKNOWN')
            WRITE(623,4400)
            WRITE(623,6200)
          ENDIF            
          
        ENDIF

      ENDIF
      
      FNWPRS=FN(1:LENCH)//'.wprs.plt'
      OPEN(690, FILE=FNWPRS, STATUS='UNKNOWN')
!Allen Du 12/22/2017 output cp at the control points
      OPEN(7042, FILE="control_wprs.plt", STATUS='UNKNOWN')
!s---Allen Du 01/11/2018 output the cp_min
      OPEN(7043, FILE="cpmin.plt", STATUS='UNKNOWN')
      WRITE(7043,7044)
      WRITE(7043,7045)
 7044 FORMAT(1x,'VARIABLES = "T","Cp_min"')
 7045 FORMAT(1x,'ZONE T="fully wetted"')
!e---Allen Du 01/11/2018 output the cp_min

      IF(ISCAV .EQ. 0) THEN
        WRITE(690,1400)
        WRITE(690,6700)
!s---Allen Du 12/22/2017 output cp at the control points
        WRITE(7042,1400)
        WRITE(7042,6701)
!e---Allen Du 12/22/2017 output cp at the control points
      ENDIF
      
C/S.N.KIM - Cp Blade Plotting for unsteady runs.
      IF (IAN.EQ.2.AND.IVISC.EQ.1.AND.ISTEADY.NE.0) THEN
        DO I = 1, MR
          WRITE(FILENAME3,'(F6.5)') HRZP(1,I) 
          FNCPBUSTD='uns-wprs-r0'//trim(FILENAME3)//'.plt'
          FNCPBUSTD_BL='uns-bl-r0'//trim(FILENAME3)//'.plt'
          OPEN(1999+I, FILE=FNCPBUSTD, STATUS='UNKNOWN')
          OPEN(2999+I, FILE=FNCPBUSTD_BL, STATUS='UNKNOWN')
          IF(ISCAV .EQ. 0) THEN
            WRITE(1999+I,1400)
            WRITE(1999+I,6700)
            WRITE(2999+I,1600)
            WRITE(2999+I,6900)
          ENDIF
        ENDDO
      ENDIF
C/S.N.KIM - Cp Blade Plotting for usnteady runs.

      FNDET=FN(1:LENCH)//'.det'
      OPEN(720, FILE=FNDET, STATUS='UNKNOWN')      
      WRITE(720,1500)
      WRITE(720,6600)
      
C--   Title for Blade

 1100 FORMAT(1x,'TITLE = "Wetted Circulation distribution"')
 1200 FORMAT(1X,'TITLE="Total Wetted KT & KQ "')
C     1200      FORMAT(1X,'TITLE="Wetted KT & KQ per blade"') Hong temporary change
 1300 FORMAT(1x,'TITLE="Wetted Pressure & Potential on Blade"')
 1400 FORMAT(1x,'TITLE="Wetted Pressures on Blade (2D)"')
 1500 FORMAT(1x,'TITLE = "Detachment Locations"')
 1600 FORMAT(1x,'TITLE="Wetted Pressures on Blade BLayer (2D)"')

C     -- Title for HUB

 2100 FORMAT(1x,'TITLE="Wetted Pressure & Potential on Hub"')
 2200 FORMAT(1x,'TITLE="Wetted Pressures on Hub (2D)"')

C     -- Title for TUNNEL

 3100 FORMAT(1x,'TITLE="Wetted Pressure & Potential on Tunnel"')

C     -- Title for DUCT

 4100 FORMAT(1x,'TITLE="Wetted Pressures on duct (2D)"')
 4200 FORMAT(1x,'TITLE="Wetted Pressure & Potential on Duct"')

 4240 FORMAT(1x,'TITLE="Wetted Velocity,Pressure&Potential on Hub"')
 4250 FORMAT(1x,'TITLE="Wetted Velocity,Pressure&Potential on Duct"')
 4260 FORMAT(1x,'TITLE="Wetted Velocity,Pressure&Potential on Blade"')
 4270 FORMAT(1x,'TITLE="Wetted Velocity,Pressure&Potential on Tunnel"')

 4300 FORMAT(1x,'TITLE="Wetted Mean Pressures on duct (2D)"')
 4400 FORMAT(1X,'TITLE="Wetted Total Force on Duct"')

C     -- Variables

 6100 FORMAT(1x,'VARIABLES = "Angle","100*DELP/(2*PI*R*UR)"')
 6200 FORMAT(1X,'VARIABLES="ANGLE","KT","KQ","KTV","KQV"')
 6300 FORMAT(1x,'VARIABLES="X","Y","Z","-Cp","POT"')
 6301 FORMAT(1x,'VARIABLES="X","Y","Z","-Cp","POT","DPDT","Re_r","Cf"')

 6350 FORMAT(1x,'VARIABLES="X","Y","Z","Q","-Cp","POT"')
 6351 FORMAT(1x,'VARIABLES="X","Y","Z","Q","-Cp","POT","DPDU","DPDV"')

 6400 FORMAT(1x,'VARIABLES = "s/s_max","-Cp"')
 6500 FORMAT(1x,'VARIABLES = "x","-Cp(Pm)","-Cp(Potm)"')
 6600 FORMAT(1x,'VARIABLES = "Strip","Det-back","Cav Lgth","No-CavB"',
     *  ',"Det-face","Cav Lgth","No-CavF"')
 6700 FORMAT(1x,'VARIABLES = "Xd/C", "-Cp"')
 6701 FORMAT(1x,'VARIABLES = "x", "-cp"')
 6900 FORMAT(1x,'VARIABLES = "Xd/C", "-Cp"')

 6800 FORMAT(1x,'ZONE T="Time=',F8.1,'", I=',I5,', J=',I5,
     *', K=1, DATAPACKING=BLOCK, VARLOCATION=(4=CELLCENTERED',
     *',5=CELLCENTERED,6=CELLCENTERED) ')

 6801 FORMAT(1x,'ZONE T="Time=',F8.1,'", I=',I5,', J=',I5,
     *', K=1, DATAPACKING=BLOCK, VARLOCATION=(4=CELLCENTERED',
     *',5=CELLCENTERED,6=CELLCENTERED,7=CELLCENTERED',
     *',8=CELLCENTERED) ')

      RETURN
      END

