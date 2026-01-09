      SUBROUTINE OPFILEC
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      CHARACTER*50 FNVOL
C      CHARACTER*50 FNCPH
      CHARACTER*50 FN110,FN111,FN112,FN113,FNWK,FNPRS,FNCAVP
      CHARACTER*2 BLAIDX(10)
      DATA BLAIDX /'01','02','03','04','05','06','07','08','09','10'/

      IF(NTSTEP.EQ.1)THEN

C-----------------------------------------------------------------------
C       open scratch files for unsteady cavity solution
C-----------------------------------------------------------------------

         CALL CHRLEN(FNSCR,LENCH)
         CALL CHRLEN(FN,LENCH)

C.......Let's open some more files (non-scratch files)..................
C.......blade pressure from cavity solution.............................

         FNCAVP=FN(1:LENCH)//'.cavp'
         OPEN(630, FILE=FNCAVP, STATUS='UNKNOWN')
         WRITE(630,1100)
         WRITE(630,6101)

         IF(IHUB.NE.0.AND.IPHUB.EQ.1) THEN

            OPEN(18,FILE='hub-solc-3D.plt',STATUS='UNKNOWN')
            OPEN(310,FILE='hub-prs.plt',STATUS='unknown')

            WRITE(18,2100)
            WRITE(18,6100)

            WRITE(310,2200)
            WRITE(310,6200)

         END IF

         IF(ITUN .NE. 0) THEN
            OPEN(632, FILE='tun-solc-3D.plt', STATUS='UNKNOWN')
            WRITE(632,3100)
            WRITE(632,6100)
         ENDIF

         IF(ISCAV .NE. 2) THEN
            FNPRS=FN(1:LENCH)//'.prs'
            OPEN(680, FILE=FNPRS, STATUS='UNKNOWN')
            WRITE(680,1200)
            WRITE(680,6300)
         ENDIF

         IF(IDUCT .NE. 0) THEN

            OPEN(312, FILE='duct-solc-3D.plt', STATUS='UNKNOWN')
            OPEN(313, FILE='duct-prs.plt', STATUS='UNKNOWN')
            OPEN(314,FILE='duct-mean-prs.plt',STATUS='UNKNOWN')

                          
            WRITE(312,4100)
            WRITE(312,6100)

            WRITE(313,4200)
            WRITE(313,6200)
            
            WRITE(314,4300)
            WRITE(314,6400)

         ENDIF

C.......cavity volume history.........................................

         FNVOL=FN(1:LENCH)//'.vol'
         OPEN(58,FILE=FNVOL,STATUS='UNKNOWN',FORM='FORMATTED')

         WRITE(58,1300)
         IF(ISTEADY.EQ.0) THEN
            WRITE(58,6600)
         ELSE
            WRITE(58,6500)
         END IF

         WRITE(58,6700) 

      END IF

      IF(NTSTEP.EQ.1.OR.(IAN.EQ.2.AND.ISTEADY.NE.0)) THEN

         CALL CHRLEN(FNSCR,LENCH)
         FN110=FNSCR(1:LENCH)//'S110.DAT'
         OPEN(110,FILE=FN110,STATUS='UNKNOWN',FORM='UNFORMATTED')
         FN111=FNSCR(1:LENCH)//'S111.DAT'
         OPEN(111,FILE=FN111,STATUS='UNKNOWN',FORM='UNFORMATTED')
         FN112=FNSCR(1:LENCH)//'S112.DAT'
         OPEN(112,FILE=FN112,STATUS='UNKNOWN',FORM='UNFORMATTED')
         FN113=FNSCR(1:LENCH)//'S113.DAT'
         OPEN(113,FILE=FN113,STATUS='UNKNOWN',FORM='UNFORMATTED')
         
         DO 5 K=1,NBLADE
            IO=90+K
            FNWK=FNSCR(1:LENCH)//'SW1'//BLAIDX(K)//'.DAT'
            OPEN(IO,FILE=FNWK,STATUS='UNKNOWN',FORM='UNFORMATTED')
            IO=30+K
            FNWK=FNSCR(1:LENCH)//'SW2'//BLAIDX(K)//'.DAT'
            OPEN(IO,FILE=FNWK,STATUS='UNKNOWN',FORM='UNFORMATTED')
            IO=520+K
            FNWK=FNSCR(1:LENCH)//'SWD'//BLAIDX(K)//'.DAT'
            OPEN(IO,FILE=FNWK,STATUS='UNKNOWN',FORM='UNFORMATTED')
 5       CONTINUE
         
      ENDIF
      
      
 1100 FORMAT(1x,'TITLE="Cavitating Pressure & Potential on Blade"')
 1200 FORMAT(1x,'TITLE="Cavitating Pressure on Blade (2D)"')
 1300 FORMAT(1x,'TITLE="Cavity Volume vs. Blade Angle"')
      
 2100 FORMAT(1x,'TITLE="Cavitating Pressure & Potential on Hub"')
 2200 FORMAT(1x,'TITLE="Cavitating Pressure on Hub (2D)"')
      
 3100 FORMAT(1x,'TITLE="Cavitating Pressure & Potential on Tunnel"')
      
 4100 FORMAT(1x,'TITLE="Cavitating Pressures & Potential on Duct"')
 4200 FORMAT(1x,'TITLE="Cavitating Pressure on Duct (2D)"')
 4300 FORMAT(1x,'TITLE="Wetted Mean Pressures on duct (2D)"')
      
 6100 FORMAT(1x,'VARIABLES="x","y","z","-Cp","Pot"')
 6101 FORMAT(1x,'VARIABLES="x","y","z","-Cp","Pot","Re_r","Cf"')
 6200 FORMAT(1x,'VARIABLES = "s/s_max","-Cp"')
 6300 FORMAT(1x,'VARIABLES = "Xd/C", "-Cp"')
 6400 FORMAT(1x,'VARIABLES = "x","-Cp(Pm)","-Cp(Potm)"')
 6500 FORMAT(1X,'VARIABLES="Blade Angle(Deg)","Vol-Tot",',
     *     '"Vol-Back","Vol-Face","Vol-Super"')
 6600 FORMAT(1X,'VARIABLES="NTSTEP","Vol-Tot",',
     *     '"Vol-Back","Vol-Face","Vol-Super"')
 6700 FORMAT(1x,'ZONE T="Cavity Volume"')
 
      
      RETURN
      END
      






