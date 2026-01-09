C =======================================
      SUBROUTINE INDPOT2
C
C     This subroutine calculates the induced potentials for
C     the cavitating analysis in the case of 
C     full wake alignment.
C =======================================
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      CHARACTER*50 FN41,FN42,FN43,FN45,FN46,FN47,FN48,FN50,FN51,FN52
     *     ,       FN53
      CHARACTER*2 BLAIDX(10)
      DATA BLAIDX /'01','02','03','04','05','06','07','08','09','10'/

      NWMIN=999
      DO M=1,MRP
         NWMIN=MIN(NWMIN,NSW(M))
      END DO
      NWMIN=NWMIN-1
!S.N.KIM add following line for ian=6 & 8 Aug. 2018.
      if ((ian.eq.6).or.(ian.eq.8)) NWMIN = max(NWMIN,NWPANEL)
      NWMIN1=NWMIN-1+INT(DTPROP)*2
      NRECL1=MAX(NPANEL,NWMIN1*MR)*8

      N1=NRECL1/1024
      NLEFT=NRECL1-N1*1024
      IF(NLEFT.EQ.0) THEN
        NRECL=1024*N1
      ELSE
        NRECL=1024*(N1+1)
      END IF
C      IF(NTSTEP.EQ.0) THEN
C        WRITE(*,'(A)') ' ' 
C        WRITE(*,'('' ......Direct access file record length: '',I8)')
C     *        NRECL
C      END IF

C-----------------------------------------------------------------------
C     Generate file names and open them for solving the matrix 
C-----------------------------------------------------------------------
      CALL CHRLEN(FNSCR,LENCH)

      FN41=FNSCR(1:LENCH)//'S41.DAT'
      FN42=FNSCR(1:LENCH)//'S42.DAT'
      FN45=FNSCR(1:LENCH)//'S45.DAT'
      FN46=FNSCR(1:LENCH)//'S46.DAT'
      FN47=FNSCR(1:LENCH)//'S47.DAT'
      FN48=FNSCR(1:LENCH)//'S48.DAT'
      FN50=FNSCR(1:LENCH)//'S50.DAT'
      FN51=FNSCR(1:LENCH)//'S51.DAT'
      FN52=FNSCR(1:LENCH)//'S52.DAT'
      FN53=FNSCR(1:LENCH)//'S53.DAT'

      OPEN(41,FILE=FN41,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(42,FILE=FN42,STATUS='UNKNOWN',FORM='UNFORMATTED')

      if(iscav .ne. 0) then
      OPEN(45,FILE=FN45,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      OPEN(46,FILE=FN46,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      OPEN(47,FILE=FN47,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      OPEN(48,FILE=FN48,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      endif

      OPEN(50,FILE=FN50,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(51,FILE=FN51,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(52,FILE=FN52,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(53,FILE=FN53,STATUS='UNKNOWN',FORM='UNFORMATTED')
CSH--Calculate the influence coeff from the image, so NBLADE=2
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2
CSH--------------------------------------------           
      DO K=1,NBLADE+1
        IO=80+K
        FN43=FNSCR(1:LENCH)//'S'//BLAIDX(K)//'.DAT'
        OPEN(IO,FILE=FN43,STATUS='UNKNOWN',FORM='UNFORMATTED')
      ENDDO

      CALL CONPT
      CALL CLEAR(CHDK,NPANEL)
      CALL CLEAR(CHDKDT,NPANEL)
      CALL CLEAR(WDK,NPANEL)
      IF(IHUB .NE. 0) CALL INFDISK2
      CALL INFWAKS
      CALL INFWAK2_FAST

      DO  I=1,NPANEL
        DO  M=1,MR
          W(I,M)=WSTINF(I,M)
        ENDDO
      ENDDO

      DO M=1,MR
        CALL WRITE1(53,W(1,M),NPANEL)
      ENDDO

      CALL INFCOF

C-----------------------------------------------------------------------
C     Initialization
C-----------------------------------------------------------------------

C.....FILE 50 (IUMAS) -- [A] for the steady problem.....................
C.....FILE 51 (IUMAU) -- [A] for the unsteady problem
C.....FILE 52 (IUMAK) -- [A] for the p.k. iterations

      IUMAS = 50
      IUMAU = 51
      IUMAK = 52

      REWIND 41
      REWIND 81
      DO M=MR,1,-1

C.......TEMP2 here is wake influence coeff. at L=1, KK=1................
        CALL READ1(81,TEMP2,NPANEL)
        DO 90 N=1,NC
           DO 10 I=1,NPANEL
              A(I)=ZERO
              B(I)=ZERO
 10        CONTINUE
C-----------------------------------------------------------------------
C         KEY blade's influence coefficients
C-----------------------------------------------------------------------
           CALL READ1(41,TEMP1,NPANEL)

           SUMINF=0
           DO 20 I=1,NPANEL
              A(I)=A(I)+TEMP1(I)
 20        CONTINUE

C.........Write matrix [A] without Morino's Kutta condition.............
          CALL WRITE1(IUMAK,TEMP1,NPANEL)

C.........Create [A] for the unsteady problem (with linear correction)..
          IF(N.EQ.1)THEN
             DO 25 I=1,NPANEL  
                TEMP1(I)=TEMP1(I)-HALF*TEMP2(I)-WSUBIF(I,M)
 25          CONTINUE
          ELSEIF(N.EQ.NC)THEN
             DO 30 I=1,NPANEL
                TEMP1(I)=TEMP1(I)+HALF*TEMP2(I)+WSUBIF(I,M)
 30          CONTINUE
          END IF
c.........Tip seperation will be considered later
          CALL WRITE1(IUMAU,TEMP1,NPANEL)

C-----------------------------------------------------------------------
C         OTHER blades' influence coefficients
C----------------------------------------------------------------------- 
          DO 50 KK=2,NBLADE
             CALL READ1(41,TEMP1,NPANEL)
             DO 40 I=1,NPANEL
                A(I)=A(I)+TEMP1(I)
 40          CONTINUE
 50       CONTINUE

C.........Matrix for Morino's problem...................................
C.........for steady case (matrix is written on file 50, IUMAS).........
          DO 60 I=1,NPANEL
C...........Add the influence of the wake...............................
             IF(N.EQ.1) THEN
                A(I)=A(I)-WSTINF(I,M)
             ELSE IF(N.EQ.NC) THEN
                A(I)=A(I)+WSTINF(I,M)
             END IF
 60       CONTINUE

          CALL WRITE1(IUMAS,A,NPANEL)
 90    CONTINUE
      ENDDO

C-----------------------------------------------------------------------
C     Compute influence coefficients due to the hub
C-----------------------------------------------------------------------
      IF(IHUB.NE.0) THEN
        DO 200 N=1,NHBX
          DO 190 M=1,MHBT
            DO 110 I=1,NPANEL
              A(I)=ZERO
              B(I)=ZERO
  110       CONTINUE   

C-----------------------------------------------------------------------
C           KEY blade's influence coefficients
C-----------------------------------------------------------------------
            CALL READ1(41,TEMP1,NPANEL)
            DO 120 I=1,NPANEL
              A(I)=A(I)+TEMP1(I)
  120       CONTINUE

C-----------------------------------------------------------------------
C     Add the influence of the hub disk to IUMAK & IUMAU.  This is
C     necessary because otherwise the effect of the hub disk are not
C     included in the wetted or cavitating solution.
C     S.N.KIM | Dec. 2018.
C-----------------------------------------------------------------------
            IF(N.EQ.1) THEN
               DO I=1,NPANEL
                  TEMP1(I)=TEMP1(I)+CHDK(I)/MHBT
               END DO
            ELSE IF(N.EQ.NHBX) THEN
               DO I=1,NPANEL
                  TEMP1(I)=TEMP1(I)+CHDKDT(I)/MHBT
               END DO
            END IF

C...........Write matrix [A] without Morino's Kutta condition...........
            CALL WRITE1(IUMAK,TEMP1,NPANEL)
            CALL WRITE1(IUMAU,TEMP1,NPANEL)
C-----------------------------------------------------------------------
C           OTHER blades' influence coefficients
C-----------------------------------------------------------------------
            DO 150 KK=2,NBLADE
              CALL READ1(41,TEMP1,NPANEL)
              DO 140 I=1,NPANEL
                A(I)=A(I)+TEMP1(I)
  140         CONTINUE
  150       CONTINUE

CSH--ADDed the influence coeff from the image, 
C...........Add the influence of the hub disk...........................
            IF(N.EQ.1) THEN
               DO 160 I=1,NPANEL
                  A(I)=A(I)+CHDK(I)/MHBT
 160           CONTINUE
C/s S.N.KIM | Dec. 2018.             
            ELSE IF(N.EQ.NHBX) THEN
               DO I=1,NPANEL
                  A(I)=A(I)+CHDKDT(I)/MHBT
               END DO
C/e S.N.KIM | Dec. 2018.
            END IF 
            CALL WRITE1(IUMAS,A,NPANEL)
 190     CONTINUE
 200  CONTINUE
      END IF
      
C*********************************************************************
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
C*********************************************************************
cC -- Begin Tip HSLEE(10/12/99)
c
c           do n = 1 , nthx
c              do m = 1 , mcvt
cc             Hong change the initialzation of A & B        09/13/07              
cc                 a = 0.0          
cc                 b = 0.0
c              DO I=1,NPANEL
c                 A(I)=ZERO
c                 B(I)=ZERO
c              ENDDO  
c
cC --- Read tip hub influence coeff
c      
c               call read1(41,temp1,npanel)
c      
c               do i = 1, npanel
c                   a(i) = a(i) + temp1(i)
c               enddo
c
c               call write1(iumak,temp1,npanel)
c               call write1(iumau,temp1,npanel)
c
c                 call write1(iumas,a,npanel)
c
c             enddo
c          enddo
c
c           do n = 1 , ncvx
c              do m = 1 , mcvt
cc             Hong change the initialzation of A & B        09/13/07              
cc                 a = 0.0          
cc                 b = 0.0
c              DO I=1,NPANEL
c                 A(I)=ZERO
c                 B(I)=ZERO
c              ENDDO  
c
cC --- Read tip hub influence coeff
c
c                 call read1(41,temp1,npanel)
c               do i = 1 , npanel
c                 a(i) = a(i) + temp1(i)
c               enddo
c              
c               call write1(iumak,temp1,npanel)
c               call write1(iumau,temp1,npanel)
c                 call write1(iumas,a,npanel)
c             enddo
c          enddo
C**********************************************************************
C/e S.N.KIm | Aug. 2018.
C**********************************************************************

CSH--REPLACE, NBLADE=1-------------------------
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------

        IO=81+NBLADE
        REWIND IO

        DO M=MR,1,-1
           CALL READ1(IO,TEMP1,NPANEL)
           DO I=1,NPANEL 
              WSTINF(I,M)=TEMP1(I)
           ENDDO
        ENDDO

      RETURN
      END
