      SUBROUTINE INDPOT(IPPP)
************************************************************************
*     INDPOT: INDuced POTentials calculation                           *
*      --- Calculate the influence coefficients                        *
*                                                                      *
*   Date of last Revision        Revision                              *
*   ---------------------        --------                              *
*     05-09-89   Fully unsteady scheme                                 *
*     05-14-89   Multi-bladed, fully unsteady scheme                   *
*     05-04-92   NF Added scratch file 53 to house steady wake ic's    *
*     07-26-92   NF Added scratch file 47 to house unsteady dpdn's     *
*     CM010598   Changes made by HSLEE for foil case                   *
*     23-08-18   S.N.KIM modified to incorporate full wake alignment   *
*                on blade/duct wake.                                   *
*                                                                      *
************************************************************************
      use su_inp
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      CHARACTER*50 FN41,FN42,FN43,FN45,FN46,FN47,FN48,FN49,
     %             FN50,FN51,FN52,FN53, FN_FPG
      CHARACTER*2 BLAIDX(10),BLAIDX2(10)
      CHARACTER*60 :: WKCHECK3,WKCHECK4,WKCHECK5,FILENAME3
C      DIMENSION RRHS(NTZ)
      CHARACTER*29 :: INFCHECK,FILENAMEINF
      CHARACTER*29 :: BCHECK,FILENAMEB
      ALLOCATABLE :: RRHS(:)

      DATA BLAIDX /'01','02','03','04','05','06','07','08','09','10'/
      DATA BLAIDX2 /'11','12','13','14','15','16','17','18','19','20'/

CVV
      ALLOCATE(RRHS(NTZ))
CVV

CSH--Calculate the influence coeff from the image, so NBLADE=2
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2
CSH--------------------------------------------
      NWMIN=999
      DO M=1,MRP
         NWMIN=MIN(NWMIN,NSW(M))
      END DO
      NWMIN=NWMIN-1
!YE TIAN add following line for ian=6 04/15/2012
      if ((ian.eq.6).or.(ian.eq.8)) NWMIN = max(NWMIN,NWPANEL)
      NWMIN1=NWMIN-1+INT(DTPROP)*2
      NRECL1=MAX(NPANEL,NWMIN1*MR)*8

C.....Make the record length to be the multipliers of 1024..............
      N1=NRECL1/1024
      NLEFT=NRECL1-N1*1024
      IF(NLEFT.EQ.0) THEN
        NRECL=1024*N1
      ELSE
        NRECL=1024*(N1+1)
      END IF
C      IF(NTSTEP.EQ.0) THEN
c        WRITE(*,'(A)') ' '
c        WRITE(*,'('' ......Direct access file record length: '',I8)')
c     *        NRECL
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
      OPEN(45,FILE=FN45,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      OPEN(46,FILE=FN46,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      OPEN(47,FILE=FN47,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      OPEN(48,FILE=FN48,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      IF(IDUCT .NE. 0) THEN
         FN49=FNSCR(1:LENCH)//'S49.DAT'
         OPEN(49,FILE=FN49,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *                  ACCESS='DIRECT',RECL=NRECL)
      ENDIF

      OPEN(50,FILE=FN50,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(51,FILE=FN51,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(52,FILE=FN52,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(53,FILE=FN53,STATUS='UNKNOWN',FORM='UNFORMATTED')

C.....Open file to record the wake inf. functions.......................

      DO K=1,NBLADE+1
        IO=80+K
        FN43=FNSCR(1:LENCH)//'S'//BLAIDX(K)//'.DAT'
        OPEN(IO,FILE=FN43,STATUS='UNKNOWN',FORM='UNFORMATTED')
      END DO

      DO K=1,NBLADE+1
        IO=500+K
        FN43=FNSCR(1:LENCH)//'S'//BLAIDX2(K)//'.DAT'
        OPEN(IO,FILE=FN43,STATUS='UNKNOWN',FORM='UNFORMATTED')
      ENDDO


C-----------------------------------------------------------------------
C     At IPPP=1 (defined in propcav.f), cal, the influence functions
C-----------------------------------------------------------------------
!s--YE TIAN if ippp >0 then cal. the influence functions
!   ippp = 1, cal all
!   ippp = 2, cal only wake induced
      IF(IPPP .GE. 1) THEN

C.......Calculate control point locations and normal vector.............
        IF(IPPP .eq. 1) CALL CONPT

C.......Calculate the influence coefficients due to the ultimate .......
C-----------------------------------------------------------------------
C     Wake disk and hub disk
C-----------------------------------------------------------------------
        CALL CLEAR(CHDK,NPANEL)
        CALL CLEAR(CHDKDT,NPANEL)
        CALL CLEAR(WDK,NPANEL)
!YE TIAN 04/12/2012
!       IF(IAN.NE.2) THEN
        IF((IAN.NE.2).and.(IAN.NE.6)) THEN
          IF(ICON.NE.5.AND.ICON.NE.6.AND.ISP.NE.1.AND.ICON.NE.8)
     *      CALL INFDISK
          CALL INFWAKS
          CALL INFWAK
c          WRITE(*,*) 'INFWAK is called'
        ELSE
          IF(IHUB .NE. 0) CALL INFDISK2
          CALL INFWAKS
          CALL INFWAK2_FAST
c          WRITE(*,*) 'INFWAK2_FAST is called'
          IF(IDFWA.EQ.1.AND.ippp.NE.1) THEN
            CALL INFWAK_DUCT_FAST
            CALL INFWAKSD
c            WRITE(*,*) 'INFWAK_DUCT_FAST & INFWAKSD are called'
          ENDIF
        END IF

        DO I=1,NPANEL
          DO M=1,MR
            W(I,M)=WSTINF(I,M)
          END DO
        END DO

        DO M=1,MR
          CALL WRITE1(53,W(1,M),NPANEL)
        END DO

C.......Calculate the influence coefficients due to the blades, Hub.....
C.....,Duct, Tunnel Wall(If adopted), and Tip Vortex Cavity(If adopted)

        IF (ippp .EQ. 1) THEN
          IF(IDUCT .NE. 0) THEN
            CALL INFWAK_DUCT_FAST
            CALL INFWAKSD
c            WRITE(*,*) 'INFWAK_DUCT_FAST & INFWAKSD are called'
          ENDIF
          IF(INRPAN.EQ.1) THEN
            CALL INFCOF2
c            WRITE(*,*) 'INFCOF2(NRPAN) is called'
          ELSE
            CALL INFCOF
c            WRITE(*,*) 'INFCOF(RPAN) is called'
          ENDIF
        ENDIF

C.......Sum all the influence coefficients and generate the matrix......
C.......system..........................................................

        CALL INFGEN

C.......WSTINF now is the influence function of the far wake............
C.......(the strength of which is assumed to be the steady solution)....
C.......Read from 80+NBLADE+1...........................................
        IO=81+NBLADE
        REWIND IO

        DO M=MR,1,-1
           CALL READ1(IO,TEMP1,NPANEL)
           DO I=1,NPANEL
              WSTINF(I,M)=TEMP1(I)
           END DO
        END DO
      END IF
C....If ISP=1, this subroutine is only called to calculate the
C....influence coefficients. (JY112499)
      IF(ISP.EQ.1) GO TO 200
CSH--------------------------------------------------------------------
CSH     interpolate to get the inflow over the rudder control points
CSH--------------------------------------------------------------------
      if(itergb .eq. 1) CALL GBPOST(1)
CSH-----------------------------rudder----Shreenaath 05/20/03---------------
C-----------------------------------------------------------------------
C     Strength of the wake dipoles at different time steps
C-----------------------------------------------------------------------

      CALL RHSUPD(RRHS)

      DO I=1,NPANEL
         B(I)=RRHS(I)
      END DO

 200  CONTINUE
CSH--REPLACE, NBLADE=1-------------------------
      IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------

CVV
      DEALLOCATE(RRHS)
CVV
      RETURN
      END




!s--- YE TIAN --- 07/12/2013----
!  for hullfpp
      SUBROUTINE write_fpg
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      CHARACTER*50 fname

      CALL CHRLEN(FN,LENCH2)
      fname =FN(1:LENCH2)//'.fpg'

!     write(*,*) fname
!     write(*,*) NWMIN
      NREAD = NWMIN * MR
      OPEN(20131,FILE=fname,STATUS='REPLACE',FORM='UNFORMATTED')
!     write(*,*) fname
!     OPEN(20131,FILE=fname,STATUS='REPLACE')
      write(20131) NBLADE,NTPREV,MR,MRP,NC,NCP,NRECL,NHBX,MHBT
      write(20131) ADVCO,NTHX,NCVX,MCVT,NTHXP,NCVXP,MCVTP
      write(20131) NPANEL,NREAD,NPWAKS,NTRA,NWSUB,NUWDK,MPDK
      write(20131) RULT,DELK,TANBUW
!     write(*,*) 'NRECL=', NRECL


      write(20131)
     &        ((XB(N,M),N=1,NCP),M=1,MRP),((YB(N,M),N=1,NCP),M=1,MRP),
     &        ((ZB(N,M),N=1,NCP),M=1,MRP),(NSW(M),M=1,MRP),
     &        ((XW(N,M),N=1,NSW(M)),M=1,MRP),
     &        ((YW(N,M),N=1,NSW(M)),M=1,MRP),
     &        ((ZW(N,M),N=1,NSW(M)),M=1,MRP),
c     &        ((XWDK(N,K),N=1,NUWDK+1),K=1,2),
c     &        ((YWDK(N,K),N=1,NUWDK+1),K=1,2),
c     &        ((ZWDK(N,K),N=1,NUWDK+1),K=1,2),
     &        ((XWS(N,M),N=1,NTRA+1),M=1,MRP),
     &        ((YWS(N,M),N=1,NTRA+1),M=1,MRP),
     &        ((ZWS(N,M),N=1,NTRA+1),M=1,MRP)
      CLOSE(20131)

      RETURN
      END SUBROUTINE write_fpg
!e--- YE TIAN --- 07/12/2013----
