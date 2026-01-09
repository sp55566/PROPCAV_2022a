      SUBROUTINE USSOLN
************************************************************************
*     USSOLN: UnSteady SOLutioN                                        *
*      --- Solve the unsteady problem at the current timestep          *
*                                                                      *
*  Date of last revision       Revision                                *
*  --------------------        --------                                *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C      COMMON/BKUTTA/POTEMP(NTZ),WWK(NTZ,MBZ)
C      DIMENSION WPRE(NTZ)

      ALLOCATABLE :: WPRE(:)

C-----------------------------------------------------------------------
C     Set up parameters for the iterative matrix solver
C-----------------------------------------------------------------------

C.....FILE 51 [A] + 0.5*[W] (NOT including far wake) for key blade at
C                present time step.
C.....FILE 52 [A] for key blade

      IUMAU=51
      IUMAK=52

      IF(IAN .NE. 2) THEN
         IF(NTSTEP.LE.1) THEN
            IPASS=1
         ELSE
            IPASS=2
         END IF
      ELSE
         IPASS = 1
      ENDIF

      ITMAX=50
      NORD=NPANEL

      NBLOCK = MR
      IF(IHUB .NE. 0) NBLOCK=NBLOCK+NHBX
      IF(IDUCT .NE. 0) NBLOCK = NBLOCK + MDUCT
      IF(ITUN .NE. 0) NBLOCK = NBLOCK + NAXT
c      IF(IAN .EQ. 2) NBLOCK = NBLOCK+NTHX+NCVX

      DO M=1,MR
        NPERB(M)=NC
      ENDDO
      NNN = MR

      IF(IHUB .NE. 0) THEN
         DO N=1,NHBX
            NPERB(NNN+N)=MHBT
         ENDDO
         NNN = NNN + NHBX
      ENDIF

      IF(IDUCT .NE. 0) THEN
         DO N=1,MDUCT
            NPERB(NNN+N)=NDUCT
         ENDDO
         NNN = NNN + MDUCT
      ENDIF

      IF(ITUN .NE. 0) THEN
         DO N = 1 , NAXT
            NPERB(NNN+N) = MTUNEL
         ENDDO
         NNN = NNN + NAXT
      ENDIF

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
c      IF(IAN.EQ.2) THEN
c         DO N=1,NTHX
c            NPERB(NNN+N)=MCVT
c         END DO
c
c         NNN=NNN+NTHX
c         DO N=1,NCVX
c            NPERB(NNN+N)=MCVT
c         END DO
c         NNN = NNN + NCVX
c      END IF
C/e S.N.KIM | Aug. 2018.

C.....Tolerance of the matrix solution is set to be 0.000005............

      TOL=5.0E-06

C-----------------------------------------------------------------------
C     Solving the base problems at time step 1
C-----------------------------------------------------------------------
c.....Only solve the base problem if iterkt.ne.0

CVV
      ALLOCATE(WPRE(NTZ))
CVV

      IFLAG=1
      IF((NTSTEP.EQ.1) .AND. (ITERKT.NE.0)) IFLAG=2
      IF((IAN .EQ. 2) .AND. (ITERKT.NE.0)) IFLAG=2

      IF(IFLAG .EQ. 2) THEN

        IPASS=1

        DO 20 M=1,MR
          DO 10 I=1,NORD
            WWK(I,M)=0.0 
            WPRE(I)=-HALF*W(I,M)-WSUBIF(I,M)
   10     CONTINUE

          CALL BLIC1(WWK(1,M),WPRE,TOL,NORD,NBLOCK,NPERB,ITMAX
     *              ,IPASS,IUMAK)

          IPASS=2
   20   CONTINUE

      END IF

CVV
      DEALLOCATE(WPRE)
CVV

C-----------------------------------------------------------------------
C     Solve the simultaneous equation and impose an iterative pressure 
C       Kutta condition
C-----------------------------------------------------------------------

C.....Solving Morino's problem..........................................
      IPASS=1

      DO 77 I=1,NORD
         POTEMP(I)=0.0
 77   CONTINUE

      CALL BLIC1(POTEMP,B,TOL,NORD,NBLOCK,NPERB,ITMAX,IPASS,IUMAU)

      IF(ITERKT.NE.0) THEN
         CALL KUTTA(ITERKT)
      ELSE
          DO I = 1 , NORD
             POT(I) = POTEMP(I)
          ENDDO
          CALL PRSDIF(1)
      END IF

c      IF(IAN .EQ. 2) CALL CALVEL2

      IF(IDXREV.GT.0.AND.ISTEADY.NE.0.AND.ISTEADY.NE.1.
     *     AND.ISC.NE.1) CALL DPOTDT

C/s S.N.KIM | Unsteady Wake Alignment
      IF((IAN.EQ.2.and.ICAVT.EQ.ICAVMAX).or.(IAN.NE.2))THEN
        WRITE(1700,*)
     *         'ZONE T="Unsteady Wetted Circulation"'
        WRITE(1701,*)
     *         'ZONE T="Unsteady Wetted Circulation"'
        DO M=1,MR
          UR = SQRT(1.+(.7*PI/ADVCO)**2.)
          AVCIR = HUNTPI*DELP(M)/UR
          AVCIR1 = HUNTPI*DELP(M)
            WRITE(1700,*) HRZP(1,M),AVCIR
            WRITE(1701,*) HRZP(1,M),AVCIR1
        ENDDO
      ENDIF
C/e S.N.KIM | Aug. 2018.

      RETURN
      END                                                               




