       SUBROUTINE STSOLN
************************************************************************
*     STSOLN: STeady SOLutioN                                          *
*      --- Solve the steady problem                                    *
*                                                                      *
*   Date of last revision      Revision                                *
*   ---------------------      --------                                *
*   05-24-90 CYHsin  Include correction term of DELP                   *
*   02-27-98 JY      Added a new file "POT.DBG" to check difference    *
*                    in potential at T.E for each strip.               *
*   06-14-99 HL,JY   This subroutine was modified to correctly         *
*                    implement the pressure kutta conditon.            *
*                                                                      *
************************************************************************
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

C       COMMON/BKUTTA/POTEMP(NTZ),WWK(NTZ,MBZ)
C       DIMENSION WPRE(NTZ)
       ALLOCATABLE :: WPRE(:) 


CSH--Calculate the influence coeff from the image, so NBLADE=2
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2
CSH--------------------------------------------     
C-----------------------------------------------------------------------
C     Set up parameters for the iterative matrix solver
C-----------------------------------------------------------------------

C.....FILE 50 [A] + [W] (including far wake) for all blades
C.....FILE 41 [A] for all blades

       IUMAS = 50
       IUMAK = 41

       NORD=NPANEL
       NBLOCK=MR
       IF(IHUB .NE. 0) NBLOCK = NBLOCK + NHBX  ! NHBX = NHBU + NH + NHBDT
       IF(IDUCT .NE. 0) NBLOCK = NBLOCK + MDUCT
       IF(ITUN .NE. 0) NBLOCK = NBLOCK + NAXT
C       IF(IAN .EQ. 2) NBLOCK = NBLOCK + NTHX + NCVX

       DO 10 M=1,MR
          NPERB(M)=NC
 10    CONTINUE    

       NNN = MR

       IF(IHUB .NE. 0) THEN
          DO 20 N=1,NHBX
             NPERB(NNN+N)=MHBT
 20       CONTINUE
          NNN = NNN + NHBX
       ENDIF

       IF(IDUCT .NE. 0) THEN
          DO N=1,MDUCT
             NPERB(NNN+N)=NDUCT
          ENDDO
          NNN = NNN + MDUCT
       ENDIF

       IF(ITUN .NE. 0) THEN
          DO N=1,NAXT
             NPERB(NNN+N)=MTUNEL
          ENDDO
          NNN = NNN + NAXT
       ENDIF

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
c       IF(IAN.EQ.2) THEN
c          DO N=1,NTHX
c             NPERB(NNN+N)=MCVT
c          END DO
c          
c          NNN=NNN+NTHX
c          DO N=1,NCVX
c             NPERB(NNN+N)=MCVT
c          END DO
c          NNN = NNN + NCVX
c       END IF
C/e S.N.KIM | Aug. 2018.

C -- End Tip HSLEE(10/12/99)

C.....Tolerance of the matrix solution is set to be 0.000001............
       TOL=5.0E-06
       ITMAX = 100

C-----------------------------------------------------------------------
C     Solve the simultaneous equation and impose an iterative pressure 
C       Kutta condition
C-----------------------------------------------------------------------
C.....Solving Morino's problem..........................................
C --- POTEMP : Solution (Potential) of Morino's problem
       
CVV
      ALLOCATE(WPRE(NTZ))   
CVV
       IF(ITERKT.NE.0) THEN
          IPASS = 1
          DO M=1,MR
             DO I=1,NORD  ! NORD = NPANEL
                WWK(I,M) = 0.0
                WPRE(I)=-WINF(I,M)
             ENDDO
             CALL NBLIC1(WWK(1,M),WPRE,TOL,NORD,NBLOCK,NPERB,ITMAX
     *            ,IPASS,IUMAK,IHUB,IAN)
             IPASS = 2
          ENDDO
       END IF

       IPASS = 1
       DO N = 1 , NPANEL
          POTEMP(N) = 0.0
       ENDDO

       CALL BLIC1(POTEMP,B,TOL,NORD,NBLOCK,NPERB,ITMAX,IPASS,IUMAS)

       IF(ITERKT.NE.0) THEN
          CALL KUTTA(ITERKT)
       ELSE
          DO N = 1 ,NORD
             POT(N) = POTEMP(N)
          ENDDO
          CALL PRSDIF(1)
       END IF

C -- Begin Tip HSLEE(10/12/99)

c       if(ian.eq.2) call calvel2

          WRITE(1700,*) 
     %         'ZONE T="Steady Wetted Circulation"'
          WRITE(1701,*) 
     %         'ZONE T="Steady Wetted Circulation"'
       
       Gamma = delp(mr)
   
C -- End Tip (10/12/99)

       IF(ICON .NE. 5) THEN
          IF(ICAVT .LE. ICAVMAX) THEN
             DO M=1,MR
                UR=SQRT(1.+(.7*PI/ADVCO)**2.)
                AVCIR = HUNTPI*DELP(M)/UR
                AVCIR1 = HUNTPI*DELP(M)
                WRITE(1700,*) HRZP(1,M),AVCIR
                WRITE(1701,*) HRZP(1,M),AVCIR1
             END DO
          ENDIF
       ELSE
          DO M=1,MR
             AVCIR = HUNTPI*DELP(M)
             WRITE(1700,*) HRZP(1,M),AVCIR
             WRITE(1700,*) HRZP(1,M),AVCIR
          END DO
       ENDIF
          

CSH--REPLACE, NBLADE=1-------------------------
       IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1
CSH--------------------------------------------
 
CVV
      DEALLOCATE(WPRE)
CVV
       RETURN
       END
      
