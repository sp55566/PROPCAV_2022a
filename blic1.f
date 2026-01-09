      SUBROUTINE BLIC1(SIG,RHSM,TOL,NORD,NBLOCK,NPERB,ITMAX,
     * IPASS,IUMAT)

************************************************************************
*     *
*     Date of last revision                      Revision                 *
*     -----------------------                   -------------              *
************************************************************************
      USE MEMSOL  
      INCLUDE 'PARAM.INC'
      COMMON/INTG/NC,MR,MRTIP,NH,NHP,NHM,NCP,MRP,NPANB
     * ,      NPANH,NPANW,NPANEL
      
C     PARAMETER(NNR=NTZ)
C     COMMON /MEMSOL/AINF(NNR,NNR)
      DIMENSION RHSM(*),SIG(*),NPERB(*)

C-----------------------------------------------------------------------
C     Read the influence coeff. file
C-----------------------------------------------------------------------
      IF(IPASS.LE.1) THEN
       REWIND IUMAT
       DO 20 J=1,NORD
        CALL READ1(IUMAT,ALHS(1,J),NORD)
 20    CONTINUE
      END IF

C-----------------------------------------------------------------------
C     This here is to check the calculatin of the influence coeff's.
C     The sum of the dipoles over each panel should be 4*PI.    
C     JY010899
C-----------------------------------------------------------------------

      DO I=1,NORD
       SUM=0.         
       DO J=1,NORD
        SUM=SUM+ALHS(I,J)
c     if(I .eq. NORD-1) write(999,*) I, J, ALHS(I,J),SUM
       END DO
c     WRITE(999,*) I,SUM,rhsm(i)
      END DO

C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Call the accelerated matrix solver
C-----------------------------------------------------------------------
      CALL BLIC2(SIG,RHSM,TOL,NORD,NBLOCK,NPERB,ITMAX,IPASS,
     &           NTZ,NBLKMAX,101)

      RETURN
      END




