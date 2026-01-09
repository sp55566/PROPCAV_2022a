
      SUBROUTINE CBLIC1(SIG,RHSM,TOL,NORD,NBLOCK,NPERB,ITMAX,AINF,NNR)

************************************************************************
*                                                                      *
*  Date of last revision                      Revision                 *
* -----------------------                   -------------              *
************************************************************************
  
      INCLUDE 'PARAM.INC'
      COMMON/INTG/NC,MR,MRTIP,NH,NHP,NHM,NCP,MRP,NPANB
     *     ,      NPANH,NPANW,NPANEL
!     PARAMETER(NNR=NTZ)

      DIMENSION RHSM(*),SIG(*),NPERB(*)
      DIMENSION AINF(NNR,NNR)
C-----------------------------------------------------------------------
C     Call the accelerated matrix solver
C-----------------------------------------------------------------------
      CALL CBLIC2(SIG,RHSM,TOL,NORD,NBLOCK,NPERB,ITMAX,AINF,
     &            NTZ,NBLKMAX,101)

      RETURN
      END




