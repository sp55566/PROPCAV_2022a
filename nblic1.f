      SUBROUTINE NBLIC1(SIG,RHSM,TOL,NORD,NBLOCK,NPERB,ITMAX,
     *                IPASS,IUMAT,IHUB,IAN)

************************************************************************
*                                                                      *
*  Date of last revision                      Revision                 *
* -----------------------                   -------------              *
************************************************************************
      USE MEMSOL  
      INCLUDE 'PARAM.INC'
!     PARAMETER (NNR=NTZ)
      COMMON/INTG/ NC,MR,MRTIP,NH,NHP,NHM,NCP,MRP,NPANB
     *     ,      NPANH,NPANW,NPANEL
      COMMON/HUBI/NHBU,NHBDT,NHBX,MHBT
      COMMON /MTIP1/  NTHX, NCVX, MCVT,NTHXP,NCVXP,MCVTP
      COMMON /DAT2/ NX, NBLADE
      COMMON /ITUNCON1/ ITUNGEO,ITUN
      COMMON /ITUNCON2/ NAXT, MTUNEL, NPANTN 
      COMMON /IDUCTCON1/ IDUCT
C.....Hong corrected IDUCTCON2 block    07/05/2007
      COMMON /IDUCTCON2/ NDUCT, NDUCTH, NDUCTHP, MDUCT, 
     %                     NDUCTP, MDUCTP,NPAND
C      COMMON /MEMSOL/ AINF(NNR,NNR)
      DIMENSION RHSM(*),SIG(*),NPERB(*)

C      dimension temp(ntz)
      ALLOCATABLE :: TEMP(:)

      integer NNR
      NNR = NTZ


CVV   
      ALLOCATE(TEMP(NTZ))
CVV

C-----------------------------------------------------------------------
C     Read the influence coeff. file
C-----------------------------------------------------------------------
      IF(IPASS.LE.1) THEN

      REWIND IUMAT 

        DO JJ=1,NORD
           DO I=1,NORD
              ALHS(I,JJ)=0.0
           END DO
        END DO

      do m = mr, 1, -1
        do n = 1 , nc
          jj = indexb(n,m)
          do kk= 1 , nblade
             call read1(iumat,temp,nord)
             do i = 1 , nord
               alhs(i,jj)=alhs(i,jj) + temp(i)
             enddo
          enddo
        enddo
      enddo

      if(ihub.ne.0) then       
        do n = 1 , nhbx
          do m = 1, mhbt 
            jj = indexh(n,m)
            do kk= 1 , nblade
               call read1(iumat,temp,nord)
               do i = 1 , nord
                 alhs(i,jj)=alhs(i,jj) + temp(i)
               enddo
            enddo
          enddo
        enddo
      endif

      if(iduct.ne.0) then      
        do m = 1 , mduct
          do n = 1, nduct 
            jj = indexd(n,m)
            do kk= 1 , nblade
               call read1(iumat,temp,nord)
               do i = 1 , nord
                 alhs(i,jj)=alhs(i,jj) + temp(i)
               enddo
            enddo
          enddo
        enddo
      endif

      if(itun.ne.0) then       
        do n = 1 , naxt
          do m = 1, mtunel 
            jj = indextn(n,m)
            do kk= 1 , nblade
               call read1(iumat,temp,nord)
               do i = 1 , nord
                 alhs(i,jj)=alhs(i,jj) + temp(i)
               enddo
            enddo
          enddo
        enddo
      endif

C**********************************************************************
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
C**********************************************************************
c      if(ian .eq. 2) then
c        do n = 1 , nthx
c          do m = 1, mcvt 
c            jj = indext(n,m)
c            call read1(iumat,temp,nord)
c            do i = 1 , nord
c              alhs(i,jj)=alhs(i,jj) + temp(i)
c            enddo
c          enddo
c        enddo
c
c        do n = 1 , ncvx
c          do m = 1, mcvt 
c            jj = indexc(n,m)
c            call read1(iumat,temp,nord)
c            do i = 1 , nord
c              alhs(i,jj)=alhs(i,jj) + temp(i)
c            enddo
c          enddo
c        enddo
c
c      endif
C**********************************************************************
C/e S.N.KIM | Aug. 2018.
C**********************************************************************
      END IF

C-----------------------------------------------------------------------
C     This here is to check the calculatin of the influence coeff's.
C     The sum of the dipoles over each panel should be 4*PI.    
C                                                             JY010899
C-----------------------------------------------------------------------
c      DO I=1,NORD
c         SUM=0.
c         DO J=1,NORD
c            SUM=SUM+ALHS(I,J)
c         END DO
c         WRITE(99,*) ,I,SUM
c      END DO
c      STOP
C-----------------------------------------------------------------------
CVV
      DEALLOCATE(TEMP)
CVV       
C-----------------------------------------------------------------------
C     Call the accelerated matrix solver
C-----------------------------------------------------------------------

      CALL BLIC2(SIG,RHSM,TOL,NORD,NBLOCK,NPERB,ITMAX,IPASS,
     &           NTZ,NBLKMAX,101)

      RETURN
      END




