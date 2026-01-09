      SUBROUTINE GWAKE_FOIL
************************************************************************
*     GWAKE: Geometry of the WAKE for Foil                             *
*      --- Generate geometries of the transition wake and the ultimate *
*          wake                                                        *
*                                                                      *
************************************************************************

      use m_WAKGEO
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
!     PARAMETER(NWR=MBPZ)
      integer nwr
      nwr = MBPZ
!     COMMON /WAKGEO/ RTBL(11,NWR),XTBL(11),RW(NWR),VAR(NWR),VTR(NWR),
!    *     UASTAR(NWR),UTSTAR(NWR),UAUSTR(NWR),UTUSTR(NWR)

C-----------------------------------------------------------------------
C     Parameters of the transition wake
C-----------------------------------------------------------------------

      XFINAL=XULT

      CCMAX = XB(1,1)
      DO M = 1 , MRP
         CCMAX = AMAX1(CCMAX,XB(1,M))         
      ENDDO

      XXFIN = CCMAX + XFINAL

      NWK=NWPZ-20
      IF(NWK .LT. 60) NWK = 60

      DTH = HALF * PI / (NWK - 1)

      DO 60 M=1,MRP

         XW(1,M)=XB(1,M)
         YW(1,M)=YB(1,M)
         ZW(1,M)=ZB(1,M)

         XSLNG = XXFIN - XW(1,M)

         DO 40 N=2,NWK
            SSS = 1. - COS(DTH*(N-1))
            XW(N,M) = XW(1,M) + SSS * XSLNG 
            ZW(N,M)=ZW(N-1,M)
            YW(N,M)=YW(N-1,M)
 40      CONTINUE

         NSW(M) = NWK

60    CONTINUE

      DO 11 M=1,MR
         IF(NSW(M).EQ.NSW(M+1))THEN
            NWPAN(M)=NSW(M)-1
         ELSE
            NWPAN(M)=MIN(NSW(M),NSW(M+1))
            NWPAN(M)=NWPAN(M)-1
         ENDIF
 11   CONTINUE


      CALL ULTWAK

      RETURN
      END


