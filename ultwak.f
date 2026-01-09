      SUBROUTINE ULTWAK
C***********************************************************************
C     ULTWAK: ULTmate WAKe geometry
C      --- Geometry of the ultimate wake which is replaced by a source
C          disk
C
      INCLUDE 'PUFCAV.INC'
      NUWDK=MHBT
      DTH=DELK/NUWDK
      DO 10 N=1,NUWDK+1
         XWDK(N,1)=XUWDK
         XWDK(N,2)=XUWDK
         TH=DTH*(N-1)
         YWDK(N,1)=RHULT*COS(TH)
         ZWDK(N,1)=RHULT*SIN(TH)
         YWDK(N,2)=RULT*COS(TH)
         ZWDK(N,2)=RULT*SIN(TH)
10    CONTINUE

      RETURN
C))))))))))))))))))))) End of subroutine ULTWAK ((((((((((((((((((((((((
      END
