      SUBROUTINE READXTRAN
C**********************************************************************
C     This subroutine reads the transition location (x/c) at both the *
C     suction and pressure sides for the boundary layer analysis      *
C                                                                     *
C      XTRANLP(M)      M=1,...MR       for pressure side              *  
C      XTRANLS(M)      M=1,...MR       for suction side               *
C                                                                     *
C            Created by Hong Sun        December 2006                 *
C**********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFBL.INC'
       
      OPEN(UNIT=8,FILE='xtran.dat')
      READ(8,*) 
      DO M = 1, MR
         READ(8,*) I1, XTRANLP(M), XTRANLS(M) 
       WRITE(*,*) 'XTRANLS=', XTRANLS(M), 'XTRANLP=', XTRANLP(M)        
      ENDDO

      IF(IDUCT.EQ.1.AND.IDOPT.EQ.1) THEN
        READ(8,*) 
        DO M = 1, MDUCT
           READ(8,*) I1, XTRANDP(M), XTRANDS(M) 
        ENDDO      
      ENDIF 
      
      CLOSE(8)

      RETURN 
      END
