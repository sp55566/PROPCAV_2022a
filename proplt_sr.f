      SUBROUTINE PROPLT_SR   
************************************************************************
*     Plots propeller geometry for supercavitating propellers.
************************************************************************
      
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION YB1(NBPZ,MBPZ),ZB1(NBPZ,MBPZ)

C-----------------------------------------------------------------------
C     Plot the the actual blades without the separated region 
C-----------------------------------------------------------------------
 7000 Format(1x,'ZONE T="PROP-SR", N=',I5,', E=',I5,' F=FEPOINT')
 7001 FORMAT(3(1x,F16.9))
 7002 FORMAT(4(1X,I6))

      DELK=TWOPI/NBLADE

      NC0=NCOLD+1
      NC1=NC0+1
      WRITE(700,7000) MRP*NC1*NBLADE,MR*NC0*NBLADE

      DO K=1,NBLADE
         T=-DELK*FLOAT(K-1)

         DO M=1,MRP
            DO N=N0(2),N0(1)

               XX0 = XB(N,M)
               YY0 = YB(N,M)
               ZZ0 = ZB(N,M)
               IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)                  
               YB1(N,M)=YY0*COS(T)-ZB(N,M)*SIN(T)
               ZB1(N,M)=YY0*SIN(T)+ZB(N,M)*COS(T)
               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YB1(N,M))                  
               WRITE(700,7001) XX0,YB1(N,M),ZB1(N,M)
            END DO
            N=N0(2)
            XX0 = XB(N,M)
            YY0 = YB(N,M)
            ZZ0 = ZB(N,M)
            IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)                  
            YB1(N,M)=YY0*COS(T)-ZB(N,M)*SIN(T)
            ZB1(N,M)=YY0*SIN(T)+ZB(N,M)*COS(T)
            IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YB1(N,M))                  
            WRITE(700,7001) XX0,YB1(N,M),ZB1(N,M)
         END DO
      END DO

      NSTB=0
      DO IK=1,NBLADE
         DO M = 1 , MR
            DO N = N0(2),N0(1)
               J1 = CONSR(N,M,IK)
               J2 = CONSR(N+1,M,IK)
               J3 = CONSR(N+1,M+1,IK)
               J4 = CONSR(N,M+1,IK)
               WRITE(700,7002) J1,J2,J3,J4
            ENDDO
         ENDDO
      ENDDO

      RETURN
C))))))))))))))))))))) End of subroutine PROPLT ((((((((((((((((((((((((
      END

C --------------------------------------
      FUNCTION CONSR(N,M,KK0)
C ---------------------------------------
C     Connectivity for actual blade (w/o SR)
C ---------------------------------------
      INCLUDE 'PUFCAV.INC'

      NC0=NCOLD+1
      NC1=NC0+1
      CONSR=(KK0-1)*(MRP*NC1)+(M-1)*NC1+(N-N0(2)+1)

      RETURN
      END








