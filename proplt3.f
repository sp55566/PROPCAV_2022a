      SUBROUTINE PROPLT3 
      
      INCLUDE 'PUFCAV.INC'

C-----------------------------------------------------------------------
C     Plot the propeller for *.vec
C-----------------------------------------------------------------------
 5040 FORMAT(1x,'ZONE T="Propeller", N=',I5,', E=,',I5,
     *     ' F=FEPOINT, ET=QUADRILATERAL')
 5060 FORMAT(7(1X,F10.6))
 5080 FORMAT(4(1X,I5))

      ZOFF=1.0
      IADD=(MR+1)*(NH+1)
      IADD1=MR*NH

C-----------------------------------------------------------------------
C     Construct the blade mesh 
C-----------------------------------------------------------------------
      WRITE(65,5040) (MR+1)*(NH+1)*2,(MR)*NC
      
      DO M=1,MR+1         
         DO N=1,NH+1
            WRITE(65,5060) XB(N,M),YB(N,M),ZB(N,M),0.,0.,0.,0.
         END DO
      END DO
      
      DO M=1,MR+1         
         DO N=NHP,NC+1
            WRITE(65,5060) XB(N,M),YB(N,M),ZB(N,M)+ZOFF,0.,0.,0.,0.
         END DO
      END DO

      DO M=1,MR
         DO N=1,NH
            WRITE(65,5080) (M-1)*(NH+1)+N,(M-1)*(NH+1)+N+1,
     *           M*(NH+1)+N+1,M*(NH+1)+N
         END DO
      END DO
      
      DO M=1,MR
         DO N=1,NH
            WRITE(65,5080) (M-1)*(NH+1)+N+IADD,
     *           (M-1)*(NH+1)+N+1+IADD,M*(NH+1)+N+1+IADD,
     *           M*(NH+1)+N+IADD
         END DO
      END DO

      RETURN
      END
