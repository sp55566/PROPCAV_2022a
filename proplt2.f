      SUBROUTINE PROPLT2 
      
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION YB1(NBPZ,MBPZ),ZB1(NBPZ,MBPZ)

C-----------------------------------------------------------------------
C     Plot the propeller for *.plt
C-----------------------------------------------------------------------

 1001 Format(1x,'ZONE T="Propeller", N=',I5,', E=,',I5,
     *     ' F=FEPOINT, ET=QUADRILATERAL')
 1002 FORMAT(1x,F12.9,2x,F12.9,2x,F12.9)               
 1003 FORMAT(1x,I5,2x,I5,2x,I5,2x,I5)
                  
      DELK=TWOPI/NBLADE
      NTOTAL=0

      WRITE(700,1001) (MR+2)*(NC+1)*NBLADE,(MR+1)*NC*NBLADE

      DO 10 K=1, NBLADE

C-------The direction of rotation is wrong (JY011000)-------------------
CJY         T=DELK*FLOAT(K-1)
         T=-DELK*FLOAT(K-1)
         
         DO 20 M=1, MR+1
            DO 30 N=1,NC+1
               YB1(N,M)=YB(N,M)*COS(T)-ZB(N,M)*SIN(T)
               ZB1(N,M)=YB(N,M)*SIN(T)+ZB(N,M)*COS(T)
 30         CONTINUE
 20      CONTINUE

         M=1
         DO 40 N=1,NC+1
            WRITE(700,1002) XB(N,M),YB1(N,M),ZB1(N,M)            
            IF(K.EQ.1)THEN
               NTOTAL=NTOTAL+1
            ENDIF
 40      CONTINUE
         
         DO 50 M=1, MR
            DO 60 N=1, NC+1
               WRITE(700,1002) (XB(N,M)+XB(N,M+1))/2,
     *              (YB1(N,M)+YB1(N,M+1))/2,
     *              (ZB1(N,M)+ZB1(N,M+1))/2               
               IF(K.EQ.1)THEN
                  NTOTAL=NTOTAL+1
               ENDIF
 60         CONTINUE
 50      CONTINUE

         M=MR+1
         DO 70 N=1, NC+1
            WRITE(700,1002) XB(N,M), YB1(N,M), ZB1(N,M)
            IF(K.EQ.1)THEN
               NTOTAL=NTOTAL+1
            ENDIF
 70      CONTINUE

 10   CONTINUE

      DO 80 K=1, NBLADE
         DO 90 M=1, MR+1
            DO 100 N=1, NC
               WRITE(700, 1003) ((M-1)*(NC+1)+N)+(K-1)*NTOTAL,
     *              ((M-1)*(NC+1)+N+1)+(K-1)*NTOTAL,
     *              (M*(NC+1)+N+1)+(K-1)*NTOTAL,
     *              (M*(NC+1)+N)+(K-1)*NTOTAL
 100        CONTINUE
 90      CONTINUE
 80   CONTINUE

      RETURN
C))))))))))))))))))))) End of subroutine PROPLT ((((((((((((((((((((((((
      END











