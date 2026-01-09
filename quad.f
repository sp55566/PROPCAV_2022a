C    
C [BQUAD.FOR] INTERPOLATION BY 2ND ODER POLYNOMIALS    
C    
      SUBROUTINE BQUAD(X,Y,XX,IM,COA,COB,COC,YY,N)     
      dimension X(*),Y(*)   
C------- X      : ARRAY  
C------- Y      : ARRAY  
C------- XX     : POINT TO BE INTEPOLATED    
C------- N      : NO. OF X OR Y    
C------- YY     : INTERPOLATED VALUE    
      IM=0     
      IF(XX.LE.X(2)) IM=2     
      IF(XX.GT.X(N-2)) IM=N-1 
      IF(XX.GT.X(N-1)) IM=N-1 
      IF(IM.NE.0) GOTO 200    
      DO 100 I=2,N-2     
         IF(XX.GT.X(I).AND.XX.LE.X(I+1)) THEN     
            IM=I  
            GOTO 200     
         ENDIF 
 100      CONTINUE  
 200      IA=IM-1   
      IB=IM    
      IC=IM+1  
      A=X(IA)  
      B=X(IB)  
      C=X(IC)  
      YA=Y(IA) 
      YB=Y(IB) 
      YC=Y(IC) 
      P=(A-B)*(B-C)*(C-A)     
      COA=-(YA*(B-C)+YB*(C-A)+YC*(A-B))/P    
      COB=(YA*(B**2-C**2)+YB*(C**2-A**2)+YC*(A**2-B**2))/P  
      COC=-(YA*(B-C)*B*C+YB*(C-A)*C*A+YC*(A-B)*A*B)/P  
      YY=COA*XX**2+COB*XX+COC 
      RETURN   
      END 
C    
C [DQUAD.FOR] SLOPE USING BQUAD    
C    
      FUNCTION DQUAD(X,Y,XX,N)   
      dimension X(*),Y(*)     
C------- X      : ARRAY  
C------- Y      : ARRAY  
C------- XX     : POINT TO BE EVALUATED 
C------- N      : NO. OF X OR Y    
C------- YY     : SLOPE  
      CALL BQUAD(X,Y,XX,IM,COA,COB,COC,YY,N) 
      DQUAD = 2.0 * COA * XX + COB    
      RETURN   
      END 
C    
C [QUADITG.FOR] INTEGRATION USING 2ND ORDER POLYNOMIAL 
C    
      FUNCTION QUADITG(X,Y,BU,BL,N)     
      dimension X(*),Y(*)   
C------- X      : ARRAY  
C------- Y      : ARRAY  
C------- BU     : UPPER BOUNDARY   
C------- BL     : LOWER BOUNDARY   
C------- N      : NO. OF X OR Y    
      CALL BQUAD(X,Y,BL,IML,COAL,COBL,COCL,YBL,N) 
      CALL BQUAD(X,Y,BU,IMU,COAU,COBU,COCU,YBU,N) 
      JM=IMU-IML    
      IF(JM.EQ.0.OR.JM.EQ.1) THEN  
         QUADITG=COAL/3.0*(BU**3-BL**3)+COBL/2.0*(BU**2-BL**2) 
     #                 +COCL*(BU-BL)    
         RETURN     
      ENDIF    
      QUADITG=COAL/3.0*(X(IML)**3-BL**3)+COBL/2.0*(X(IML)**2-BL**2) 
     #                 +COCL*(X(IML)-BL)     
      DO 100 J=1,JM-1    
         XXL=X(IML+J-1)  
         XXU=X(IML+J)    
         CALL BQUAD(X,Y,XXU,IMX,COAX,COBX,COCX,YYU,N)  
         QUADITG=QUADITG+COAX/3.0*(XXU**3-XXL**3)     
     #               +COBX/2.0*(XXU**2-XXL**2)+COCX*(XXU-XXL)   
 100      CONTINUE  
      QUADITG=QUADITG+COAU/3.0*(BU**3-X(IMU-1)**3)    
     #            +COBU/2.0*(BU**2-X(IMU-1)**2)+COCU*(BU-X(IMU-1))   
      RETURN   
      END 
