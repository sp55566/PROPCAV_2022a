

      FUNCTION DPDXQ2(X1,X2,P0,P1,P2,IFR)
************************************************************************
*     DPDXQ (2nd order):                                               *
*      --- Forward difference formula (IFR=1)                          *
*                                                                      *
*                   0      1      2                                    *
*         DP/DX if  +------+------+  ---> X                            *
*                      x1    x2                                        *
*                                                                      *
*      --- Backward difference formula (IFR=2)                         *
*                                                                      *
*                           2      1      0                            *
*         DP/DX if  X <---  +------+------+                            *
*                              x2    x1                                *
************************************************************************
      
      DNOM=X1**2.*X2+X1*X2**2.
      DUM1=X2*(2.*X1+X2)/DNOM
      DUM2=(X1+X2)**2./DNOM
      DUM3=X1**2./DNOM
      IF(IFR.EQ.1) THEN
         DPDXQ2=-DUM1*P0+DUM2*P1-DUM3*P2
      ELSE IF(IFR.EQ.2) THEN
         DPDXQ2=DUM1*P0-DUM2*P1+DUM3*P2
      END IF
      RETURN
C)))))))))))))))))))))) End of function DPDXQ ((((((((((((((((((((((((((
      END
