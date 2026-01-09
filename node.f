       SUBROUTINE NODE(IP,N,M,XN)
************************************************************************
*      This subroutine determines the node index of element I.  The    *
*      nodes are numbered in the same fashion as in CONPT.             
*                                                                      *
*      Input:                                                          *
*      --------------------------------------------------------------  *
*      IP = 1 : Blade                                                  *
*           2 : HUB                                                    *
*           3 : Wake                                                   *
*                                                                      *
*      Output:                                                         *
*      --------------------------------------------------------------  *
*      XN(I,K) : k-coordinate of the I node for panel (N,M)            *
************************************************************************
       INCLUDE 'PUFCAV.INC'
       DIMENSION XN(4,3)

       IF(IP.EQ.1) THEN
C........nodes on the blade
          XN(1,1)=XB(N,M)
          XN(1,2)=YB(N,M)
          XN(1,3)=ZB(N,M)

          XN(2,1)=XB(N,M+1)
          XN(2,2)=YB(N,M+1)
          XN(2,3)=ZB(N,M+1)

          XN(3,1)=XB(N+1,M+1)
          XN(3,2)=YB(N+1,M+1)
          XN(3,3)=ZB(N+1,M+1)

          XN(4,1)=XB(N+1,M)
          XN(4,2)=YB(N+1,M)
          XN(4,3)=ZB(N+1,M) 
       ELSE IF(IP.EQ.2) THEN
C........nodes on the hub
          XN(1,1)=XH(N,M+1)
          XN(1,2)=YH(N,M+1)
          XN(1,3)=ZH(N,M+1)

          XN(2,1)=XH(N,M)
          XN(2,2)=YH(N,M)
          XN(2,3)=ZH(N,M)

          XN(3,1)=XH(N+1,M)
          XN(3,2)=YH(N+1,M)
          XN(3,3)=ZH(N+1,M)

          XN(4,1)=XH(N+1,M+1)
          XN(4,2)=YH(N+1,M+1)
          XN(4,3)=ZH(N+1,M+1) 
       ELSE IF(IP.EQ.3) THEN
C........nodes on the wake
          XN(1,1)=XW(N,M)
          XN(1,2)=YW(N,M)
          XN(1,3)=ZW(N,M)

          XN(2,1)=XW(N,M+1)
          XN(2,2)=YW(N,M+1)
          XN(2,3)=ZW(N,M+1)

          XN(3,1)=XW(N+1,M+1)
          XN(3,2)=YW(N+1,M+1)
          XN(3,3)=ZW(N+1,M+1)

          XN(4,1)=XW(N+1,M)
          XN(4,2)=YW(N+1,M)
          XN(4,3)=ZW(N+1,M) 
       END IF

       RETURN
       END
