       SUBROUTINE VISFRBL(M,N1,N2,FX,FY,FZ,TMX,TMY,TMZ)
************************************************************************
*      This subroutine applies the viscous correction to the blade     *
*      forces.                                                         *
*                                                                      *
*      Date           Revisions or Comments                            *
*      ------         -----------------------                          *
*      XM Yu 2/20/2012  For viscous force and torque  
************************************************************************

       INCLUDE 'PUFCAV.INC'

       DO 10 N=N1,N2
          L=INDEXB(N,M)
          IF(ISP.NE.1) THEN
             TX=UXTOT(N,M)/SQRT(VTOTS(L))
             TY=UYTOT(N,M)/SQRT(VTOTS(L))
             TZ=UZTOT(N,M)/SQRT(VTOTS(L))
          ELSE
             TX=UL(L,1)
             TY=UL(L,2)
             TZ=UL(L,3)
          END IF
          CFA=cfskin(n,m)*SS(L,1)
          FX=FX+CFA*TX
          FY=FY+CFA*TY
          FZ=FZ+CFA*TZ
          TMX=TMX+CFA*(XCT(L,2)*TZ-XCT(L,3)*TY)
          TMY=TMY+CFA*(XCT(L,3)*TX-XCT(L,1)*TZ)
          TMZ=TMX+CFA*(XCT(L,1)*TY-XCT(L,2)*TX)
 10    CONTINUE
       
       RETURN
       END



