       SUBROUTINE VECPLT2
************************************************************************
*                                                                      *
*  Subroutine WETted VECtor plots out the wetted total velocity        *
*  vectors on the blade for the last wetted revolution.  The data      *
*  are arranged for TECPLOT format.                                    *
*                                                                      *
*                           ----------        ----------               *
*                          |          |      |          |              *
*                          |  PROPCAV |> --->|  WETVEC  |              *
*                          |          |      |          |-->           *
*                           -----^----        ----------   |           *
*                                |                         |           *
*                                 --<---------<------------v           *
*                                                                      *
*  Author: Julie Young 100598                                          *
*  Date:      Revision/comments                                        *
*  -----      ---------------                                          *
*  JY092099   Modified routine to also plot vectors in wake (subpanel) *
*             region.                                                  *
*  JY080200   Modified subroutine to plot vectors on the face and back *
*             side of the blades separately.  Also, vectors on wake    *
*             will no longer be plotted.                               *
*                                                                      *
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       OPEN(965,FILE='check.vec',STATUS='UNKNOWN')
 5020  FORMAT(1X,'VARIABLES="X","Y","Z","U","V","W","QC"')
       WRITE(965,5020)

 5040  FORMAT(1x,'ZONE T="Propeller", N=',I5,', E=,',I5,
     *      ' F=FEPOINT, ET=QUADRILATERAL')
 5060  FORMAT(7(1X,F10.6))
 5080  FORMAT(4(1X,I5))
          
       ZOFF=1.0
       IADD=(MR+1)*(NH+1)
       IADD1=MR*NH

C-----------------------------------------------------------------------
C     Construct the blade mesh 
C-----------------------------------------------------------------------
       WRITE(965,5040) (MR+1)*(NH+1)*2,(MR)*NC
             
       DO M=1,MR+1         
          DO N=1,NH+1
             WRITE(965,5060) XB(N,M),YB(N,M),ZB(N,M),0.,0.,0.,0.
          END DO
       END DO
       
       DO M=1,MR+1         
          DO N=NHP,NC+1
             WRITE(965,5060) XB(N,M),YB(N,M),ZB(N,M)+ZOFF,0.,0.,0.,0.
          END DO
       END DO

       DO M=1,MR
          DO N=1,NH
             WRITE(965,5080) (M-1)*(NH+1)+N,(M-1)*(NH+1)+N+1,
     *            M*(NH+1)+N+1,M*(NH+1)+N
          END DO
       END DO
          
       DO M=1,MR
          DO N=1,NH
             WRITE(965,5080) (M-1)*(NH+1)+N+IADD,
     *            (M-1)*(NH+1)+N+1+IADD,M*(NH+1)+N+1+IADD,
     *            M*(NH+1)+N+IADD
          END DO
       END DO
       
C-----------------------------------------------------------------------
C      Plot velocity vectors on the blade
C-----------------------------------------------------------------------
       WRITE(965,5100) TT(IDXREV),MR*NH*2,(MR-1)*(NH-1)*2
 5100  FORMAT(1x,'ZONE T="T=',F4.0,'", N=',I5,', E=,',I5,
     *      ' F=FEPOINT, ET=QUADRILATERAL')

       DO M=1,MR
          DO N1=1,NH
             J1=INDEXB(N1,M)
C             DU=UXTOT(N1,M)
C             DV=UYTOT(N1,M)
C             DW=UZTOT(N1,M)
C             DQC=SQRT(DU**2.+DV**2.+DW**2.)      
             UT=VOX1(J1)*VL(J1,1)+VOY1(J1)*VL(J1,2)+VOZ1(J1)*VL(J1,3)
             DQC=UT+DPDVB(N1,M)
             DU=DQC*VL(J1,1)
             DV=DQC*VL(J1,2)
             DW=DQC*VL(J1,3)
             WRITE(965,5060) XCT(J1,1),XCT(J1,2),XCT(J1,3),DU,DV,DW,DQC
          END DO
       END DO

       DO M=1,MR
          DO N1=NH+1,NC
             J1=INDEXB(N1,M)
C             DU=UXTOT(N1,M)
C             DV=UYTOT(N1,M)
C             DW=UZTOT(N1,M)
C             DQC=SQRT(DU**2.+DV**2.+DW**2.)  
             UT=VOX1(J1)*VL(J1,1)+VOY1(J1)*VL(J1,2)+VOZ1(J1)*VL(J1,3)
             DQC=UT+DPDVB(N1,M)
             DU=DQC*VL(J1,1)
             DV=DQC*VL(J1,2)
             DW=DQC*VL(J1,3)              
             WRITE(965,5060) XCT(J1,1),XCT(J1,2),XCT(J1,3)+ZOFF,
     *            DU,DV,DW,DQC
          END DO
       END DO

       DO M=1,MR-1
          DO N=1,NH-1
             WRITE(965,5080) (M-1)*(NH)+N,(M-1)*(NH)+N+1,
     *            M*(NH)+N+1,M*(NH)+N
          END DO
       END DO
       
       DO M=1,MR-1
          DO N=1,NH-1
             WRITE(965,5080) (M-1)*(NH)+N+IADD1,
     *            (M-1)*(NH)+N+1+IADD1,M*(NH)+N+1+IADD1,
     *            M*(NH)+N+IADD1
          END DO
       END DO

       CLOSE(965)

       RETURN
       END
