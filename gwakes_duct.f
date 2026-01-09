      SUBROUTINE GWAKES_DUCT
************************************************************************
*     GWAKES: Geometry of the WAKE Subpanels                           *
*     ---  Generate geometries of the wake subpanels                   *
*          wake                                                        *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     CM010898 Added straight wake for hydrofoil case                  *
*     JY041198 Corrected disagreement between wake subpanels and wake  *
*              macropanels.                                            *
*     JY120498 Made many changes in the grid generation of the         *
*              subpanels.  I also merged gwakescs.f into this routine  *
*              so all we need is one subroutine to generate the sub-   *
*              panels for both fully wetted and cavitating case.       *
*                                                                      *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      NDSUB=5
C.....Change NSUB, for viscous boundary layer analysis   by Hong Sun
      IF(IVISC.EQ.1.AND.ICAVT.LE.1) THEN
        WRITE(*,*) ' PROPCAV> INPUT NSUB for Bare Duct Case'
        WRITE(*,*)'          (no. of macro-panels in the wake', 
     *          ' to subdivide):'
        READ(*,*) NDSUB
      ENDIF
        
      IF(IAN.NE.2) THEN
         NWDSUB1=4
C.....Change NWSUB1, for viscous boundary layer analysis  by Hong Sun
        IF(IVISC.EQ.1.AND.ICAVT.LE.1) THEN
          WRITE(*,*) ' PROPCAV> INPUT NWSUB1 for Bare Duct Case'   
          WRITE(*,*) '          (no. of panels each wake macro-panel', 
     *               ' to be subdivided into):'
          READ(*,*) NWDSUB1
        ENDIF 
      ELSE
         NWDSUB1=2
      END IF
         
      NWDSUB=NDSUB*NWDSUB1
      
      MN=MDUCTP
         
      DO M=1,MN
         
c........Initial coordinates of the transition wake.....................
         XWSD(1,M)=XDW(1,M)
         YWSD(1,M)=YDW(1,M)
         ZWSD(1,M)=ZDW(1,M)
         
         DO N=1,NDSUB
            DWX=(XDW(N+1,M)-XDW(N,M))/FLOAT(NWDSUB1)
            DWY=(YDW(N+1,M)-YDW(N,M))/FLOAT(NWDSUB1)
            DWZ=(ZDW(N+1,M)-ZDW(N,M))/FLOAT(NWDSUB1)
            DO L=1,NWDSUB1
               NIDX=(N-1)*NWDSUB1+L+1
               XWSD(NIDX,M)=XWSD(NIDX-1,M)+DWX
               YWSD(NIDX,M)=YWSD(NIDX-1,M)+DWY
               ZWSD(NIDX,M)=ZWSD(NIDX-1,M)+DWZ
            ENDDO
         ENDDO
      ENDDO
 
      RETURN
      END



