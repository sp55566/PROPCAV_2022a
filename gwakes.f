!s---YE TIAN ---- 07/24/2012
!     SUBROUTINE GWAKES
      SUBROUTINE GWAKES(NSUB_inp,NWSUB1_inp)
!e---YE TIAN ----
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

C-----------------------------------------------------------------------
C     Parameters of the transition wake
C-----------------------------------------------------------------------
      MN=MRP
!debug ye tian
!     write(*,*) 'calling gwakes'
!     write(*,*) 'nwpanel=', nwpanel
!     write(*,*) 'NSUB_inp=', NSUB_inp
!     do i = 1, mrp
!       write(*,*) 'nsw',i,nsw(i)
!     end do
      NSUB = NSUB_inp
      NWSUB1 = NWSUB1_inp
!     write(*,*) 'NSUB=',NSUB 
!     write(*,*) 'NWSUB1=',NWSUB1 
!debug ye tian

C-----------------------------------------------------------------------
C     Generate the wake subpanels
C-----------------------------------------------------------------------
C.....set NSUB=NWMIN and NWSUB1=1 so that we actually do not sub-devide
C.....any of the wake panels.  Also, this will allow the cavity to 
C.....grow all the way to NWMIN.  (JY061499)

C.....NSUB is the number of macro-panels to subdivide...................
C-----Change NSUB, for viscous boundary layer analysis---by Hong Sun---
!s--YE TIAN move this out to propcav.f  07/24/2012
!     NSUB=5
!     IF(IVISC.EQ.1) THEN
!       WRITE(*,*) 
!       WRITE(*,*) ' PROPCAV> INPUT NSUB'
!       WRITE(*,*)'          (no. of macro-panels in the wake', 
!    *          ' to subdivide):'
!       READ(*,*) NSUB
!     ENDIF
!e--YE TIAN move this out to propcav.f 07/24/2012

C.....NWSUB1 is the number of panels each macro-panel will be devided...
C.....into..............................................................
C
!s--YE TIAN move this out to propcav.f  07/24/2012
!C.....Special adjustment for unsteady wake alignment (IAN=2).  
!C.....HSL 10/12/01
!      IF(IAN.NE.2) THEN
!         NWSUB1=4
!C-----Change NWSUB1, for viscous boundary layer analysis--by Hong Sun--
!        IF(IVISC.EQ.1) THEN  
!          WRITE(*,*) ' PROPCAV> INPUT NWSUB1'  
!          WRITE(*,*) '          (no. of panels each wake macro-panel', 
!     *               ' to be subdivided into):'
!          READ(*,*) NWSUB1  
!        ENDIF
!      ELSE
!         NWSUB1=2
!      END IF
!e--YE TIAN move this out to propcav.f  07/24/2012

C.....NWSUB is the total number of subpanels............................
!     NSUB=5
!     NWSUB1=4
      NWSUB=NSUB*NWSUB1

C-----------------------------------------------------------------------
C     Generate the transition wake geometry
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     I'm changing the way it's defining the coordinates of the sub-
C     panels here because before the pitch of the subpanels don't match
C     with that of the macropanels in the wake.                 JY041198
C-----------------------------------------------------------------------
!     write(*,*) 'NSUB =', NSUB
      
!      IF(NSUB.EQ.1.and.NWSUB1.EQ.1) THEN
!         DO 9999  M=1,MN
!           nw_tmp = nsw(m)
!           if(ian.eq.6) nw_tmp = nwpanel+1
!           do 999  n=1,nw_tmp
!           xwvs3d(n,m)=xw(n,m)
!           ywvs3d(n,m)=yw(n,m)
!           zwvs3d(n,m)=zw(n,m)
! 999       CONTINUE
! 9999    CONTINUE
!      ELSE
       DO 100 M=1,MN           
c........Initial coordinates of the transition wake.....................
         XWS(1,M)=XW(1,M)
         YWS(1,M)=YW(1,M)
         ZWS(1,M)=ZW(1,M)

         DO 110 N=1,NSUB
            DWX=(XW(N+1,M)-XW(N,M))/FLOAT(NWSUB1)
            DWY=(YW(N+1,M)-YW(N,M))/FLOAT(NWSUB1)
            DWZ=(ZW(N+1,M)-ZW(N,M))/FLOAT(NWSUB1)
            DO 120 L=1,NWSUB1
               NIDX=(N-1)*NWSUB1+L+1
               XWS(NIDX,M)=XWS(NIDX-1,M)+DWX
               YWS(NIDX,M)=YWS(NIDX-1,M)+DWY
               ZWS(NIDX,M)=ZWS(NIDX-1,M)+DWZ
 120        CONTINUE
 110     CONTINUE

c XM YU 12/2011------------------------------------------------
c Regenerate the wake geometry including subpanels for the influence
c coefficients calculation of viscous run
         do 130 n=1,nwsub
            xwvs3d(n,m)=xws(n,m)
            ywvs3d(n,m)=yws(n,m)
            zwvs3d(n,m)=zws(n,m)
 130     continue
         
         nw_tmp = nsw(m)
         if(ian.eq.6) nw_tmp = nwpanel+1
C/s S.N.KIM | Unsteay wake alignment
         if(ian.eq.2) nw_tmp = nwpanel+1
C/e S.N.KIM | Aug. 2018.

         do 140 n=nsub+1,nw_tmp
            nidx=n+nwsub-nsub
            xwvs3d(nidx,m)=xw(n,m)
            ywvs3d(nidx,m)=yw(n,m)
            zwvs3d(nidx,m)=zw(n,m)
 140     continue

!! S.N.Kim : Above Process is writing the wake geometry including subpanels into 
!! xwvs3d, ywvs3d, and zwvs3d. 

c test   
c         open(171,file='subpanel')
c         open(176,file='wakepanel')
c         open(181,file='totalwakep')
c         do n=1,nwsub 
c            write(171,*) xws(n,m),yws(n,m),zws(n,m)
c         enddo
c         do n=1,nsw(m)
c            write(176,*) xw(n,m),yw(n,m),zw(n,m)
c         enddo
c         do n=1,nsw(m)+nwsub-nsub
c            write(181,*) xwvs3d(n,m),ywvs3d(n,m),zwvs3d(n,m)
c         enddo
c         write(*,*) nsw(m),nwsub,nsub,nsw(m)+nwsub-nsub
c         stop
c test
c XM YU 12/2012 -------------------------------------------------
 100     CONTINUE
!     ENDIF
      
c      OPEN(8541, FILE="Subpanel.dat", STATUS='unknown')
c      do j=1, MRP
c       do i=1, nsw(j)
c        write(8541,*) XWS(i,j),YWS(i,j),ZWS(i,j),XW(i,j),YW(i,j),ZW(i,j)
c       ! write(8541,*) undefined(i,j),undefined(i,j),undefined(i,j)
c       enddo
c      enddo
c      CLOSE(8541)



      RETURN
      END



