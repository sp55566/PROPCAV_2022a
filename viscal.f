c
c    *****************   SUBROUTINE VISCAL     *******************
c
c    Couple with the BL equations through wall transpiration model     
c
c
c XM YU 2/2012-----change name, cdvis to cfvis
      SUBROUTINE VISCAL(NPAN1,NPANW,SPAN,XC,YC,XTRANS,XTRANP,RE,
     & AOA,UI,APANIJ,NDIM,UEDGE,DSTARS,DSTARP,THETAS,THETAP,ICAV,
     & ILEP,ITEP,IICLE,IICTE,XNV,YNV,RVCRIT,MAXIT,EPSL,cfvis,
     & LVISCON, BETA_BL_INP)
c XM YU 2/2012

C------------------------------------------------------------------
C     Converges viscous operating point
C------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)   
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'

      REAL  SPAN(ndim),XC(ndim),YC(ndim),UI(ndim),XNV(ndim),
     &      YNV(ndim),APANIJ(ndim,ndim)
      REAL  RE,RVCRIT,AOA,XTRANS,XTRANP,EPSL

C      LOGICAL LFORCE, LSYMB, LREF       
C      INTEGER NBLMAX(2)
      INTEGER  LVISCON
      DIMENSION  UEDGE(ndim),DSTARS(ndim),DSTARP(ndim)
     &     ,     THETAS(ndim),THETAP(ndim)
     &     ,     X0(ndim),Y0(ndim),XN(ndim),YN(ndim),DS(ndim)
      DIMENSION  HH(280)
cc XM YU 2/17/2012
      dimension cfvis(ndim)
cc XM YU
!s--- YE TIAN ---- debug
      LOGICAL FL_ERR
      integer:: viscal_count=0


      BETA_BL = BETA_BL_INP
      viscal_count = viscal_count + 1
!     write(*,*) 'viscal.f:36'
!     write(*,*) 'vis_count=', viscal_count

!     write(1236,*) 'viscal_count=', viscal_count
!     write(1236,*) NPAN1, NPANW
!     do i = 1, NPAN1+NPANW
!       write(1236,*) i, SPAN(i), UI(i)
!     end do

!e--- YE TIAN ----


c
c      write(*,*) 'start of viscal'
      Nvis = NPAN1
      NWvis = NPANW

        IWRT=0
        ICV=ICAV
C      ICLE=IICLE
C      ICTE=IICTE
1010      LALFA =.true.

      DTOR = 3.141592654D0/180.D0

      ARAD = DBLE(AOA)*DTOR                  
        EPS1 = DBLE(EPSL)
        ACRIT = DBLE(RVCRIT)
  
      do i=1,nvis+nwvis
            s(i)=DBLE(span(i))
            x(i)=DBLE(xc(i))
            y(i)=DBLE(yc(i))
            xn(i)=DBLE(xnv(i))
            yn(i)=DBLE(ynv(i))
            qinv(i)=DBLE(ui(i))
c        write(71,110) i, s(i),x(i),y(i), qinv(i)     
      enddo
c 110            format(I4,2x,4f10.5)          
 110            format(5x,f10.5)

c       Added 'ds' calculation by Hong Sun
        DO J = 1, Nvis+NWvis-1
          DS(J) = S(J+1)-S(J)
        ENDDO
        DS(Nvis+NWvis)=DS(Nvis+NWvis-1)

c  obtain DX/DS,DY/DS for SINVRT routine
c        write(*,*) 'before segspl'        
      call segspl(x,xp,s,nvis)
      call segspl(y,yp,s,nvis)

c  leading edge arc length SLE
      call lefind(sle,x,xp,y,yp,s,nvis) 
c       do i=1, nvis                                     
c         write(71,*) i, xp(i), yp(i)                  
c        enddo

        do i=1,nvis+nwvis
           do j=1,nvis+nwvis
                dij(i,j)=DBLE(apanij(i,j))
           enddo
c          WRITE(71,111) (dij(I,J),J=1,Nvis+NWvis)       
        enddo
c 111      FORMAT(1X,7F10.4)

      xtr1(1)=DBLE(xtrans)                  ! @@
      xtr1(2)=DBLE(xtranp)                  ! @@
c  force transition at the beginning of the wake
      if(xtr1(1).gt.1.d0)xtr1(1)=1.D0
      if(xtr1(2).gt.1.d0)xtr1(2)=1.D0

      reinf=dble(re)
c        write(71,*) 'aoa =', ARAD, 'sle =', sle  
c        write(71,*) 'EPS1 =', EPS1
c      write(71,*)'re=',reinf,' nt+1 =',nvis,'ntw+1=',nwvis+nvis
c      write(71,*)'xtrans =',xtr1(1),'xtranp=',xtr1(2)

C---- convergence tolerance
c      EPS1 = 1.0E-4
c
C---- calculate wake trajectory from current inviscid solution if necessary
cccc      IF(.NOT.LWAKE) CALL XYWAKE
C
C---- set velocities on wake from airfoil vorticity for alpha=0, 90
cccc      CALL QWCALC
C
C---- set velocities on airfoil and wake for initial alpha
cccc      CALL QISET
C
cccc      IF(.NOT.LIPAN) THEN
C
cccc      IF(LBLINI) CALL GAMQV
c        write(71,*) 'before gamqi (Qinv)'
      call GAMQI
c      do i=1,nvis+nwvis
c          write(71,110) gam(i)
c      enddo
C
C----- locate stagnation point arc length position and panel index
c       write(71,*) 'before stfind'
       CALL STFIND
c      write(71,*)'ist = ',ist,'sst = ',sst
c      write(71,*)'sst_go', sst_go,'sst_gp', sst_gp
C
C----- set  BL position -> panel position pointers
       CALL IBLPAN
c      write(71,*)'iblte(1)=',iblte(1),'iblte(2)=',iblte(2)

C----- calculate surface arc length array for current stagnation point location
C      Sharp t.e. or not
      sharp = .true.
c       write(*,*) 'before xicalc'
       CALL XICALC
C
C----- set  BL position -> system line  pointers
c       write(*,*) 'before iblsys'
       CALL IBLSYS
C
cccc      ENDIF
C
C---- set inviscid BL edge velocity UINV from QINV
c       write(*,*) 'before uicalc'
       CALL UICALC
C
cccc      IF(.NOT.LBLINI) THEN
C
C----- set initial Ue from inviscid Ue

       DO 5 IBL=1, NBL(1)
         UEDG(IBL,1) = UINV(IBL,1)
!        WRITE(71,*) UEDG(IBL,1)
    5  CONTINUE

!        write(*,*) 'NBL(1)=', NBL(1)
!        write(*,*) 'stop:viscal.f:158'
!        stop


       DO 6 IBL=1, NBL(2)
         UEDG(IBL,2) = UINV(IBL,2)
c         WRITE(71,*) UEDG(IBL,2)
    6  CONTINUE
      LVCONV = .FALSE.
      lblini =.false.

C   FIND CAVITY LE,TE INDEX FOR BL SOLUTION BY HONG SUN
C   ILEP,ITEP: Use Foil Index Numbering (anti-clockwise)
       DO IB=1,NBL(1)
          if(ILEP.eq.IPAN(IB,1))then
            ICLE=IB
          IICLE=ICLE
          elseif(ITEP.eq.IPAN(IB,1))then
            ICTE=IB
          IICTE=ICTE
          endif
       ENDDO 
!s-- YE TIAN 09/02/2013 ----
!--- The ICTE calculation is messed up ......
!temporarily set as 0
       ICTE = 0
       IICLE = 0
!e-- YE TIAN 09/02/2013 ----
C
C---- Newton iteration for entire BL solution
C      WRITE(*,*) ' '
C      WRITE(*,*) 'Solving BL system ...'
C      WRITE(*,*) ' maxit =',    maxit

      LVISCON=0
      DO 1000 ITER=1, maxit
C 
C------- fill Newton system for BL variables
c         write(*,*) 'filling Newton system for BL variables'
         CALL SETBL(XN,YN,DS,FL_ERR)
         if (FL_ERR) then
           LVISCON = 1
           return
         end if
C
C------- solve Newton system with custom solver
c         write(*,*) 'solving Newton system...'
         CALL BLSOLV
C
C------- update BL variables
c         write(*,*) 'updating BL variables...'
         CALL UPDATE
C
cccc        IF(.NOT. LALFA) THEN
C------- set new inviscid speeds QINV and UINV for new alpha
cccc         CALL QISET
cccc         CALL UICALC
cccc        ENDIF
C
C------ calculate edge velocities QVIS(.) from UEDG(..)
c         write(*,*) 'calculating edge velocities...'
         CALL QVFUE
C
C------ set GAM distribution from QVIS
         CALL GAMQV
C
C------ relocate stagnation point
         CALL STMOVE
C
C------ set CL for CL-dependent Reynolds number for prescribed alfa
cccc        CALL CLCALC
C
C------ display changes and test for convergence
        IF(RLX.LT.1.D0) 
     &   WRITE(*,2000) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL,RLX
        IF(RLX.EQ.1.D0) 
     &   WRITE(*,2010) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL
c------- wbrewer -- 4/10/95 -- comment out convergence history of forces
c         WRITE(*,2020) ARAD/DTOR, CL, CM, cdvis
        write(*,*)
C
        IF(RMSBL .LT. EPS1) THEN
         LVCONV = .TRUE.
         AVISC = ARAD
         GO TO 90
        ENDIF
C
 1000 CONTINUE
C 
      LVISCON=1
!     WRITE(*,'(/A)') 'VISCAL:  Convergence failed'
cccc      CALL CPCALC
cccc      IF(LFLAP) CALL MHINGE
cccc      RETURN
C
   90 CONTINUE
cccc      CALL CPCALC
cccc      IF(LFLAP) CALL MHINGE
c
c      capture duplicate point at TE in x coordinate
c  store solution variables to be passed back to pan2d

c---------------------     suction side first      -------------------

      is=1
        uecle=0.0d0
        icle=0

       IWRT = IWRT + 1
      do 50 ibl=2,iblte(is)
            i=ipan(ibl,is)
            uedge(i)=uedg(ibl,is)
            thetas(i)=thet(ibl,is)
            dstars(i)=dstr(ibl,is)
            hh(i)=dstars(i)/thetas(i)

cs           TEMPU = ABS( UEDGE(I)-UINV(IBL,IS) )  
cs           IF(TEMPU.GT.BLepsl) THEN
cs             BLepsl = TEMPU
cs           ENDIF 

       IF (IWRT.EQ.1) THEN
      write(11,2001)x(i),dstars(i),thetas(i),hh(i),
     &                uedge(i),uinv(ibl,is)
       ENDIF
C        write(97,2001)x(i),dstars(i),thetas(i),hh(i),
C     &                uedge(i),uinv(ibl,is) 
c       write(71,*) i, uedge(i), dstars(i), uedge(i)*dstars(i),
c     &               uinv(ibl,is), vti(ibl,is), s(i)
                if (uecle.lt.uedge(i)) then
                  uecle=uedge(i)
                  icle=ibl
                endif
 50      continue

C        write(71,*) icle, uecle

c----------------      pressure side next      -----------------------
      is=2
C        write(71,*) is, iblte(is), nbl(is)
      do 54 ibl=2,iblte(is)
            i=ipan(ibl,is)
            uedge(i)=uedg(ibl,is)
            thetap(i)=thet(ibl,is)
            dstarp(i)=dstr(ibl,is)
            hh(i)=dstarp(i)/thetap(i)

cs           TEMPU = ABS( UEDGE(I)-UINV(IBL,IS) )  
cs           IF(TEMPU.GT.BLepsl) THEN
cs             BLepsl = TEMPU
cs           ENDIF 
       IF (IWRT.EQ.1) THEN
      write(10,2001)x(i),dstarp(i),thetap(i),hh(i),
     &                uedge(i),uinv(ibl,is)
       ENDIF
c       write(71,*) i,uedge(i),dstarp(i),uedge(i)*dstarp(i),
c     &                uinv(ibl,is), vti(ibl,is), s(i)
 54      continue
c

c---------------------      wake output      -------------------------
      do 56 ibl=iblte(is)+1,nbl(is)
            i=ipan(ibl,is)
            uedge(i)=uedg(ibl,is)
            thetap(i)=thet(ibl,is)
            dstarp(i)=dstr(ibl,is)
            hh(i)=dstarp(i)/thetap(i)

cs           TEMPU = ABS( UEDGE(I)-UINV(IBL,IS) )  
cs           IF(TEMPU.GT.BLepsl) THEN
cs             BLepsl = TEMPU
cs           ENDIF 
       IF (IWRT.EQ.1) THEN
      write(10,2001)x(i),dstarp(i),thetap(i),hh(i),
     &                uedge(i),uinv(ibl,is)
      write(11,2001)x(i),dstarp(i),thetap(i),hh(i),
     &                uedge(i),uinv(ibl,is)
       ENDIF
c      write(*,*) i
 56      continue
 2001      format(6f10.5)

C XM YU---Cfskin calculation 2/17/2012
         que=0.5*qinf**2
         do is=1,2
            do ibl=2,iblte(is)
               i=ipan(ibl,is)
               cfvis(i)=tau(ibl,is)/que
            enddo
         enddo
c XM YU 

C  Find the displacemented foil geometry
      NT1 = NVIS 
      NSTAG = IST
C      write(*,*) ' NSTAG = ', NSTAG

      DO I = 1, Nvis+NWvis
        IF(I.LE.Nvis) THEN
          X0(I) = X(Nvis+1-I)
          Y0(I) = Y(Nvis+1-I)
        ELSE
          X0(I) = X(I)
          Y0(I) = Y(I)          
        ENDIF
       ENDDO

      RDST = DSTARS(1) / (DSTARP(NT1) + DSTARS(1))
       DO 58 I = NT1+NWvis-2, NT1+2, -1
         DFOIX = X0(I)
         DFOIY = Y0(I)
     *         + DSTARP(I)*(YN(I+1)+YN(I))*0.5D0*RDST
c     *          - (DSTARP(NT1+1)-DSTARS(1)+DSTARP(NT1))*0.5D0
C         WRITE(14,*) DFOIX,DFOIY,X0(I),Y0(I)
 58   CONTINUE
       DFOIX = X0(1) + DSTARS(1)*(XN(1)+0.0D0)*0.5D0
       DFOIY = Y0(1) + DSTARS(1)*(YN(1)+1.0D0)*0.5D0
C       WRITE(14,*) DFOIX,DFOIY,X0(1),Y0(1)
       DO 59 I = 2,NT1+NWvis-2
        IF (I.LE.NSTAG) THEN
          DFOIX = X0(I) + DSTARS(I)*(XN(I)+XN(I-1))*0.5D0
          DFOIY = Y0(I) + DSTARS(I)*(YN(I)+YN(I-1))*0.5D0
        ELSE IF (I.LE.NT1) THEN
          DFOIX = X0(I) + DSTARP(I)*(XN(I)+XN(I-1))*0.5D0
          DFOIY = Y0(I) + DSTARP(I)*(YN(I)+YN(I-1))*0.5D0
        ELSE
          DFOIX = X0(I)
          DFOIY = Y0(I)
     *          - DSTARP(I)*(YN(I+1)+YN(I))*0.5D0*(1.0D0-RDST)
c     *          - (DSTARP(NT1+1)-DSTARS(1)+DSTARP(NT1))*0.5D0
        ENDIF
C        WRITE(14,*) DFOIX,DFOIY,X0(I),Y0(I)
 59   CONTINUE
 
C
c  output inviscid and viscid cp plot, and both upper and lower surface
c  shape parameter (hh) plots

      do 60 i=1,nvis
c           IF(I.EQ.1) THEN           
c             UEDGE(I) = UEDGE(I+1) + (S(I)-S(I+1))*
c     &          (UEDGE(I+1)-UEDGE(I+2))/(S(I+1)-S(I+2))
c           ELSE IF(I.EQ.NVIS) THEN
c             UEDGE(I) = -UEDGE(1)
c           ENDIF
            cpii=DBLE(ui(i))**2-1.d0
            cpvv=uedge(i)**2-1.d0
c       IF (IWRT.EQ.1) THEN
c            write(*,*) before write
            write(12,2002)x(i),cpvv,cpii
c            write(*,*)x(i),cpvv,cpii,hh(i)
c            write(*,*) after write
c            stop
          
c       ENDIF
c            write(16,2102)x(i),cpii
c            write(18,2102)x(i),cpvv
 60      continue
2002      format(4f10.5)
c  2102      format(2f10.5)

C
c   compute drag and output

c XM YU comment 02/2012
c      cd=2*thetap(nvis+nwvis)
c        cdvis=cd
c XM YU

c      write(13,2004)cd
c        open(65,file='cd.dat',status='unknown')
c        write(65,2125)cd
c        close(65)
c 2125   format(F7.3)
c 2004      format(' DRAG COEFFICIENT =',F8.5)


      RETURN
C.....................................................................
 2000   FORMAT
     &   (1X,I3,'   rms: ',E11.4,'   max: ',E11.4,3X,A1,' at ',I4,I3,
     &    '   RLX:',F6.3)
 2010   FORMAT
     &   (1X,I3,'   rms: ',E11.4,'   max: ',E11.4,3X,A1,' at ',I4,I3)
 2020   FORMAT
     &   (1X,3X,'   a =' , F7.3,'       CL =' , F8.4  /
     &    1X,3X,'  Cm = ', F7.4, '      CD = ', F8.5 )
      END 
c******************* End of Routine VISCAL *********************
c
c




c
C    *************   SUBROUTINE  SEGSPL   ************
      SUBROUTINE SEGSPL(X,XS,S,N)
C-----------------------------------------------
C     Splines X(S) array just like SPLINE,      |
C     but allows derivative discontinuities     |
C     at segment joints.  Segment joints are    |
C     defined by identical successive S values. |
C-----------------------------------------------
      implicit real*8 (a-h,o-z)
      DIMENSION X(N), XS(N), S(N)
C
      IF(S(1).EQ.S(2)  ) STOP 'SEGSPL:  First input point duplicated'
      IF(S(N).EQ.S(N-1)) STOP 'SEGSPL:  Last  input point duplicated'
C
      ISEG0 = 1
      DO 10 ISEG=2, N-2
        IF(S(ISEG).EQ.S(ISEG+1)) THEN
         NSEG = ISEG - ISEG0 + 1
         CALL SPLINEV(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG)
         ISEG0 = ISEG+1
        ENDIF
   10 CONTINUE
C
      NSEG = N - ISEG0 + 1
      CALL SPLINEV(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG)
C
      RETURN
      END ! SEGSPL
c
c
c
C    *************   SUBROUTINE LEFIND    ************
      SUBROUTINE LEFIND(SLE,X,XP,Y,YP,S,N)
      implicit real*8 (a-h,o-z)
      REAL*8 X(N),Y(N),S(N),XP(N),YP(N)
C
C**** Locates leading edge arc length value SLE
C
C---- set trailing edge point coordinates
      XTE = 0.5D0*(X(1) + X(N))
      YTE = 0.5D0*(Y(1) + Y(N))
C
C---- get first guess for SLE
      DO 10 I=3, N-2
        DXTE = X(I) - XTE
        DYTE = Y(I) - YTE
        DX = X(I+1) - X(I)
        DY = Y(I+1) - Y(I)
        DOTP = DXTE*DX + DYTE*DY
        IF(DOTP .LT. 0.D0) GO TO 11
   10 CONTINUE
C
   11 SLE = S(I)
C
C---- check for sharp LE case
      IF(S(I) .EQ. S(I-1)) RETURN
C
C---- Newton iteration to get exact SLE value
      DO 20 ITER=1, 50
        XLE  = SEVAL(SLE,X,XP,S,N)
        YLE  = SEVAL(SLE,Y,YP,S,N)
        DXDS = DEVAL(SLE,X,XP,S,N)
        DYDS = DEVAL(SLE,Y,YP,S,N)
        DXDD = D2VAL(SLE,X,XP,S,N)
        DYDD = D2VAL(SLE,Y,YP,S,N)
C
        XCHORD = XLE - XTE
        YCHORD = YLE - YTE
        CHDLEN = DSQRT(XCHORD**2+YCHORD**2)
C
C------ drive dot product between chord line and LE tangent to zero
        RES  = XCHORD*DXDS + YCHORD*DYDS
        RESS = DXDS  *DXDS + DYDS  *DYDS
     &       + XCHORD*DXDD + YCHORD*DYDD
C
C------ Newton delta for SLE 
        DSLE = -RES/RESS
C
        DSLE = DMAX1( DSLE , -0.02D0*DABS(XCHORD+YCHORD) )
        DSLE = DMIN1( DSLE ,  0.02D0*DABS(XCHORD+YCHORD) )
        SLE = SLE + DSLE
        IF(DABS(DSLE) .LT. 1.0D-5) RETURN
   20 CONTINUE
      WRITE(*,*) 'LEFIND:  LE point not found.  Continuing...'
      SLE = S(I)
      RETURN
C
      END !LEFIND

c
c
C *****************    FUNCTION D2VAL     *************
      FUNCTION D2VAL(SS,X,XS,S,N)
      implicit real*8 (a-h,o-z)
      DIMENSION X(N), XS(N), S(N)
C--------------------------------------------------
C     Calculates d2X/dS2(SS)                       |
C     XS array must have been calculated by SPLINE |
C--------------------------------------------------
      ILOW = 1
      I = N
C
   10 IF((I-ILOW) .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      D2VAL = (6.D0*T-4.D0)*CX1 + (6.D0*T-2.D0)*CX2
      D2VAL = D2VAL/DS**2
      RETURN
      END ! D2VAL
c
c
c
C   ************** SUBROUTINE QISET   ***********************
      SUBROUTINE QISET
C-------------------------------------------------------
C     Sets inviscid panel tangential velocity for
C     current alpha.
C-------------------------------------------------------
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
      COSA = DCOS(ARAD)
      SINA = DSIN(ARAD)
C
      DO 5 I=1, Nvis+NWvis
        QINV  (I) =  COSA*QINVU(I,1) + SINA*QINVU(I,2)
        QINV_A(I) = -SINA*QINVU(I,1) + COSA*QINVU(I,2)
    5 CONTINUE
C
      RETURN
      END ! QISET
C
C
C   ************** SUBROUTINE GAMQI    ***********************
      SUBROUTINE GAMQI
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
      DO 10 I=1, Nvis
        GAM(I) = Qinv(I)
   10 CONTINUE
C
      RETURN
      END ! GAMQI
C
C
C   ************** SUBROUTINE GAMQV    ***********************
      SUBROUTINE GAMQV
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
      DO 10 I=1, Nvis
        GAM(I) = QVIS(I)
   10 CONTINUE
C
      RETURN
      END  ! GAMQV
C
C
C   ************** SUBROUTINE  STFIND  ********************

      SUBROUTINE STFIND
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
      DO 10 I=1, Nvis-1
        IF(GAM(I).GE.0.D0 .AND. GAM(I+1).LT.0.D0) GO TO 11
   10 CONTINUE
C
      WRITE(*,'(/A)')'STFIND: Stagnation point not found. Continuing...'
      I = Nvis/2
C
   11 CONTINUE
C
      IST = I
      DGAM = GAM(I+1) - GAM(I)
      DS = S(I+1) - S(I)
C
      SST = S(I) - GAM(I)*DS/DGAM
C
      SST_GO = -GAM(I)*DS/DGAM**2 - DS/DGAM
      SST_GP =  GAM(I)*DS/DGAM**2
C
      RETURN
      END  !STFIND
c
c
c


C   ************** SUBROUTINE  IBLPAN  ********************

      SUBROUTINE IBLPAN
C-------------------------------------------------------------
C     Sets  BL location -> panel location  pointer array IPAN
C-------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
C---- suction surface first
      IS = 1
C
      IBL = 1
      DO 10 I=IST, 1, -1
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = 1.D0
   10 CONTINUE
C
      IBLTE(IS) = IBL
      NBL(IS) = IBL
c      write(71,*)'nbl(1)=',nbl(is),'iblte(1)=',iblte(is)
C
C---- pressure surface next
      IS = 2
C
      IBL = 1
      DO 20 I=IST+1, Nvis
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = -1.D0
   20 CONTINUE
C
C---- wake
      IBLTE(IS) = IBL
C
      DO 25 I=Nvis+1, Nvis+NWvis
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = -1.D0
   25 CONTINUE
C
      NBL(IS) = IBL
c      write(71,*)'nbl(2)=',nbl(is),'iblte(2)=',iblte(2)
C
      IF(NBL(1).GT.IVX .OR. NBL(2).GT.IVX) 
     & WRITE(*,*) ' ***  BL array overflow.  Increase IVX.  *** '
C
      LIPAN = .TRUE.
      RETURN
      END !IBLPAN
c
c
C   ************** SUBROUTINE  XICALC  ********************
      SUBROUTINE XICALC
C-------------------------------------------------------------
C     Sets BL arc length array on each airfoil side and wake
C-------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
      IS = 1
C
      XSSI(1,IS) = 0.D0
C
c      write(71,*) 'IS = 1 '
      DO 10 IBL=2, IBLTE(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = SST - S(I)
c      write(71,*)xssi(ibl,is)
   10 CONTINUE
C
C
      IS = 2
C
      XSSI(1,IS) = 0.D0
C
c      write(71,*) 'IS = 2 '
      DO 20 IBL=2, IBLTE(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = S(I) - SST
c      write(71,*)xssi(ibl,is)
   20 CONTINUE
C
      IBL = IBLTE(IS) + 1
      XSSI(IBL,IS) = XSSI(IBL-1,IS)
C
      DO 25 IBL=IBLTE(IS)+2, NBL(IS)
        I = IPAN(IBL,IS)
c        XSSI(IBL,IS) = XSSI(IBL-1,IS)
c     &               + DSQRT((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
       xssi(ibl,is)=s(i)-sst
c      write(71,*)xssi(ibl,is)
   25 CONTINUE
C
C---- trailing edge flap length to TE gap ratio
      TELRAT = 2.50D0
C
C---- set up parameters for TE flap cubics
C     Hong corrected DWDXTE according to XFOIL 6.96
cc      DWDXTE = YP(1)/XP(1) + YP(Nvis)/XP(Nvis)    !!! BUG  2/2/95
C
      CROSP = (XP(1)*YP(Nvis) - YP(1)*XP(Nvis))
     &      / DSQRT(  (XP(1)**2 + YP(1)**2)
     &               *(XP(Nvis)**2 + YP(Nvis)**2) )
      DWDXTE = CROSP / DSQRT(1.D0 - CROSP**2)
C
C---- limit cubic to avoid absurd TE gap widths
      DWDXTE = MAX(DWDXTE,-3.D0/TELRAT)
      DWDXTE = MIN(DWDXTE, 3.D0/TELRAT)

      AA =  3.D0 + TELRAT*DWDXTE
      BB = -2.D0 - TELRAT*DWDXTE
C
      IF(SHARP) THEN
       DO 30 IW=1, NWvis
         WGAP(IW) = 0.D0
   30  CONTINUE
      ELSE
C----- set TE flap (wake gap) array
       IS = 2
       DO 35 IW=1, NWvis
         IBL = IBLTE(IS) + IW
         ZN = 1.D0 - (XSSI(IBL,IS)-XSSI(IBLTE(IS),IS)) / (TELRAT*ANTE)
         WGAP(IW) = 0.D0
         IF(ZN.GE.0.D0) WGAP(IW) = ANTE * (AA + BB*ZN)*ZN**2
   35  CONTINUE
      ENDIF
C
      RETURN
      END !XICALC
C
C
C
C   ************** SUBROUTINE  UICALC  ********************
      SUBROUTINE UICALC
C--------------------------------------------------------------
C     Sets inviscid Ue from panel inviscid tangential velocity
C--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
      DO 10 IS=1, 2
c      write(71,*)'ZONE T=" IS =', IS, ' " '
        UINV  (1,IS) = 0.D0
        UINV_A(1,IS) = 0.D0
c      write(71,*) 1, uinv(1,is), UINV_A(1,IS)
        DO 110 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          UINV  (IBL,IS) = VTI(IBL,IS)*QINV  (I)
          UINV_A(IBL,IS) = VTI(IBL,IS)*QINV_A(I)
c      write(71,*) ibl,uinv(ibl,is),uinv_a(ibl,is)
  110   CONTINUE
   10 CONTINUE
C
      RETURN
      END !UICALC
C
C
C
C   ************** SUBROUTINE QVFUE   ********************  
      SUBROUTINE QVFUE
C--------------------------------------------------------------
C     Sets panel viscous tangential velocity from viscous Ue
C--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
      DO 1 IS=1, 2
c      write(91,*)'iside = ',is
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          QVIS(I) = VTI(IBL,IS)*UEDG(IBL,IS)
c      WRITE(64,*) UEDG(IBL,IS)
c        write(91,*)qvis(i),UEDG(IBL,IS)  
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END   ! QVFUE
C
C
C
C   **************    SUBROUTINE STMOVE    **************************
      SUBROUTINE STMOVE
      implicit real*8 (a-h,o-z)
      INCLUDE 'XFOIL.INC'
C
C---- locate new stagnation point arc length SST from GAM distribution
      ISTOLD = IST
      SSTOLD = SST
      CALL STFIND
C
C---- change in stagnation point arc length
CCC      DSST = SST - SSTOLD
C
      IF(ISTOLD.EQ.IST) THEN
C
C----- recalculate new arc length array
       CALL XICALC
C
      ELSE
C
CCC       WRITE(*,*) 'STMOVE: Resetting stagnation point'
C
C----- set new BL position -> panel position  pointers
       CALL IBLPAN
C
C----- set new inviscid BL edge velocity UINV from QINV
       CALL UICALC
C
C----- recalculate new arc length array
       CALL XICALC
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS
C
       IF(IST.GT.ISTOLD) THEN
C------ increase in number of points on suction side (IS=1)
        IDIF = IST-ISTOLD
C
        ITRAN(1) = ITRAN(1) + IDIF
        ITRAN(2) = ITRAN(2) - IDIF
C
C------ move suction side BL variables downstream
        DO 110 IBL=NBL(1), IDIF+2, -1
          CTAU(IBL,1) = CTAU(IBL-IDIF,1)
          THET(IBL,1) = THET(IBL-IDIF,1)
          DSTR(IBL,1) = DSTR(IBL-IDIF,1)
          UEDG(IBL,1) = UEDG(IBL-IDIF,1)
  110   CONTINUE            
C
C------ set BL variables between old and new stagnation point
        DUDX = UEDG(IDIF+2,1)/XSSI(IDIF+2,1)
        DO 115 IBL=IDIF+1, 2, -1
          CTAU(IBL,1) = CTAU(IDIF+2,1)
          THET(IBL,1) = THET(IDIF+2,1)
          DSTR(IBL,1) = DSTR(IDIF+2,1)
          UEDG(IBL,1) = DUDX * XSSI(IBL,1)
  115   CONTINUE
C
C------ move pressure side BL variables upstream
        DO 120 IBL=2, NBL(2)
          CTAU(IBL,2) = CTAU(IBL+IDIF,2)
          THET(IBL,2) = THET(IBL+IDIF,2)
          DSTR(IBL,2) = DSTR(IBL+IDIF,2)
          UEDG(IBL,2) = UEDG(IBL+IDIF,2)
  120   CONTINUE            
C
       ELSE
C------ increase in number of points on pressure side (IS=2)
        IDIF = ISTOLD-IST
C
        ITRAN(1) = ITRAN(1) - IDIF
        ITRAN(2) = ITRAN(2) + IDIF
C
C------ move pressure side BL variables downstream
        DO 210 IBL=NBL(2), IDIF+2, -1
          CTAU(IBL,2) = CTAU(IBL-IDIF,2)
          THET(IBL,2) = THET(IBL-IDIF,2)
          DSTR(IBL,2) = DSTR(IBL-IDIF,2)
          UEDG(IBL,2) = UEDG(IBL-IDIF,2)
  210   CONTINUE            
C
C------ set BL variables between old and new stagnation point
        DUDX = UEDG(IDIF+2,2)/XSSI(IDIF+2,2)
        DO 215 IBL=IDIF+1, 2, -1
          CTAU(IBL,2) = CTAU(IDIF+2,2)
          THET(IBL,2) = THET(IDIF+2,2)
          DSTR(IBL,2) = DSTR(IDIF+2,2)
          UEDG(IBL,2) = DUDX * XSSI(IBL,2)
  215   CONTINUE
C
C------ move suction side BL variables upstream
        DO 220 IBL=2, NBL(1)
          CTAU(IBL,1) = CTAU(IBL+IDIF,1)
          THET(IBL,1) = THET(IBL+IDIF,1)
          DSTR(IBL,1) = DSTR(IBL+IDIF,1)
          UEDG(IBL,1) = UEDG(IBL+IDIF,1)
  220   CONTINUE            
       ENDIF
C
      ENDIF
C
C---- set new mass array since Ue has been tweaked
      DO 50 IS=1, 2
        DO 510 IBL=2, NBL(IS)
c        WRITE(64,*) DSTR(IBL,IS)
c        WRITE(67,*) UEDG(IBL,IS)
          MASS(IBL,IS) = DSTR(IBL,IS)*UEDG(IBL,IS)
  510   CONTINUE
   50 CONTINUE
C
      RETURN
      END  !  STMOVE 
c
c 
c
C    *************   SUBROUTINE  SPLINEV   ************ 

      SUBROUTINE SPLINEV(X,XS,S,N)
      implicit real*8 (a-h,o-z)
      DIMENSION X(N),XS(N),S(N)
      PARAMETER (NMAX=300)
      DIMENSION A(NMAX),B(NMAX),C(NMAX)
C-------------------------------------------------------
C     Calculates spline coefficients for X(S).          |
C     Zero 2nd derivative end conditions are used.      |
C     To evaluate the spline at some value of S,        |
C     use SEVAL and/or DEVAL.                           |
C                                                       |
C     S        independent variable array (input)       |
C     X        dependent variable array   (input)       |
C     XS       dX/dS array                (calculated)  |
C     N        number of points           (input)       |
C                                                       |
C-------------------------------------------------------
      IF(N.GT.NMAX) STOP 'SPLINEV: array overflow, increase NMAX'
C     
      DO 1 I=2, N-1
        DSM = S(I) - S(I-1)
        DSP = S(I+1) - S(I)
        B(I) = DSP
        A(I) = 2.D0*(DSM+DSP)
        C(I) = DSM
        XS(I) = 3.D0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
    1 CONTINUE
C
C---- set zero second derivative end conditions
      A(1) = 2.D0
      C(1) = 1.D0
      XS(1) = 3.D0*(X(2)-X(1)) / (S(2)-S(1))
      B(N) = 1.D0
      A(N) = 2.D0
      XS(N) = 3.D0*(X(N)-X(N-1)) / (S(N)-S(N-1))
C
C---- solve for derivative array XS
      CALL TRISOL(A,B,C,XS,N)
C
      RETURN
      END ! SPLINEV  
C
C    
C    *************   SUBROUTINE   TRISOL  ************     

      SUBROUTINE TRISOL(A,B,C,D,KK)
      implicit real*8 (a-h,o-z)
      DIMENSION A(KK),B(KK),C(KK),D(KK)
C-----------------------------------------
C     Solves KK long, tri-diagonal system |
C                                         |
C             A C          D              |
C             B A C        D              |
C               B A .      .              |
C                 . . C    .              |
C                   B A    D              |
C                                         |
C     The righthand side D is replaced by |
C     the solution.  A, C are destroyed.  |
C-----------------------------------------
C
      DO 1 K=2, KK
        KM = K-1
        C(KM) = C(KM) / A(KM)
        D(KM) = D(KM) / A(KM)
        A(K) = A(K) - B(K)*C(KM)
        D(K) = D(K) - B(K)*D(KM)
    1 CONTINUE
C
      D(KK) = D(KK)/A(KK)
C
      DO 2 K=KK-1, 1, -1
        D(K) = D(K) - C(K)*D(K+1)
    2 CONTINUE
C
      RETURN
      END ! TRISOL

