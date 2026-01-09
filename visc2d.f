      SUBROUTINE CAVBL2D(ITA)

C**********************************************************************
C     Coupled with 2D Boundary Layer Analysis Code XFOIL              *
C     (Wetted Case)                                 *
C     *
C     By Hong Sun                   April, 2006                   *
C**********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      INCLUDE 'PUFBL.INC'

!s--- YE TIAN 08/20/2013 ------
! --- change following arrays to allocatable
!     DIMENSION XN(NZVIS),YN(NZVIS),XC(NZVIS),YC(NZVIS),UINV(NZVIS)
!    &  ,    XP(NZVIS),YP(NZVIS),DELS(NZVIS)
!    &  ,    FW(NZVIS),ATEMP(NZVIS,NZVIS),BTEMP(NZVIS,NZVIS)
c XM YU 12/2011
c add new variables for viscous run
!    &  ,    CPBLO(NBZ,MBZ),xcold(nzvis),ycold(nzvis),duvlre(nzvis)
!    &  ,    duem(nzvis),uinvsave(nzvis)
c XM YU 12/2011

!     DIMENSION RV1(NZVIS),RRTEMP(NZVIS)
C      DIMENSION DELUE(NZVIS,MBZ)
!     REAL*8    UEDGE, DSTARS, DSTARP, THETAS, THETAP
!     DIMENSION UEDGE(NZVIS), DSTARS(NZVIS), DSTARP(NZVIS)
!    &  ,    THETAS(NZVIS), THETAP(NZVIS)

c XM YU 03/2012
c use new local skin friction for the calculation of propeller forces
c      REAL*8    CDVIS
!     real*8 cfvis(nzvis)
c XM YU 03/2012

!     DIMENSION UXIV(NZVIS),UETAV(NZVIS),UEDGV(NZVIS)


      real,allocatable:: XN(:),YN(:),XC(:),YC(:),UINV(:)
     &  ,    XP(:),YP(:),DELS(:)
     &  ,    FW(:),ATEMP(:,:),BTEMP(:,:)
     &  ,    CPBLO(:,:),xcold(:),ycold(:),duvlre(:)
     &  ,    duem(:),uinvsave(:)

      real,allocatable :: RV1(:),RRTEMP(:)
      real*8,allocatable :: UEDGE(:), DSTARS(:), DSTARP(:)
     &  ,    THETAS(:), THETAP(:)

      real*8,allocatable :: cfvis(:)
      real,allocatable :: UXIV(:),UETAV(:),UEDGV(:)

!e--- YE TIAN 08/20/2013 ------

C.....For Duct
!s--- YE TIAN 08/20/2013 ------
! --- change following arrays to allocatable
!     DIMENSION XND(NZVISD),YND(NZVISD),XCD(NZVISD),YCD(NZVISD)
!    &  ,    XPD(NZVISD),YPD(NZVISD),UINVD(NZVISD),FWD(NZVISD)
!    &  ,    ATEMPD(NZVISD,NZVISD),BTEMPD(NZVISD,NZVISD)
!    &  ,    DELSD(NZVISD)

!     REAL*8    UEDGED,DSTARSD,DSTARPD,THETASD,THETAPD
!     DIMENSION UEDGED(NZVISD),DSTARSD(NZVISD),DSTARPD(NZVISD)
!    &  ,    THETASD(NZVISD),THETAPD(NZVISD)
!     DIMENSION UXIVD(NZVISD),UETAVD(NZVISD),UEDGVD(NZVISD)

      real,allocatable :: XND(:),YND(:),XCD(:),YCD(:)
     &  ,    XPD(:),YPD(:),UINVD(:),FWD(:)
     &  ,    ATEMPD(:,:),BTEMPD(:,:)
     &  ,    DELSD(:)

      real*8,allocatable ::UEDGED(:),DSTARSD(:),DSTARPD(:)
     &  ,    THETASD(:),THETAPD(:)
      real,allocatable ::UXIVD(:),UETAVD(:),UEDGVD(:)
!e--- YE TIAN 08/20/2013 ------

      REAL*8    CDVISD
      INTEGER LVISCON
      INTEGER ITA

!s--- YE TIAN 08/20/2013 ----
! --- allocation

      allocate(XN(NZVIS),YN(NZVIS),XC(NZVIS),YC(NZVIS),UINV(NZVIS)
     &  ,    XP(NZVIS),YP(NZVIS),DELS(NZVIS)
     &  ,    FW(NZVIS),ATEMP(NZVIS,NZVIS),BTEMP(NZVIS,NZVIS)
     &  ,    CPBLO(NBZ,MBZ),xcold(nzvis),ycold(nzvis),duvlre(nzvis)
     &  ,    duem(nzvis),uinvsave(nzvis))

      allocate(RV1(NZVIS),RRTEMP(NZVIS))
      allocate(UEDGE(NZVIS), DSTARS(NZVIS), DSTARP(NZVIS)
     &  ,    THETAS(NZVIS), THETAP(NZVIS))

      allocate(cfvis(nzvis))
      allocate(UXIV(NZVIS),UETAV(NZVIS),UEDGV(NZVIS))

      if (IDUCT .NE. 0) then
        allocate(XND(NZVISD),YND(NZVISD),XCD(NZVISD),YCD(NZVISD)
     &  ,    XPD(NZVISD),YPD(NZVISD),UINVD(NZVISD),FWD(NZVISD)
     &  ,    ATEMPD(NZVISD,NZVISD),BTEMPD(NZVISD,NZVISD)
     &  ,    DELSD(NZVISD))

        allocate(UEDGED(NZVISD),DSTARSD(NZVISD),DSTARPD(NZVISD)
     &  ,    THETASD(NZVISD),THETAPD(NZVISD))
        allocate(UXIVD(NZVISD),UETAVD(NZVISD),UEDGVD(NZVISD))
      end if
!e--- YE TIAN 08/20/2013 ----


C.....IDOPT=1 solve duct, prop separately
      IF(IDUCT.NE.0.AND.IDOPT.EQ.1) GOTO 1111
C
      WRITE(*,*)
      WRITE(*,*)'Start coupling with XFOIL......'
      WRITE(*,*)

      CALL OPFILEVIS

      CALL CLEAR(PPOTW, NSCWZ)
      CALL CLEAR(PPOTWS, NSCWZ)

c      open(1163,file='duvl-test')
C*************LOOP 1000 solves BL variables at each strip**************

      DO 1000 M = 1, MR
c XM YU comment 01/2012
c usd fixed 2d geometry to estimate RE
c        CHDM = (CHORD(M)+CHORD(M+1))/2.0
cC     Reynolds number for each strip
c        IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
c          REC(M) = REYD/2.0     ! based on D
c        ELSE
c          REC(M) = REYD*CHDM/ADVCO*SQRT(ADVCO**2+PI**2*RZP(M)**2)
c        ENDIF
c XM YU comment 01/2012

C     Calculate the inflow angle of attack (pitch angle -inflow angle)

        CALL EVALDKs(NX,1,XR,RZP(M),VAR0,VACUB)
        WROVS=PI*RZP(M)/ADVCO
        AOA = ATAN(0.5*(PITCH(M)+PITCH(M+1))/PI/RZP(M))
     &    - ATAN(VAR0/WROVS)
!s---Allen Du 02/09/2018 add the following so the aoa is not changed in
!rudder case
!       ALPHA =AOA/RAD
        IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
        ELSE
            ALPHA =AOA/RAD
        ENDIF
!e---Allen Du

        NTW(M) = NC+NWSUB-NSUB+NWMIN

C,,,,,,,,,,,,,,,,,,Calculate Geometric Information  ,,,,,,,,,,,,,,,,,,,

c XM YU 11/2011--------------------------------------------------------------
c modify the 2D geometry
c section body
      do n=1,nc+1
         xcold(n)=0.5*(xb(n,m)+xb(n,m+1))
         ymne=0.5*(yb(n,m)+yb(n,m+1))
         zmne=0.5*(zb(n,m)+zb(n,m+1))
         r2c=ymne**2+zmne**2
         rc=sqrt(r2c)
         tangc=atan2(zmne,ymne)
         ycold(n)=rc*tangc
      enddo

      tranx=xcold(nc/2+1)
      trany=ycold(nc/2+1)
      rotang=atan((ycold(nc+1)-ycold(nc/2+1))/(xcold(nc+1)
     &       -xcold(nc/2+1)))
      do n=1,nc+1
         xc(n)=(xcold(n)-tranx)*cos(rotang)+(ycold(n)-trany)*sin(rotang)
         yc(n)=-(xcold(n)-tranx)*sin(rotang)
     &         +(ycold(n)-trany)*cos(rotang)

      enddo
      chdm=(xc(nc+1)-xc(nc/2+1))/2.0
      do n=1,nc+1
         xc(n)=xc(n)-chdm
      enddo


      do n=nc+2,ntw(m)+1
         xc(n)=xivis(n,m)
         yc(n)=etavis(n,m)
      enddo

c write out 2d geo
      if (m.eq.10) then
         write(*,*) 'ntw,nwmin',ntw(m),nwmin
         open(2312,file="2dgeotest")
         do n=1,ntw(m)
           write(2312,*) xc(n),yc(n)
         enddo
      endif
C     Reynolds number for each strip

        IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
!Allen Du 12/22/2017 changed the definition of Re
!         REC(M) = REYD/2.0     ! based on D
          REC(M) = REYD*CHDM    ! based on chordlength
        ELSE
          REC(M) = REYD*CHDM/ADVCO*SQRT(ADVCO**2+PI**2*RZP(M)**2)
        ENDIF

C XM YU 11/2011-------------------------------------------------------------

        DO 20 N = 1, NTW(M)
          XP(N) = (XC(N+1) + XC(N))/2.0
          YP(N) = (YC(N+1) + YC(N))/2.0
          DX1 = XC(N+1) - XC(N)
          DY1 = YC(N+1) - YC(N)
          DS1 = SQRT(DX1*DX1 + DY1*DY1)
          XN(N) = -1.0*DY1/DS1
          YN(N) =  DX1/DS1
          DELTZ(N,M) = DS1
 20     CONTINUE

        DO 30 N = 2, NTW(M)
          CALL ARCLEN( SS1,SS2,XC(N-1),YC(N-1),XC(N),YC(N)
     &      ,XC(N+1),YC(N+1),DELTZ(N-1,M),DELTZ(N,M) )
          IF(N.EQ.2) THEN
            DELTS(1,M)=SS1
          ELSE
            DELTS(N-1,M)=HALF*(SS1+SPRE)
          END IF
          SPRE=SS2
 30     CONTINUE
        DELTS(NTW(M),M)=SS2

        SC(1,M) = 0.0
        DO 40 N = 1, NTW(M)
          SC(N+1,M) = SC(N,M) + DELTS(N,M)
          SP(N,M) = SC(N,M) + HALF*DELTS(N,M)
          DELS(N) = DELTS(N,M)
 40     CONTINUE


C'''''''''''''Calculate Matrices for each strip ''''''''''''''''''''''
C     write(*,*) 'Calculate Matrices for each strip'

        CALL SETUP1(NC,NC,NC,BTEMP,ATEMP,XC,YC,XP,YP,XN,YN,
     &    DELS,NZVIS)
C
C     FW is the wake influence matrix
        DO N = 1, NC
          IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
            FW(N) = - ALPHA*RAD - ATAN2( YP(N) - YC(1), XC(1) - XP(N) )
          ELSE
            FW(N) = - ATAN2( YP(N)-YC(1), XC(1)-XP(N) ) - AOA
          ENDIF
        ENDDO
c XM Yu 10/2011--------------------------------------------------------------
c Replace the 2D coefficients with 3D coefficients
c        DO IM = 1, NC
c          DO IN = 1, NC
c            IF(IM.EQ.1) THEN
c              AV(IN,IM) = ATEMP(IN,IM) - FW(IN)
c            ELSE IF(IM.EQ.NC) THEN
c              AV(IN,IM) = ATEMP(IN,IM) + FW(IN)
c            ELSE
c              AV(IN,IM) = ATEMP(IN,IM)
c            ENDIF
c          END DO
c        END DO                  ! AV

c calculate the update for inviscid velocity (the effects of
c other strips)
        call duvlc(m)

        do im=1,nc
           do in=1,nc
            av(in,im)=avl(in,im,m)
           enddo
        enddo

        CALL INVRSE(AV,AVINV,NC,NZVIS) !AVINV

        NWv = NTW(M)-NC

c        CALL SETUP2(NWv,NC,NC,CS,AWV,XC,YC,XP,YP,XN,YN,
c     &    DELS,NZVIS)           !AWV

c        CALL SETUP1(NC,NTW(M),NC,BETA,ATEMP,XC,YC,XP,YP,XN,YN,
c     &    DELS,NZVIS)           !BETA

        do i=1,nwv
           do j=1,nc
              awv(i,j)=awv3(i,j,m)
           enddo
        enddo

        do im=1,nc
           do in=1,ntw(m)
               beta(im,in)=bv(im,in,m)
           enddo
        enddo

        do i=1,nwv
           do j=1,ntw(m)
              csigma(i,j)=csig(i,j,m)
           enddo
        enddo
c XM Yu 10/2011--------------------------------------------------------------

        CALL MAT2DOT(AVINV,BETA,CA,NC,NZVIS,NC,NZVIS,NTW(M),NZVIS)

        I1=NC/4
        I2=3*I1
        NSTAG=I1
        CPMIN=1.0
        DO I=I1,I2
          CP = -CPBN(I,M)*ADVCO**2. !CPB
          IF(CP.LT.CPMIN)THEN
            NSTAG = I
            CPMIN = CP
          ENDIF
        ENDDO
        CPM = -CPBN(NSTAG-1,M)*ADVCO**2.
        CPP = -CPBN(NSTAG+1,M)*ADVCO**2.
        IF(CPM.LT.CPP) NSTAG=NSTAG-1

        DO 140 I=2,NC
          DO 130 J=1, NTW(M)
            DH(I,J)=(CA(I,J)-CA(I-1,J))/(SP(I,M)-SP(I-1,M))
            IF(I.LT.NSTAG)THEN
              DH(I,J)=-DH(I,J)
            ENDIF
 130      CONTINUE
 140    CONTINUE


c XM YU comment 10/2011
c        CALL SETUP2(NWv,NTW(M),NC,CSIGMA,ATEMP,XC,YC,XP,YP,XN,YN,
c     &    DELS,NZVIS)           !CSIGMA
c XM YU

        CALL MAT2DOT(AWV,AVINV,CB,NWv,NWZ,NC,NBZ,NC,NBZ)

        CALL MAT2DOT(CB,BETA,DH1,NWv,NWZ,NC,NBZ,NTW(M),NZVIS)

        DO 160 I=1,NWv
          DO 150 J=1,NTW(M)
c XM YU 10/2011-----------------------------------------------------------
c refer to 3D influence coefficients formulation
c            DH1(I,J)=(DH1(I,J)+CSIGMA(I,J))/2./PI
            DH1(I,J)=DH1(I,J)+CSIGMA(I,J)
c XM YU 10/2011-----------------------------------------------------------
 150      CONTINUE
 160    CONTINUE

        DO 180 I=2,NWv
          DO 170 J=1,NTW(M)
            DH(NC+I,J)=(DH1(I,J)-DH1(I-1,J))/(SP(NC+I,M)-SP(NC+I-1,M))
 170      CONTINUE
 180    CONTINUE

        DO 190 J=1,NTW(M)
          DH(1,J)   = DH(2,J)+(DH(2,J)-DH(3,J))*(SC(1,M)-SC(2,M))
     &      / (SC(2,M)-SC(3,M))
          DH(NC+1,J)= DH(NC,J)+(DH(NC,J)-DH(NC-1,J))
     &      * (SC(NC+1,M)-SC(NC,M))/(SC(NC,M)-SC(NC-1,M))
          IF(J.EQ.1)   DH(1,J)=-DH(1,J)
          IF(J.EQ.NC)  DH(NC+1,J)=-DH(NC+1,J)
          IF(J.LT.NSTAG) DH(NC+1,J)=-DH(1,J)
          IF(J.GE.NSTAG) DH(1,J)=-DH(NC+1,J)
 190    CONTINUE
c==   %==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%=
c-----------------check symmetry of H matrix-------------------------
C     CALL MAT2DOT(AVINV,AV,BTEMP,NC,NZVIS,NC,NZVIS,NC,NZVIS)
C     IF(M.EQ.1) THEN
C     WRITE(126,*) 'I,J,AV,AVINV,ATEMP,CB,BTEMP'
C     DO I = 1, NC
C     DO J = 1, NC
C     WRITE(126,'(2I6,5F10.4)')I,J,AV(I,J),AVINV(I,J),
C     &       ATEMP(I,J),CB(I,J),BTEMP(I,J)
C     ENDDO
C     ENDDO
C
C     WRITE(127,*) 'I,J,DH1(I,J)'
C     DO I = 1, NWv
C     DO J = 1, NTW(M)
C     IF(J.LE.NC) THEN
C     WRITE(127,'(2I6,4F10.4)') I,J,DH1(I,J),CSIGMA(I,J)
C     &                                  ,AWV(I,J)
C     ELSE
C     WRITE(127,'(2I6,4F10.4)') I,J,DH1(I,J),CSIGMA(I,J),0.
C     ENDIF
C     ENDDO
C     ENDDO
C
C     WRITE(128,*) 'I,J,CA(I,J),beta(i,j)'
C     DO I = 1, NC
C     DO J = 1, NTW(M)
C     WRITE(128,'(2I6,2F10.4)') I,J,CA(I,J),beta(i,j)
C     ENDDO
C     ENDDO
C     ENDIF

c-------------------check symmetry of H matrix------------------------

        DO 210 I=1,NTW(M)-1
          DO 200 J=2,NTW(M)-1
            DD(I,J) = DH(I,J-1)/(SC(J,M)-SC(J-1,M))
     &        - DH(I,J)/(SC(J+1,M)-SC(J,M))
            IF(I.LT.NSTAG.AND.J.LE.NC+1) THEN
              DD(I,J) = -DD(I,J)
            ENDIF
            IF(J.GT.NC+1) DD(1,J)=-DD(1,J)
 200      CONTINUE
 210    CONTINUE

        DO 220 I=1,NTW(M)-1
          IF(I.GT.NC+1)THEN
            DD(I,1)=DD(I,NC+1)
          ELSE
            DD(I,1)=DD(NC+2-I,NC+1)
          ENDIF
 220    CONTINUE

C'''''''''''''''''Calculate UINV  for each strip '''''''''''''''''''''
C     write(*,*)'Calculate UINV  for each strip'

C.......CALCULATE PERTUBATION POTENTIAL IN THE WAKE

        JM = M
        CALL VISPOT(JM)

        DO 50 N = 2, NC
          N1 = N-1
          L  = INDEXB(N,M)
          L1 = INDEXB(N1,M)
          UXI = DPDUB(N,M)+VXIB(N,M)
          UXIM = DPDUB(N1,M)+VXIB(N1,M)
c      VXI1 = VXIB(N,M)
c      VXI2 = VXIB(N1,M)
c      VXIA = (VXI1+VXI2)/2.0
c      UINV(N)= (POT(L)-POT(L1))/(SP(N,M)-SP(N-1,M))-VXIA
          UINV(N)= -(UXI+UXIM)/2.0
          VINFBA = SQRT((VINFSB(N,M)+VINFSB(N1,M))/2.0)
          UINV(N) = UINV(N)/VINFBA
 50     CONTINUE
c         if(m.eq.10) then
c         write(*,*) 'xiangming',vinfsb(10,m),vinfsb(11,m)
c         endif


        write(138,*) 'Zone T="M= ', M, ' " '
        DO 60 IN = 2,NWSUB
          L = NWSUB*(M-1)+IN
          L1 = NWSUB*(M-1)+IN-1
          VXI1 = ABS(VXIWs(IN,M))
          VXI2 = ABS(VXIWs(IN-1,M))
          VXIWA = (VXI1+VXI2)/2.0
          UINV(NC+IN) = (PPOTWs(L)-PPOTWs(L1)) /
     &      (SP(NC+IN,M)-SP(NC+IN-1,M))+VXIWA
          VINFWsA = SQRT((VINFSWs(IN,M)+VINFSWs(IN-1,M))/2.0)
          UINV(NC+IN) = UINV(NC+IN)/VINFWsA
          write(138,*) NC+IN,PPOTWs(L1)
 60     CONTINUE

        write(139,*) 'Zone T="M= ', M, ' " '
        DO 70 IN = NSUB+2, NWMIN
c XM YU 04/13/2012
          L = IDXWAK(IN,M)
          L1 = IDXWAK(IN-1,M)
c          l=(mr-m)*nwmin+in
c          l1=(mr-m)*nwmin+in-1
c          write(*,*)' idx,cal'
c          write(*,*) lin,l
          VXI1 = ABS(VXIW(IN,M))
          VXI2 = ABS(VXIW(IN-1,M))
          VXIWA = (VXI1+VXI2)/2.0
          UINV(NC+IN+NWSUB-NSUB) = (PPOTW(L)-PPOTW(L1)) /
     &      (SP(NC+IN+NWSUB-NSUB,M)-SP(NC+IN-1+NWSUB-NSUB,M))+VXIWA
          VINFWA = SQRT((VINFSW(IN,M)+VINFSW(IN-1,M))/2.0)
          UINV(NC+IN+NWSUB-NSUB) = UINV(NC+IN+NWSUB-NSUB)/VINFWA
          write(139,*) NC+IN+NWSUB-NSUB,PPOTW(L1)
 70     CONTINUE

        L2 = NWSUB*M
        L3 = IDXWAK(NSUB+1,M)
c        l3=(mr-m)*nwmin+nsub+1
        VXI1 = ABS(VXIW(NSUB+1,M))
        VXI2 = ABS(VXIW(NSUB,M))
        VXIWA = (VXI1+VXI2)/2.0
        UINV(NC+NWSUB+1) = (PPOTW(L3)-PPOTWs(L2)) /
     &    (SP(NC+NWSUB+1,M)-SP(NC+NWSUB,M))+VXIWA
        VINFWA = SQRT((VINFSW(NSUB+1,M)+VINFSW(NSUB,M))/2.0)
        UINV(NC+NWSUB+1) =  UINV(NC+NWSUB+1)/VINFWA

c     L2 = IDXWAK(NSUB,M)
c     L3 = IDXWAK(NSUB+1,M)
c     VXI1 = ABS(VXIW(NSUB+1,M))
c     VXI2 = ABS(VXIW(NSUB,M))
c     VXIWA = (VXI1+VXI2)/2.0
c     UINV(NC+NWSUB+1) = (PPOTW(L3)-PPOTW(L2)) /
c     &            (SP(NC+NWSUB+1,M)-SP(NC+NWSUB-NWSUB1,M))+VXIWA
c     VINFWA = SQRT((VINFSW(NSUB+1,M)+VINFSW(NSUB,M))/2.0)
c     UINV(NC+NWSUB+1) =  UINV(NC+NWSUB+1)/VINFWA

C     UINV(1) = UINV(2) + (SC(1,M)-SC(2,M))*
C     &              (UINV(2)-UINV(3))/(SC(2,M)-SC(3,M))
C     UINV(NC+1) = - UINV(1)

        UINV(1) = -UINV(NC+2)
        UINV(NC+1) = UINV(NC+2)

c code test
        if (m.eq.10) then
          open(2111,file='uinv.dat')
          do i=1,ntw(m)
             write(2111,*) xc(i),uinv(i)
          enddo
        endif
C*********************PUT TO VISCAL CONVENTIONS   *****************
C     NON-DIMENSIONALIZATION OF THE MATRICIES BY CHORD/PROPELLER RADIUS
C**********************************************************************
C     NON-DIMENSIONALIZATION

        DO I = 1, NTW(M)
          XC(I)= XC(I)/CHDM/2. + 0.5
          YC(I)= YC(I)/CHDM/2.
          SC(I,M)= SC(I,M)/CHDM/2.
        ENDDO

C     PUT TO VISCAL SIGN CONVENTIONS
        DO 230 I=1,NTW(M)-1
          IF(I.LE.NC+1) THEN
            UVL(I)= UINV(NC+2-I)
            SVL(I)= SC(NC+1,M)-SC(NC+2-I,M)
            XNVIS(I)= XN(NC+2-I)
            YNVIS(I)= YN(NC+2-I)
            XCVIS(I)= XC(NC+2-I)
            YCVIS(I)= YC(NC+2-I)
          ELSE
            UVL(I)= -UINV(I)
            SVL(I)= SC(I,M)
            XNVIS(I)= XN(I)
            YNVIS(I)= YN(I)
            XCVIS(I)= XC(I)
            YCVIS(I)= YC(I)
             ENDIF
 230    CONTINUE

c XM YU 11/2011----------------------------------------------------------
c update the inviscid velocity (effects of other strips)
c test
c        if (m.eq.10) then
c           open(1161,file='duvl')
c           do i=1,ntw(m)-1
c               write(1161,*) duvl(i,m)
c           enddo
c           stop
c        endif
c test

        do i=1,ntw(m)-1
          uinvsave(i)=uvl(i)
          if(i.le.nc/2) then
             uvl(i)=uvl(i)+duvl(i,m)
c          else if (i.le.nc/2) then
c             uvl(i)=uvl(i)
          else
             uvl(i)=uvl(i)-duvl(i,m)
          endif
        enddo
c XM YU 11/2011----------------------------------------------------------

C     NON-DIMENSIONALIZATION
C     CHDVIS = 2.*XCVIS(1)
C     DO I = 1, NTW(M)
C     XCVIS(I)= XCVIS(I)/CHDVIS + 0.5
C     YCVIS(I)= YCVIS(I)/CHDVIS
C     SVL(I)= SVL(I)/CHDVIS
C     DO J=1,NTW(M)
C     DD(I,J)= DD(I,J)*CHDVIS
C     ENDDO
C     ENDDO

C     CALCULATE DVL
        DO 250 I = 1,NTW(M)-1
          DO 240 J = 1,NTW(M)-1
            IF(I.LE.NC+1.AND.J.LE.NC+1)THEN
              DVL(I,J)=DD(NC+2-I,NC+2-J)
            ELSE IF(I.LE.NC+1.AND.J.GT.NC+1)THEN
              DVL(I,J)=DD(NC+2-I,J)
            ELSE IF(I.GT.NC+1.AND.J.LE.NC+1)THEN
              DVL(I,J)=DD(I,NC+2-J)
            ELSE IF(I.GT.NC+1.AND.J.GT.NC+1)THEN
              DVL(I,J)=DD(I,J)
            ENDIF
 240      CONTINUE
 250    CONTINUE

C     PUT TO VISCAL SIGN CONVENTION
        DO 270 I=1,NTW(M)
          DO 260 J=1,NTW(M)
            IF(I.NE.NC+1) DVL(I,J)=-DVL(I,J)
            IF(I.GT.NC+1.AND.J.LE.NC) DVL(I,J)=-DVL(I,J)
            IF(J.GT.NC+1.AND.I.LE.NSTAG) DVL(I,J)=-DVL(I,J)
            IF(J.GT.NC+1.AND.I.EQ.NC+1) DVL(I,J)=-DVL(I,J)
 260      CONTINUE
 270    CONTINUE

        XCVIS(NC+2) = XCVIS(NC+1)
        YCVIS(NC+2) = YCVIS(NC+1)
        SVL(NC+2) = SVL(NC+1)
        UVL(NC+2) = UVL(NC+1)

        DO 280 J=1,NTW(M)
          DVL(NC+2,J)= DVL(NC+1,J)
 280    CONTINUE
        DVL(1,NC+1)=-DVL(1,NC+1)
        DVL(NC+1,NC+1)=-DVL(NC+1,NC+1)
        DVL(NC+2,NC+1)=-DVL(NC+2,NC+1)

        DO 290 J=NSTAG+1,NC
          DVL(1,J)=-DVL(1,J)
          DVL(NC+1,J)=-DVL(NC+1,J)
          DVL(NC+2,J)=-DVL(NC+2,J)
 290    CONTINUE

C-------------------------Code Test -----------------------------
c     IREAD = 0                         !!!
c     if (IREAD.EQ.1) then
c     NWv = 30
c     open(123, file='121.dat')
c     read(123,*)
c     do I = 1, NC+NWv-1
c     read(123,*) IA,UVL(I),SVL(I),XNVIS(I),YNVIS(I)
c     &     ,          XC(I),YC(I)
c     enddo
c     close(123)
c
c     open(124, file='122.dat')
c     DO I = 1, NC+NWv
c     READ(124,293) (DVL(I,J),J=1,NC+NWv)
c     ENDDO
c     close(124)
c     endif
C---------------------------------------------------------------

        WRITE(121,*)'Zone T= "r/R=',RZP(M),', AOA=', AOA/RAD,
     &    ', RE=', REC(M), '"'
        DO I = 1, NC+NWv-1
          WRITE(121,292) I,UVL(I),SVL(I),XNVIS(I),YNVIS(I)
     &      ,          XCVIS(I),YCVIS(I)
        ENDDO
C     DO I = NC+1, 1, -1
C     WRITE(121,*) XCVIS(I),YCVIS(I)
C     ENDDO
 292    FORMAT(1X,I3,7F12.6)

C     WRITE(122,*) 'Zone T= "M = ', M, '"'
C     WRITE(122,*) 'RE = ', REC(M), 'r/R = ', RZP(M)
C     DO I = 1, NC+NWv
C     WRITE(122,293) (DVL(I,J),J=1,NC+NWv)
C     ENDDO
 293    FORMAT(1X,7F10.4)

c==   %==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%=
C     IF(M.EQ.1) THEN
C     WRITE(125,*) 'DH(I,J),DD(I,J),DVL(I,J),DVL(J,I)'
C     DO I = 1, NTW(M)
C     DO J = 1, NTW(M)
C     WRITE(125,221) DH(I,J),DD(I,J),DVL(I,J),DVL(J,I)
C     ENDDO
C     ENDDO
C     221    FORMAT(1X,4F10.4)
C     ENDIF

C     GOTO 1000
C     VISCAL ANALYSIS AT EACH STRIP
C**********************************************************************
        WRITE(10,*) 'Zone T= "r/R = ', RZP(M), '"'
        WRITE(11,*) 'Zone T= "r/R = ', RZP(M), '"'
        WRITE(12,*) 'Zone T= "r/R = ', RZP(M), '"'
C     WRITE(14,*) 'Zone T= "r/R = ', RZP(M), '"'
C     WRITE(97,*) 'Zone T= "r/R = ', RZP(M), '"'

C     WRITE(*,*) 'BEFORE CALL VISCAL ROUTINE'
        NZ = NZVIS
        WRITE(*,*) 'M=', M, 'r/R=',RZP(M)

        XTRANS = XTRANLS(M)
        XTRANP = XTRANLP(M)
        WRITE(*,*) 'XTRANS =', XTRANS, 'XTRANP =', XTRANP
        WRITE(*,*)
c XM YU 10/2011 ---------------------------------------------------------
        do i=1,ntw(m)-1
           do j=1,ntw(m)-1
               dvl(i,j)=-dvl(i,j)
           enddo
        enddo
c XM YU 10/2011 ---------------------------------------------------------

c XM YU 2/2012 change name, cdvis to cfvis
        write(1236,*) 'strip ', M
        CALL VISCAL(NC+1,NWv-2,SVL,XCVIS,YCVIS,XTRANS,XTRANP,REC(M),
     &    ALPHA,UVL,DVL,NZ,UEDGE,DSTARS,DSTARP,THETAS,THETAP,ICV,
     &    ILEP,ITEP,ICLE,ICTE,XNVIS,YNVIS,RVCRIT,MAXIT,EPS1,cfvis,
     &    LVISCON, BETA_BL_INP)
        write(1236,*) 'after call viscal', M, LVISCON
c XM YU 2/2012

        LVFLAG(M) = LVISCON

!Allen Du 12/22/2017 add a file to check if the viscous run converged
        IF(LVISCON==1) THEN
          WRITE(704,'(A2,I3,3X,A1)')"M=",m,"N"
        ELSE
          WRITE(704,'(A2,I3,3X,A1)')"M=",m,"Y"
        ENDIF
!Allen Du 12/22/2017 add a file to check if the viscous run converged

!       CDV(M)= SNGL(CDVIS)
!       cdv(m)=sngl(cfvis(m))
        CDV(M) = 2.0*THETAP(NC+1+NWv-2)
        IF(LVFLAG(M).EQ.1)  CDV(M)=XCDF

C XM YU -- Local skin friction 2/2012
        if (lvflag(m).eq.1) then
          do i=1,nc+1
             cfskin(i,m)=xcdf*vinfba**2
          enddo
        else
          do i=1,nc+1
             cfskin(i,m)=sngl(cfvis(nc+2-i))*vinfba**2
          enddo
        endif
        WRITE(2011,*)'ZONE T="M=', M,' T=', ITSTEP,'"'
C.......Writing wetted pressure for pressure side.
        DO NN=1,NC/2
          WRITE(2011,*) SBP(NC/2-NN+1),cfskin(NN,M)
        ENDDO
C.......Writing wettedpressure for suction side.
        DO NN=1,NC/2
          WRITE(2011,*) SBP(NN),cfskin(NC/2+NN,M)
        ENDDO

c        if (m.eq.1) then
c           open(733,file='cfskin-p')
c           do i=1,nc/2+1
c             write(733,*) xc(i),cfskin(i,m)
c           enddo
c        endif
c XM YU 02/2012

        WRITE(22,'(I5,2F12.8)') M,  RZP(M), CDV(M)

C**********************************************************************
C     CALCULATE 3-D VISCOUS PRESSURE DISTRIBUTION AT PANEL CENTROID
C     BY USING VISCOUS EDGE VELOCITY

C
C     PUT BACK TO PROPELLER COORDINATE SYSTEM

        DO I = 1,NC+1
          UEDGV(I)=SNGL(UEDGE(NC+2-I))
        ENDDO

C     FIND DU CORRECTION AT BLADE T.E.(REFER TO HUFFORD, 1992)

        L1 = INDEXB(1,M)
        LNC = INDEXB(NC,M)
        UXIL = ( ABS(UEDGV(2))+ABS(UEDGV(1)) )/2.
        UXIU = ( ABS(UEDGV(NC+1))+ABS(UEDGV(NC)) )/2.
        UXIL = UXIL*SQRT(VINFSB(1,M)) - VXIB(1,M)
        UXIU = UXIU*SQRT(VINFSB(NC,M)) - VXIB(NC,M)
        ULOWER = ABS(VXIB(1,M)+UXIL)
        UUPPER = ABS(VXIB(NC,M)+UXIU)
        UBLTE = 0.5*(UUPPER+ULOWER)
        UETAL = (DPDVB(1,M)-DPDUB(1,M)*SINPHI(L1))/COSPHI(L1)
        UETAU = (DPDVB(NC,M)-DPDUB(NC,M)*SINPHI(LNC))/COSPHI(LNC)
        WLOWER = VETAB(1,M)+UETAL
        WUPPER = VETAB(NC,M)+UETAU
!       write(*,*) 'visc2d.f:702', VXIB(1,M),VXIB(NC,M),UXIL,UXIU
!       write(*,*) 'visc2d.f:702', ULOWER, UUPPER
        if (abs(UBLTE) .gt. 1e-5) then
          DU = (WLOWER**2-WUPPER**2)/(4.*UBLTE)
        else
          DU = 0.0
        end if

        DO I =1, NC
          L = INDEXB(I,M)
          UXIV(I) = ( ABS(UEDGV(I+1))+ABS(UEDGV(I)) )/2.
          UXIV(I) = UXIV(I)*SQRT(VINFSB(I,M))-VXIB(I,M)
          UETAV(I) = (DPDVB(I,M)-DPDUB(I,M)*SINPHI(L))/COSPHI(L)
          IF(I.EQ.1) THEN
            UXITV = UBLTE - DU
          ELSE IF(I.EQ.NC) THEN
            UXITV = UBLTE + DU
          ELSE
            UXITV = VXIB(I,M)+UXIV(I)
          END IF
          UETATV = VETAB(I,M)+UETAV(I)

          VTOTS_V(L) = UXITV**2 + UETATV**2
          CPBL(I,M) = CPBN(I,M)+VTOTS(L)-VTOTS_V(L)
          CPBLO(I,M) = CPB(I,M)+VTOTS(L)-VTOTS_V(L)

          UXTOT_V(I,M)=UXITV*DIR(L,1,1)+UETATV*DIR(L,2,1)
          UYTOT_V(I,M)=UXITV*DIR(L,1,2)+UETATV*DIR(L,2,2)
          UZTOT_V(I,M)=UXITV*DIR(L,1,3)+UETATV*DIR(L,2,3)

          IF(LVFLAG(M).EQ.1) THEN
            CPBL(I,M)=CPBN(I,M)
            CPBLO(I,M) = CPB(I,M)
            VTOTS_V(L) = VTOTS(L)
            UXTOT_V(I,M)=UXTOT(I,M)
            UYTOT_V(I,M)=UYTOT(I,M)
            UZTOT_V(I,M)=UZTOT(I,M)
          ENDIF
        ENDDO

        WRITE(20,*)'ZONE T="M=', M,' T=', ITSTEP,'"'
!s--- Allen Du 12/22/2017 output the pressure at the control points
        WRITE(7041,*)'ZONE T="M=', M,' T=', ITSTEP,'"'    
C.......Writing wetted pressure for pressure side.
        DO NN=1,NC/2
          LTEMP = indexb(nn,m)   
          write(7041,*) xct(ltemp,1),-CPBLO(NN,M)*ADVCO**2.
          WRITE(20,*) SBP(NC/2-NN+1),-CPBLO(NN,M)*ADVCO**2.
        ENDDO
C.......Writing wettedpressure for suction side.
        DO NN=1,NC/2
          LTEMP = indexb(nn+nc/2,m) 
          write(7041,*) xct(ltemp,1),-CPBLO(NC/2+NN,M)*ADVCO**2.
          WRITE(20,*) SBP(NN),-CPBLO(NC/2+NN,M)*ADVCO**2.
        ENDDO
!e--- Allen Du 12/22/2017 output the pressure at the control points

C/s S.N.KIM - Cp Blade Plotting
        IF(IAN.EQ.2.AND.ITA.EQ.3) THEN
          WRITE(2999+M,*) 'ZONE T= "TIME STEP =',ITSTEP,'"'
          DO NN = 1, NC/2
            L = INDEXB(NN,M)
            WRITE(2999+M,*) XCT(L,1),-CPBLO(NN,M)*ADVCO**2.
          ENDDO
          DO NN = 1, NC/2
            L = INDEXB(NC/2+NN,M)
            WRITE(2999+M,*) XCT(L,1),-CPBLO(NC/2+NN,M)*ADVCO**2.
          ENDDO
        ENDIF
C/e S.N.KIM - Cp Blade Plotting

C**********************************************************************


C
C     COMPUTE MASS DEFECT ALONG EACH STRIP
C

c        if (m.eq.13) then
c        open(1173,file='dstar')
c        do i=1,nc/2
c           write(1173,*) xcvis(i),dstars(i)
c        enddo
c        do i=nc/2+1,ntw(m)-1
c           write(1173,*) xcvis(i),dstarp(i)
c        enddo
c        endif

        DO 310  I = 1, NC/2
          RV1(I) = SNGL(UEDGE(I)*DSTARS(I))*(CHDM*2.0)
 310    CONTINUE
        DO 315  I = NC/2+1, NC+NWv-1
          RV1(I) = SNGL(UEDGE(I)*DSTARP(I))*(CHDM*2.0)
 315    CONTINUE

c XM YU 12/2012 --------------------------------------------------------------
c update the mass defect
        do i=1,ntw(m)-1
           updatem(i,m)=rv1(i)
        enddo
c XM YU 12/2012 --------------------------------------------------------------
c XM YU code test ------------------------------------------
c        sum=0.0
c        do i=1,ntw(m)-1
c          if (i.le.nc/2) then
c             coei=1
c          else
c             coei=-1
c          endif
c          do j=1,ntw(m)-1
c             if (j.le.nc/2) then
c                coej=1
c             else
c                coej=-1
c             endif
c             sum=sum-coei*coej*dvl(i,j)*rv1(j)
c          enddo
c          duem(i)=sum
c          sum=0.0
c         enddo
c
c        if (m.eq.10) then
c        open(1165,file='vel_convergence')
c        do i=1,ntw(m)-1
c         if (i.le.nc/2) then
c          write(1165,*)uedge(i),uinvsave(i),duem(i),duvl(i,m)
c         else
c          write(1165,*)uedge(i),-uinvsave(i),duem(i),duvl(i,m)
c         endif
c        enddo
c        endif

c XM YU code test ---------------------------------------------------

C     PUT BACK TO PROPCAV CONVENTION
        DO 320  I = 1, NC
          RRTEMP(I) = RV1(NC+2-I)
 320    CONTINUE
        DO 325  I = 1, NC
          RV1(I) = RRTEMP(I)
 325    CONTINUE

c     if(M.EQ.5) then
c     do i = 1, nc+nwv+1
c     write(131,*) I, RV1(I)
c     enddo
c     endif
C
C
C     COMPUTE BLOWING SOURCE STRGTH (SIGMA) ALONG EACH STRIP
C
        DO 330  I = 1, NC
          L = INDEXB(I,M)
          BSRCB(L,1) = -(RV1(I+1)-RV1(I))/(SC(I+1,M)-SC(I,M))
          IF(I.GT.NC/2) THEN
            BSRCB(L,1) = - BSRCB(L,1)
          ENDIF

c     if(M.EQ.5) then
c     write(131,*) I, BSRCB(L,1)
c     endif
 330    CONTINUE

        DO 340 I = 1, NWSUB
          L = NWSUB*(M-1)+I
          BSRCWS(L,1)=(RV1(NC+I+1)-RV1(NC+I))/(SC(NC+I+1,M)-SC(NC+I,M))

c     if(M.EQ.5) then
c     write(131,*) I, BSRCWS(L,1)
c     endif
 340    CONTINUE

        DO 350 I = NSUB+1, NWMIN
          L = IDXWAK(I,M)
          BSRCW(L,1)=(RV1(NC+I+NWSUB-NSUB+1)-RV1(NC+I+NWSUB-NSUB))
     &      /(SC(NC+I+NWSUB-NSUB+1,M)-SC(NC+I+NWSUB-NSUB,M))

c     if(M.EQ.5) then
c     write(131,*) I, BSRCW(L,1)
c     endif
 350    CONTINUE


 1000 CONTINUE                  !!     LOOP 1000 ENDED

      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
C     CLOSE(14)
      CLOSE(22)
C     CLOSE(97)

cCRyan Added to ploy x,y,z,viscous pressure wprs_BL.plt
c      OPEN(6911,FILE='wprs_BL.plt',STATUS='UNKNOWN')
c      DO M=1,MR
c         WRITE(6911,5322) M,0
c         DO NN=1,NC/2
c            L = INDEXB(NN,M)
c            RTEMP = SQRT(XCT(L,2)**2+XCT(L,3)**2)
c            WRITE(6911,5335) XCT(L,1),XCT(L,2),XCT(L,3),
c     &               RTEMP,-CPBLO(NN,M)*ADVCO**2
c         END DO
c         DO NN=1,NC/2
c            L = INDEXB(NC/2+NN,M)
c            RTEMP = SQRT(XCT(L,2)**2+XCT(L,3)**2)
c            WRITE(6911,5335) XCT(L,1),XCT(L,2),XCT(L,3),
c     &                    RTEMP,-CPBLO(NC/2+NN,M)*ADVCO**2
c         END DO
c      END DO
c      CLOSE(6911)
c 5335 FORMAT(1X,3(F10.6,1X),F6.3,1X,F12.6)
c 5322 FORMAT('ZONE T="M=',I3,' T=',I3,'"')
cCRyan finished



 1111 CONTINUE


C***********************************************************************
C     INCLUDE DUCT CASE                             *
C***********************************************************************

C     IF(IDUCT.EQ.1) THEN
      IF(IDUCT.EQ.1.AND.IDOPT.EQ.1) THEN ! For Bare Duct case

        WRITE(*,*)
        WRITE(*,*)'Start coupling with XFOIL......DUCT RUN'
        WRITE(*,*)

        CALL OPFILEVIS_D

        DO 2000 M = 1, MDUCT
c     write(*,*) 'M= ', M
C     Reynolds number for each strip
          RECD(M) = REYD/2.0*DCHORD ! based on D
          NTWD(M) = NDUCT+NWSUB-NSUB+NDWK

C     Calculate the inflow angle of attack

C     For Bare Duct case
C     CALL EVALDK(NX,1,XR,YLED,VAR0,VACUB)
C     WROVS=PI*YLED/ADVCO
C     AOA = ATAN(0.5*DUCTPT/PI/YLED)
C     &     - ATAN(VAR0/WROVS)
          AOA = ATAN(0.5*DUCTPT/PI/YLED)
          ALPHAD =AOA/RAD


C,,,,,,,,,,,,,,,,,,Calculate Geometric Information  ,,,,,,,,,,,,,,,,,,,
C     write(*,*) ' Calculate Geometric Information '
          DO 1010 N = 1, NTWD(M)+1
            XCD(N) = XIVISD(N,M)
            YCD(N) = ETAVISD(N,M)
 1010     CONTINUE

          DO 1020 N = 1, NTWD(M)
            XPD(N) = (XCD(N+1) + XCD(N))/2.0
            YPD(N) = (YCD(N+1) + YCD(N))/2.0
            DX1 = XCD(N+1) - XCD(N)
            DY1 = YCD(N+1) - YCD(N)
            DS1 = SQRT(DX1*DX1 + DY1*DY1)
            XND(N) = -1.0*DY1/DS1
            YND(N) =  DX1/DS1
            DELTZD(N,M) = DS1
 1020     CONTINUE

          DO 1030 N = 2, NTWD(M)
            CALL ARCLEN( SS1D,SS2D,XCD(N-1),YCD(N-1),XCD(N),YCD(N)
     &        ,XCD(N+1),YCD(N+1),DELTZD(N-1,M),DELTZD(N,M) )
            IF(N.EQ.2) THEN
              DELTSD(1,M)=SS1D
            ELSE
              DELTSD(N-1,M)=HALF*(SS1D+SPRED)
            END IF
            SPRED=SS2D
 1030     CONTINUE
          DELTSD(NTWD(M),M)=SS2D

          SCD(1,M) = 0.0
          DO 1040 N = 1, NTWD(M)
            SCD(N+1,M) = SCD(N,M) + DELTSD(N,M)
            SPD(N,M) = SCD(N,M) + HALF*DELTSD(N,M)
            DELSD(N) = DELTSD(N,M)
 1040     CONTINUE


C'''''''''''''Calculate Matrices for each strip ''''''''''''''''''''''
C     write(*,*) 'Calculate Matrices for each strip'

          CALL SETUP1(NDUCT,NDUCT,NDUCT,BTEMPD,ATEMPD,XCD,YCD,XPD,YPD,
     &      XND,YND,DELSD,NZVISD)

C     FWD is the wake influence matrix
          DO N = 1, NDUCT
            FWD(N) = -DUCTANG*RAD - ATAN2(YPD(N)-YCD(1), XCD(1)-XPD(N))
     &        -ALPHAD*RAD
          ENDDO

          DO IM = 1, NDUCT
            DO IN = 1, NDUCT
              IF(IM.EQ.1) THEN
                AVD(IN,IM) = ATEMPD(IN,IM) - FWD(IN)
              ELSE IF(IM.EQ.NDUCT) THEN
                AVD(IN,IM) = ATEMPD(IN,IM) + FWD(IN)
              ELSE
                AVD(IN,IM) = ATEMPD(IN,IM)
              ENDIF
            END DO
          END DO                ! AVD

          CALL INVRSE(AVD,AVDINV,NDUCT,NZVISD) !AVDINV

          NWDv = NTWD(M)-NDUCT

          CALL SETUP2(NWDv,NDUCT,NDUCT,CSD,AWVD,XCD,YCD,XPD,YPD,XND,YND,
     &      DELSD,NZVISD)       !AWVD

          CALL SETUP1(NDUCT,NTWD(M),NDUCT,BETAD,ATEMPD,XCD,YCD,XPD,YPD,
     &      XND,YND,DELSD,NZVISD) !BETAD

         CALL MAT2DOT(AVDINV,BETAD,CAD,NDUCT,NZVISD,NDUCT,NZVISD,NTWD(M)
     &       ,NZVISD)


          I1=NDUCT/4
          I2=3*I1
          NSTAGD=I1
          CPMIND=1.0
          DO I=I1,I2
            CP = -CPD(I,M)*ADVCO**2.
            IF(CP.LT.CPMIND)THEN
              NSTAGD = I
              CPMIND = CP
            ENDIF
          ENDDO
          CPM = -CPD(NSTAGD-1,M)*ADVCO**2.
          CPP = -CPD(NSTAGD+1,M)*ADVCO**2.
          IF(CPM.LT.CPP) NSTAGD=NSTAGD-1
          WRITE(*,*)'M=',M, 'NSTAGD=',NSTAGD,'CPMIND=',CPMIND

C     NSTAGD=NDUCT/2+1

          DO 1140 I=2,NDUCT
            DO 1130 J=1, NTWD(M)
              DHD(I,J)=(CAD(I,J)-CAD(I-1,J))/(SPD(I,M)-SPD(I-1,M))
              IF(I.LT.NSTAGD)THEN
                DHD(I,J)=-DHD(I,J)
              ENDIF
 1130       CONTINUE
 1140     CONTINUE

          CALL SETUP2(NWDv,NTWD(M),NDUCT,CSIGMAD,ATEMPD,XCD,YCD,XPD,YPD,
     &      XND,YND,DELSD,NZVISD) !CSIGMAD

          CALL MAT2DOT(AWVD,AVDINV,CBD,NWDv,NZVISD,NDUCT,NZVISD,
     &      NDUCT,NZVISD)

          CALL MAT2DOT(CBD,BETAD,DH1D,NWDv,NZVISD,NDUCT,NZVISD,
     &      NTWD(M),NZVISD)

          DO 1160 I=1,NWDv
            DO 1150 J=1,NTWD(M)
              DH1D(I,J)=(DH1D(I,J)+CSIGMAD(I,J))/2./PI
 1150       CONTINUE
 1160     CONTINUE

          DO 1180 I=2,NWDv
            DO 1170 J=1,NTWD(M)
              DHD(NDUCT+I,J) = (DH1D(I,J)-DH1D(I-1,J))/
     &          (SPD(NDUCT+I,M)-SPD(NDUCT+I-1,M))
 1170       CONTINUE
 1180     CONTINUE

          DO 1190 J=1,NTWD(M)
            DHD(1,J) = DHD(2,J)+(DHD(2,J)-DHD(3,J))*(SCD(1,M)-SCD(2,M))
     &        / (SCD(2,M)-SCD(3,M))
            DHD(NDUCT+1,J) = DHD(NDUCT,J)+(DHD(NDUCT,J)-DHD(NDUCT-1,J))
     &     *(SCD(NDUCT+1,M)-SCD(NDUCT,M))/(SCD(NDUCT,M)-SCD(NDUCT-1,M))
            IF(J.EQ.1)   DHD(1,J)=-DHD(1,J)
            IF(J.EQ.NC)  DHD(NDUCT+1,J)=-DHD(NDUCT+1,J)
            IF(J.LT.NSTAGD) DHD(NDUCT+1,J)=-DHD(1,J)
            IF(J.GE.NSTAGD) DHD(1,J)=-DHD(NDUCT+1,J)
 1190     CONTINUE
c==   %==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%=
c-----------------check symmetry of H matrix-------------------------
C     CALL MAT2DOT(AVDINV,AVD,BTEMPD,NDUCT,NZVISD,NDUCT,NZVISD,
C     &NDUCT,NZVISD)
C     IF(M.EQ.1) THEN
C     WRITE(126,*) 'I,J,AVD,AVDINV,ATEMPD,CBD,BTEMPD'
C     DO I = 1, NDUCT
C     DO J = 1, NDUCT
C     WRITE(126,'(2I6,5F10.4)')I,J,AVD(I,J),AVDINV(I,J),
C     &       ATEMPD(I,J),CBD(I,J),BTEMPD(I,J)
C     ENDDO
C     ENDDO
C
C     WRITE(127,*) 'I,J,DH1D(I,J)'
C     DO I = 1, NWDv
C     DO J = 1, NTWD(M)
C     IF(J.LE.NDUCT) THEN
C     WRITE(127,'(2I6,4F10.4)') I,J,DH1D(I,J),CSIGMAD(I,J),
C     &                                  AWVD(I,J)
C     ELSE
C     WRITE(127,'(2I6,4F10.4)') I,J,DH1D(I,J),CSIGMAD(I,J),0
C     ENDIF
C     ENDDO
C     ENDDO
C
C     WRITE(128,*) 'I,J,CAD(I,J),betaD(i,j)'
C     DO I = 1, NDUCT
C     DO J = 1, NTWD(M)
C     WRITE(128,'(2I6,2F10.4)') I,J,CAD(I,J),betaD(i,j)
C     ENDDO
C     ENDDO
C     ENDIF

c-------------------check symmetry of H matrix------------------------

          DO 1210 I=1,NTWD(M)-1
            DO 1200 J=2,NTWD(M)-1
              DD_D(I,J) = DHD(I,J-1)/(SCD(J,M)-SCD(J-1,M))
     &          - DHD(I,J)/(SCD(J+1,M)-SCD(J,M))
              IF(I.LT.NSTAGD.AND.J.LE.NDUCT+1) THEN
                DD_D(I,J) = -DD_D(I,J)
              ENDIF
              IF(J.GT.NDUCT+1) DD_D(1,J)=-DD_D(1,J)
 1200       CONTINUE
 1210     CONTINUE

          DO 1220 I=1,NTWD(M)-1
            IF(I.GT.NDUCT+1)THEN
              DD_D(I,1)=DD_D(I,NDUCT+1)
            ELSE
              DD_D(I,1)=DD_D(NDUCT+2-I,NDUCT+1)
            ENDIF
 1220     CONTINUE

C'''''''''''''''''Calculate UINVD  for each strip '''''''''''''''''''''
C     write(*,*)'Calculate UINVD  for each strip'

C.......CALCULATE PERTUBATION POTENTIAL IN THE DUCT WAKE

          CALL VISPOT_D(M)

          DO 1050 N = 2, NDUCT
            N1 = N-1
            L  = INDEXD(N,M)
            L1 = INDEXD(N1,M)
            UXI = DPDUD(N,M)+VXID(N,M)
            UXIM = DPDUD(N1,M)+VXID(N1,M)
            VXI1 = VXID(N,M)
            VXI2 = VXID(N1,M)
            VXIA = (VXI1+VXI2)/2.0
            UINVD(N)= (POT(L)-POT(L1))/(SPD(N,M)-SPD(N1,M))-VXIA
C     UINVD(N)= -(UXI+UXIM)/2.0
            VINFBA = SQRT((VINFSD(N,M)+VINFSD(N1,M))/2.0)
            UINVD(N) = UINVD(N)/VINFBA
 1050     CONTINUE

          DO 1060 IN = 2,NWSUB
            L = NWSUB*(M-1)+IN
            L1 = NWSUB*(M-1)+IN-1
            VXI1 = ABS(VXIDWs(IN,M))
            VXI2 = ABS(VXIDWs(IN-1,M))
            VXIWA = (VXI1+VXI2)/2.0
            UINVD(NDUCT+IN) = (PPOTDWs(L)-PPOTDWs(L1)) /
     &        (SPD(NDUCT+IN,M)-SPD(NDUCT+IN-1,M)) + VXIWA
            VINFWsA = SQRT((VINFSDWs(IN,M)+VINFSDWs(IN-1,M))/2.0)
            UINVD(NDUCT+IN) = UINVD(NDUCT+IN)/VINFWsA
 1060     CONTINUE

          DO 1070 IN = NSUB+2, NDWK
            L = INDEXWD(IN,M)
            L1 = INDEXWD(IN-1,M)
            VXI1 = ABS(VXIDW(IN,M))
            VXI2 = ABS(VXIDW(IN-1,M))
            VXIWA = (VXI1+VXI2)/2.0
            IN1 = NDUCT+IN+NWSUB-NSUB
            UINVD(IN1) = (PPOTDW(L)-PPOTDW(L1)) /
     &        (SPD(IN1,M)-SPD(IN1-1,M)) + VXIWA
            VINFWA = SQRT((VINFSDW(IN,M)+VINFSDW(IN-1,M))/2.0)
            UINVD(IN1)=UINVD(IN1)/VINFWA
 1070     CONTINUE

c     L2 = NWSUB*M
c     L3 = INDEXWD(NSUB+1,M)
c     VXI1 = ABS(VXIDW(NSUB+1,M))
c     VXI2 = ABS(VXIDW(NSUB,M))
c     VXIWA = (VXI1+VXI2)/2.0
c     UINVD(NDUCT+NWSUB+1) = (PPOTDW(L3)-PPOTDWs(L2)) /
c     &           (SPD(NDUCT+NWSUB+1,M)-SPD(NDUCT+NWSUB,M))+VXIWA
c     VINFWA = SQRT((VINFSDW(NSUB+1,M)+VINFSDW(NSUB,M))/2.0)
c     UINVD(NDUCT+NWSUB+1) =  UINVD(NDUCT+NWSUB+1)/VINFWA

          L2 = INDEXWD(NSUB,M)
          L3 = INDEXWD(NSUB+1,M)
          VXI1 = ABS(VXIDW(NSUB+1,M))
          VXI2 = ABS(VXIDW(NSUB,M))
          VXIWA = (VXI1+VXI2)/2.0
          UINVD(NDUCT+NWSUB+1) = (PPOTDW(L3)-PPOTDW(L2)) /
     &      (SPD(NDUCT+NWSUB+1,M)-SPD(NDUCT+NWSUB-NWSUB1,M))+VXIWA
          VINFWA = SQRT((VINFSDW(NSUB+1,M)+VINFSDW(NSUB,M))/2.0)
          UINVD(NDUCT+NWSUB+1) =  UINVD(NDUCT+NWSUB+1)/VINFWA

c     UINVD(1) = UINVD(2) + (SCD(1,M)-SCD(2,M))*
c     &              (UINVD(2)-UINVD(3))/(SCD(2,M)-SCD(3,M))
c     UINVD(NDUCT+1) = - UINVD(1)
          UINVD(1) = -UINVD(NDUCT+2)
          UINVD(NDUCT+1) = UINVD(NDUCT+2)

C*********************PUT TO VISCAL CONVENTIONS   *****************
C     NON-DIMENSIONALIZATION OF THE MATRICIES BY CHORD/PROPELLER RADIUS
C**********************************************************************
C     PUT TO VISCAL SIGN CONVENTIONS
C     write(*,*)' PUT TO VISCAL SIGN CONVENTIONS '
          DO 1230 I=1,NTWD(M)-1
            IF(I.LE.NDUCT+1) THEN
              UVLD(I)= UINVD(NDUCT+2-I)
              SVLD(I)= SCD(NDUCT+1,M)-SCD(NDUCT+2-I,M)
              XNVISD(I)= XND(NDUCT+2-I)
              YNVISD(I)= YND(NDUCT+2-I)
              XCVISD(I)= XCD(NDUCT+2-I)
              YCVISD(I)= YCD(NDUCT+2-I)
            ELSE
              UVLD(I)= -UINVD(I)
              SVLD(I)= SCD(I,M)
              XNVISD(I)= XND(I)
              YNVISD(I)= YND(I)
              XCVISD(I)= XCD(I)
              YCVISD(I)= YCD(I)
            ENDIF
 1230     CONTINUE
C     write(*,*) 'DCHORD=', DCHORD
C     NON-DIMENSIONALIZATION
          DO I = 1, NTWD(M)
            XCVISD(I)= XCVISD(I)/DCHORD
            YCVISD(I)= YCVISD(I)/DCHORD
            XCD(I)= XCD(I)/DCHORD
            YCD(I)= YCD(I)/DCHORD
            SVLD(I)= SVLD(I)/DCHORD
C     DO J=1,NTWD(M)
C     DD_D(I,J)= DD_D(I,J)*DCHORD
C     ENDDO
          ENDDO

C     write(*,*)' NON-DIMENSIONALIZATION '
C     CALCULATE DVLD
          DO 1250 I = 1,NTWD(M)-1
            DO 1240 J = 1,NTWD(M)-1
              IF(I.LE.NDUCT+1.AND.J.LE.NDUCT+1)THEN
                DVLD(I,J) = DD_D(NDUCT+2-I,NDUCT+2-J)
              ELSE IF(I.LE.NDUCT+1.AND.J.GT.NDUCT+1)THEN
                DVLD(I,J) = DD_D(NDUCT+2-I,J)
              ELSE IF(I.GT.NDUCT+1.AND.J.LE.NDUCT+1)THEN
                DVLD(I,J) = DD_D(I,NDUCT+2-J)
              ELSE IF(I.GT.NDUCT+1.AND.J.GT.NDUCT+1)THEN
                DVLD(I,J) = DD_D(I,J)
              ENDIF
 1240       CONTINUE
 1250     CONTINUE

C     PUT TO VISCAL SIGN CONVENTION
          DO 1270 I=1,NTWD(M)
            DO 1260 J=1,NTWD(M)
              IF(I.NE.NDUCT+1) DVLD(I,J)=-DVLD(I,J)
              IF(I.GT.NDUCT+1.AND.J.LE.NDUCT) DVLD(I,J)=-DVLD(I,J)
              IF(J.GT.NDUCT+1.AND.I.LE.NSTAGD) DVLD(I,J)=-DVLD(I,J)
              IF(J.GT.NDUCT+1.AND.I.EQ.NDUCT+1) DVLD(I,J)=-DVLD(I,J)
 1260       CONTINUE
 1270     CONTINUE

c     s      XCVISD(NDUCT+2) = XCVISD(NDUCT+1)
c     s      YCVISD(NDUCT+2) = YCVISD(NDUCT+1)
c     s      SVLD(NDUCT+2) = SVLD(NDUCT+1)

          DO 1280 J=1,NTWD(M)
            DVLD(NDUCT+2,J)= DVLD(NDUCT+1,J)
 1280     CONTINUE

          DVLD(1,NDUCT+1)=-DVLD(1,NDUCT+1)
          DVLD(NDUCT+1,NDUCT+1)=-DVLD(NDUCT+1,NDUCT+1)
          DVLD(NDUCT+2,NDUCT+1)=-DVLD(NDUCT+2,NDUCT+1)

          DO 1290 J=NSTAGD+1,NDUCT
            DVLD(1,J)=-DVLD(1,J)
            DVLD(NDUCT+1,J)=-DVLD(NDUCT+1,J)
            DVLD(NDUCT+2,J)=-DVLD(NDUCT+2,J)
 1290     CONTINUE

          WRITE(121,*) 'Zone T= "M(D)= ', M, '"'
          DO I = 1, NDUCT+NWDv-1
            WRITE(121,1292) I,UVLD(I),SVLD(I),XNVISD(I),YNVISD(I)
     &        ,          XCD(I),YCD(I)
          ENDDO
 1292     FORMAT(1X,I3,7F12.6)

C     WRITE(122,*) 'Zone T= "M = ', M, '"'
C     DO I = 1, NDUCT+NWDv
C     WRITE(122,1293) (DVLD(I,J),J=1,NDUCT+NWDv)
C     ENDDO
C     1293  FORMAT(1X,7F10.4)

c==   %==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%=
c     IF(M.EQ.1) THEN
c     WRITE(125,*) 'DHD(I,J),DD_D(I,J),DVLD(I,J),DVLD(J,I)'
c     DO I = 1, NTWD(M)
c     DO J = 1, NTWD(M)
c     WRITE(125,221) DHD(I,J),DD_D(I,J),DVLD(I,J),DVLD(J,I)
c     ENDDO
c     ENDDO
c     221    FORMAT(1X,4F10.4)
c     ENDIF

c     GOTO 2000
C     VISCAL ANALYSIS AT EACH STRIP
C**********************************************************************
          WRITE(10,*) 'Zone T= "M = ', M, '"'
          WRITE(11,*) 'Zone T= "M = ', M, '"'
          WRITE(12,*) 'Zone T= "M = ', M, '"'

C     WRITE(14,*) 'Zone T= "M = ', M, '"'
C     WRITE(97,*) 'Zone T= "M = ', M, '"'

          WRITE(*,*) 'BEFORE CALL VISCAL ROUTINE'
          NZD = NZVISD


          XTRANS = XTRANDS(M)
          XTRANP = XTRANDP(M)

          WRITE(*,*) 'XTRANS =', XTRANS, 'XTRANP =', XTRANP
          CALL VISCAL(NDUCT+1,NWDv-2,SVLD,XCVISD,YCVISD,XTRANS,XTRANP,
     &      RECD(M),ALPHAD,UVLD,DVLD,NZD,UEDGED,DSTARSD,DSTARPD,
     &      THETASD,THETAPD,ICV,ILEP,ITEP,ICLE,ICTE,XNVISD,YNVISD,
     &      RVCRIT,MAXIT,EPS1,CDVISD, LVISCON,BETA_BL_INP)

          LVFLAGD(M) = LVISCON
C
          CDVD(M)=SNGL(CDVISD)
          WRITE(22,*) M,  CDVISD

C**********************************************************************
C     CALCULATE 3-D VISCOUS PRESSURE DISTRIBUTION AT PANEL CENTROID
C     BY USING VISCOUS EDGE VELOCITY

C
C     PUT BACK TO PROPELLER COORDINATE SYSTEM

C     CALCULATE VISCOUS PRESSURE DISTRIBUTION

          DO I = 1,NDUCT+1
            UEDGVD(I)=SNGL(UEDGED(NDUCT+2-I))
          ENDDO

C     FIND DU CORRECTION AT BLADE T.E.(REFER TO HUFFORD, 1992)

          L1 = INDEXD(1,M)
          LNC = INDEXD(NDUCT,M)
          UXIL = ( ABS(UEDGVD(2))+ABS(UEDGVD(1)) )/2.
          UXIU = ( ABS(UEDGVD(NDUCT+1))+ABS(UEDGVD(NDUCT)) )/2.
          UXIL = UXIL*SQRT(VINFSD(1,M)) - VXID(1,M)
          UXIU = UXIU*SQRT(VINFSD(NDUCT,M)) - VXID(NDUCT,M)
          ULOWER = ABS(VXID(1,M)+UXIL)
          UUPPER = ABS(VXID(NDUCT,M)+UXIU)
          UBLTE = 0.5*(UUPPER+ULOWER)
          UETAL = (DPDVD(1,M)-DPDUD(1,M)*SINPHI(L1))/COSPHI(L1)
          UETAU=(DPDVD(NDUCT,M)-DPDUD(NDUCT,M)*SINPHI(LNC))/COSPHI(LNC)
          WLOWER = VETAD(1,M)+UETAL
          WUPPER = VETAD(NDUCT,M)+UETAU
          DU = (WLOWER**2-WUPPER**2)/(4.*UBLTE)

          WRITE(13,*)'ZONE T="M=', M,' T=', ITSTEP,'"'

          DO I =1, NDUCT
            L = INDEXD(I,M)
            UXIVD(I) = ( ABS(UEDGVD(I+1))+ABS(UEDGVD(I)) )/2.
            UXIVD(I) = UXIVD(I)*SQRT(VINFSD(I,M)) - VXID(I,M)
            UETAVD(I) = (DPDVD(I,M)-DPDUD(I,M)*SINPHI(L))/COSPHI(L)
            IF(I.EQ.1) THEN
              UXITV = UBLTE - DU
            ELSE IF(I.EQ.NDUCT) THEN
              UXITV = UBLTE + DU
            ELSE
              UXITV = VXID(I,M)+UXIVD(I)
            END IF
            UETATV = VETAD(I,M)+UETAVD(I)

            VTOTS_V(L) = UXITV**2 + UETATV**2

            CPBLD(I,M)=CPD(I,M)+VTOTS(L)-VTOTS_V(L)

            UXDTOT_V(I,M)=UXITV*DIR(L,1,1)+UETATV*DIR(L,2,1)
            UYDTOT_V(I,M)=UXITV*DIR(L,1,2)+UETATV*DIR(L,2,2)
            UZDTOT_V(I,M)=UXITV*DIR(L,1,3)+UETATV*DIR(L,2,3)
            WRITE(13,*) XPD(I)/DCHORD, -CPBLD(I,M), -CPD(I,M)
          ENDDO

C.......WRITE VISCOUS PRESSURE ON DUCT

          WRITE(20,*)'ZONE T="M=', M,' T=', ITSTEP,'"'

          SUM1 = 0.0
          SUM2 = 0.0

          DO N = 1 , NDUCTH
            N1 = NDUCTH - N + 1
            N2 = NDUCTH + N
            L1 = INDEXD(N1,M)
            L2 = INDEXD(N2,M)

            SUM1 = SUM1 + DELU(L1)
            SUM2 = SUM2 + DELU(L2)
          ENDDO

          DUM = 0.0
          DO N = 1 , NDUCTH
            N1 = NDUCTH - N + 1
            IF(N .EQ. 1) THEN
              L1 = INDEXD(N1,M)
              DUM = DUM + 0.5*DELU(L1) / SUM1
            ELSE
              L1 = INDEXD(N1,M)
              N10 = N1+1
              L10 = INDEXD(N10,M)
              DUM = DUM + 0.5*(DELU(L10) + DELU(L1)) /SUM1
            ENDIF
            WRITE(20,*) DUM, -CPBLD(N1,M)
          ENDDO

          DUM = SUM2
          DO N =  NDUCTH, 1 , -1
            N1 = NDUCTH + N
            IF(N .EQ. NDUCTH) THEN
              L1 = INDEXD(N1,M)
              DUM = DUM - 0.5*DELU(L1) / SUM2
            ELSE
              L1 = INDEXD(N1,M)
              N10 = N1 + 1
              L10 = INDEXD(N10,M)
              DUM = DUM - 0.5*(DELU(L10) + DELU(L1)) /SUM2
            ENDIF
            WRITE(20,*) DUM, -CPBLD(N1,M)
          ENDDO

C**********************************************************************

 2000   CONTINUE

        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
C     CLOSE(14)
        CLOSE(22)
C     CLOSE(97)

      ENDIF                     !IDUCT=1


      RETURN
      END




