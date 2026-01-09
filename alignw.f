C     -------------------------------------------------
      SUBROUTINE ALIGNW
C     -------------------------------------------------
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      
C     DIMENSION DUMV(MRMAX),DUMX(MRMAX),
C     % DUMY(MRMAX),DUMZ(MRMAX)
C      DIMENSION DUMV1(MRMAXP),DUMX1(MRMAXP),
C     % DUMY1(MRMAXP),DUMZ1(MRMAXP)
C      DIMENSION DUMVX(4*MRMAX),DUMVY(4*MRMAX),
C     % DUMVZ(4*MRMAX)
C      DIMENSION VXYZ(NWMAX,MRMAXP,3)
C      DIMENSION XW2(NWMAXP,MRMAXP),YW2(NWMAXP,MRMAXP),
C     % ZW2(NWMAXP,MRMAXP)

CVV
      ALLOCATABLE :: DUMV(:),DUMX(:),DUMY(:),DUMZ(:)
      ALLOCATABLE :: DUMV1(:),DUMX1(:),DUMY1(:),DUMZ1(:)
      ALLOCATABLE :: DUMVX(:),DUMVY(:),DUMVZ(:)
      ALLOCATABLE :: VXYZ(:,:,:)
      ALLOCATABLE :: XW2(:,:),YW2(:,:)
      ALLOCATABLE :: ZW2(:,:)
CVV

CVV
      ALLOCATE(DUMV(MRMAX),DUMX(MRMAX),DUMY(MRMAX),DUMZ(MRMAX))
      ALLOCATE( DUMV1(MRMAX),DUMX1(MRMAX),DUMY1(MRMAX),DUMZ1(MRMAX))
      ALLOCATE( DUMVX(4*MRMAX),DUMVY(4*MRMAX),DUMVZ(4*MRMAX))
      ALLOCATE( VXYZ(NWMAX,MRMAXP,3))
      ALLOCATE( XW2(NWMAXP,MRMAXP),YW2(NWMAXP,MRMAXP))
      ALLOCATE( ZW2(NWMAXP,MRMAXP))
CVV

      DTPROP=DELTAT/RAD
      DTF=DTPROP*ADVCO/180.0
      DTI=DTF
      NWPP= NWPANEL+1

      CALL CLEAR(VXYZ,MRMAXP*NWMAX*3)

C     -- Calculate total velocity at the center of free vortex grid 
C     on wake surface
C     
      DO N = 1 , NWPANEL
       XX1 = HALF*(XVC(N+1,1) + XVC(N,1))
       YY1 = HALF*(YVC(N+1,1) + YVC(N,1))
       ZZ1 = HALF*(ZVC(N+1,1) + ZVC(N,1))

       XX2 = 0.0
       YY2 = 0.0
       ZZ2 = 0.0

       DO M = 1 , MCVT
        XX2 = XX2 + HALF*(XVC(N+1,M)+XVC(N,M))/REAL(MCVT)
        YY2 = YY2 + HALF*(YVC(N+1,M)+YVC(N,M))/REAL(MCVT)
        ZZ2 = ZZ2 + HALF*(ZVC(N+1,M)+ZVC(N,M))/REAL(MCVT)
       ENDDO
       
       DDX1 = XX2 - XX1
       DDY1 = YY2 - YY1
       DDZ1 = ZZ2 - ZZ1
       RRCC = SQRT(DDX1**2+DDY1**2+DDZ1**2)

       DUMV1(1) = 0.0

       DO M = 1 , MR
        IF(M .EQ. 1) THEN
         DV1 = 0.0
        ELSE
         DX1 = 0.25*(XW(N+1,M)+XW(N,M)
     %    -XW(N+1,M-1)-XW(N,M-1))
         DY1 = 0.25*(YW(N+1,M)+YW(N,M)
     %    -YW(N+1,M-1)-YW(N,M-1))
         DZ1 = 0.25*(ZW(N+1,M)+ZW(N,M)
     %    -ZW(N+1,M-1)-ZW(N,M-1))
         DV1 = SQRT(DX1**2+DY1**2+DZ1**2)
        ENDIF
        DX2 = 0.25*(XW(N+1,M+1)+XW(N,M+1)
     %   -XW(N+1,M)-XW(N,M))
        DY2 = 0.25*(YW(N+1,M+1)+YW(N,M+1)
     %   -YW(N+1,M)-YW(N,M))
        DZ2 = 0.25*(ZW(N+1,M+1)+ZW(N,M+1)
     %   -ZW(N+1,M)-ZW(N,M))
        DV2 = SQRT(DX2**2+DY2**2+DZ2**2)
        
        IF(M .EQ. 1) THEN
         DUMV(M) = DV1+DV2
         DUMV1(M+1) = 2.*DV2
        ELSE
         DUMV(M) = DUMV(M-1)+ DV1+DV2
         DUMV1(M+1) = DUMV1(M) + 2.*DV2
        ENDIF

       ENDDO

       DUMV(MR) = DUMV(MR) + DV2 + RRCC

       DO M = 1 , MR
        IF(M .EQ. MR) THEN
         DUMX(M) = VXMEAN(N)
         DUMY(M) = VYMEAN(N)
         DUMZ(M) = VZMEAN(N)
        ELSE
         I1 = INDEXW2(N,M)
         DUMX(M) = XIND(I1)
         DUMY(M) = YIND(I1)
         DUMZ(M) = ZIND(I1)
        ENDIF
       ENDDO

       CALL UGLYDK(MR,1,1,DUMV,DUMX,0.0,0.0,DUMVX)
       CALL UGLYDK(MR,1,1,DUMV,DUMY,0.0,0.0,DUMVY)
       CALL UGLYDK(MR,1,1,DUMV,DUMZ,0.0,0.0,DUMVZ)

       CALL EVALDK(MR,MRP,DUMV,DUMV1,DUMX1,DUMVX)
       CALL EVALDK(MR,MRP,DUMV,DUMV1,DUMY1,DUMVY)
       CALL EVALDK(MR,MRP,DUMV,DUMV1,DUMZ1,DUMVZ)

       DO M = 1 , MRP
        VXYZ(N,M,1) = DUMX1(M)
        VXYZ(N,M,2) = DUMY1(M)
        VXYZ(N,M,3) = DUMZ1(M)
       ENDDO

      ENDDO             

C     -- Calculate new coorninate of wake vertex points

      IST = NWPANEL-4
      NITER = 4*ICAVT

      NN1 = 1
      NSWW = IST

      IF(NTSTEP .EQ. 0 .AND. NITER .LT. IST) THEN
       NN1 = NITER - 3
       NSWW = NITER
      ENDIF

      DO N = 1 , NWPP
       DO M = 1 , MRP
        XW2(N,M) = XW(N,M)
        YW2(N,M) = YW(N,M)
        ZW2(N,M) = ZW(N,M)
       ENDDO
      ENDDO

      IF(IHUB .NE. 0) NHBUB = NHBU + NH

      DO M = MRP,1,-1
       DO N = NN1, NSWW
        IF(NTSTEP .EQ.0 .OR. IHUB .EQ. 0 
     %   .OR. IHUB .EQ. 1 .or. ihub .eq. 2 
     %   .or. IHUB .eq. 4) 
     %   XW2(N+1,M) = XW2(N,M)+abs(VXYZ(N,M,1))*DTI
        YW2(N+1,M) = YW2(N,M) + VXYZ(N,M,2)*DTI
        ZW2(N+1,M) = ZW2(N,M) + VXYZ(N,M,3)*DTI
       ENDDO
      ENDDO

      IF(IHUB .EQ. 0) GO TO 777

C     --- Correction for m = 1

      MRH1 = MR / 3
      MRH = MR - MRH1
      MRH1P = MRH1 + 1
      MRHP = MRH + 1

      DO N = NN1 , NSWW

       DO M = 1 , MRHP
        DUMV(M) = SQRT(YW2(N+1,M+MRH1)**2+ZW2(N+1,M+MRH1)**2)
        DUMY(M) = DANGLE(ZW2(N+1,M+MRH1),YW2(N+1,M+MRH1))
        IF(M .NE. 1) THEN
         DDD = DUMY(M-1) - DUMY(M)
         IF(ABS(DDD) .GT. PI) THEN
          IF(DDD .LT. 0.0) THEN
           DUMY(M-1) = DUMY(M-1) + 2.*PI
          ELSE
           DUMY(M) = DUMY(M) + 2.*PI
          ENDIF                
         ENDIF
        ENDIF

       ENDDO

       IF(IHUB .EQ. 0 .or. IHUB .EQ. 4) THEN
        RR1 = RHUB
       ELSEIF(IHUB .EQ. 1) THEN
        RR1 = RHUB
        IF(NTSTEP .NE. 0) THEN
         RR1 = SQRT(YW2(N+1,1)**2 + ZW2(N+1,1)**2)
        ENDIF
        
       ELSEIF(IHUB .EQ. 2) THEN
        RR1 = RHUB
       ELSEIF(IHUB .EQ. 5) THEN
        RNOS = RNOSET( (XHBFD-XW2(N+1,1))/XHBT, RHUB,RHULT)
        RR1 = RNOS
        IF(RNOS .LT. RHULT) RR1 = RHULT
        IF(XW2(N+1,1) .GT. XHBFD) RR1 = RHULT
       ENDIF

       RR2 = SQRT(YW2(N+1,MRH1P)**2+ZW2(N+1,MRH1P)**2)
       DRR = (RR2 - RR1) / REAL(MRH1) 

       
C     CALL BQUAD(DUMV,DUMY,RR1,IM,COA,COB,COC,
C     %                           THE,MRHP)

       DTHE= DUMY(2) - DUMY(1)
       DRRR = DUMV(2) - DUMV(1)

       SLP = DTHE / DRRR

       THE = SLP * ( RR1 - DUMV(1))+ DUMY(1)

       YW2(N+1,1) = RR1 * COS(THE)
       ZW2(N+1,1) = RR1 * SIN(THE)

       DO M = MRHP+1,1,-1
        IF(M .EQ. 1) THEN
         DUMV(M) = RR1
         DUMY(M) = THE
        ELSE
         DUMV(M) = DUMV(M-1)
         DUMY(M) = DUMY(M-1)
        ENDIF
       ENDDO

       DO M = 1, MRH1-1

        RRR = RR1 + DRR * M

        CALL BQUAD(DUMV,DUMY,RRR,IM,COA,COB,COC,
     %   THE,MRHP+1)

        THE = SLP * ( RRR - DUMV(1))+ DUMY(1)

        YW2(N+1,M+1) = RRR * COS(THE)
        ZW2(N+1,M+1) = RRR * SIN(THE)

       ENDDO      

      ENDDO

 777  CONTINUE

      IF(IHUB .EQ. 0 .OR. IHUB .EQ. 1 .or. 
     % IHUB .eq. 2 .or. IHUB .EQ. 4) THEN
       DO M = MRP,1,-1
        RR1 = SQRT(YW2(NSWW+1,M)**2+ZW2(NSWW+1,M)**2)
        RR0 = SQRT(YW2(NSWW,M)**2+ZW2(NSWW,M)**2)
        DDX = XW2(NSWW+1,M) - XW2(NSWW,M)
        
        DDRR = RR1 - RR0

        THE = DANGLE(ZW2(NSWW+1,M),YW2(NSWW+1,M))
        THEOLD = DANGLE(ZW2(NSWW,M),YW2(NSWW,M))

        DTII = DELTAT

        IF(THEOLD .GT. THE) THEN
C     DTII = -DELTAT
         IF(THE - THEOLD .LT. -PI) DTII = DELTAT
        ENDIF
        
        XX0 = DDX
        DO N = NSWW+2, NWPANEL+1
         THE = THE + DTII
         XW2(N,M) = XW2(N-1,M) + DTI

         XX0 = XX0 + DTI

         RR2 = RR1
         
         YW2(N,M) = RR2 * COS(THE)
         ZW2(N,M) = RR2 * SIN(THE)
        ENDDO
       ENDDO
      ELSEIF(IHUB .EQ. 5) THEN

       DO M = MRP,1,-1
        THE = DANGLE(ZW2(NSWW+1,M), YW2(NSWW+1,M))
        THEOLD = DANGLE(ZW2(NSWW,M),YW2(NSWW,M))

        DTII = DELTAT

        IF(THEOLD .GT. THE) THEN
C     DTII = -DELTAT
         IF(THE - THEOLD .LT. -PI) DTII = DELTAT
        ENDIF

        
        DO N = NSWW+2, NWPP

         THE = THE + DTII

         XW2(N,M) = XW2(N-1,M) + DTI             
         IF(M .EQ. 1) THEN
          RNOS = RNOSET((XHBFD-XW2(N,1))/XHBT, RHUB,RHULT)
          IF(RNOS .LT. RHULT) THEN
           RR1 = RHULT
          ELSE
           RR1 = RNOS
          ENDIF
          IF(XW2(N,1) .GT. XHBFD) RR1 = RHULT
          if(nsww .eq. ist) then
           RR1 = RHULT
          endif
         ELSE
          RR1 = SQRT(YW2(N-1,M)**2+ZW2(N-1,M)**2)
         ENDIF
         
         YW2(N,M) = RR1 * COS(THE)
         ZW2(N,M) = RR1 * SIN(THE)
        ENDDO
       ENDDO

      ENDIF

 1111 do m = 1 , mrp
       do n = 1 , nwpanel+1
        xw(n,m) = xw2(n,m)
        yw(n,m) = yw2(n,m)
        zw(n,m) = zw2(n,m)
       enddo
      enddo

!8888 call gwakes
 8888 call gwakes(5,2)


      IF(IHUB .EQ. 1) THEN         
       CALL GHUB1
      ELSEIF(IHUB .EQ. 2) THEN         
       CALL GHUB2
      ELSEIF(IHUB .EQ. 4) THEN         
       CALL GHUB4
      ELSEIF(IHUB .EQ. 5) THEN
       CALL GHUB5
      ENDIF

CVV
      DEALLOCATE( DUMV,DUMX,DUMY,DUMZ )
      DEALLOCATE( DUMV1,DUMX1,DUMY1,DUMZ1 )
      DEALLOCATE( DUMVX,DUMVY,DUMVZ )
      DEALLOCATE( VXYZ )
      DEALLOCATE( XW2,YW2 )
      DEALLOCATE( ZW2 )
CVV
      return
      end
      
