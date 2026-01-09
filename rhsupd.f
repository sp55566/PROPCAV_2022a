      SUBROUTINE RHSUPD(RRHS)
************************************************************************
*     RHSUPD: Right-Hand-Side UPDate                                   *
*      --- Compute current right-hand-side                             *
*  Date of last Revision         Revision                              *
*  ---------------------         --------                              *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION RRHS(*)

C-----------------------------------------------------------------------
C     Initialization
C-----------------------------------------------------------------------
      CALL CLEAR(RRHS,NTZ)
CSH--Calculate the influence coeff from the image, so NBLADE=2
            IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=2
CSH--------------------------------------------           
      REWIND 41 

      IDDUCT2 = NPANB + NPANH 
      IDDUCT3 = NPANB + NPANH + NPAND

      IF(NTSTEP.EQ.0)GO TO 155
      IF(NBLADE.EQ.1)GO TO 105

C-----------------------------------------------------------------------
C     Compute dipole influence from the other blades
C-----------------------------------------------------------------------
      DO 20 KK=2,NBLADE
C.......NOTE: the angular position of the blade and wake................
         IREC=NTPOS(KK)
         CALL READ2(45,IREC,STRGTH(1,KK),NPANEL)
 20   CONTINUE

      do m = mr, 1, -1
         do n = 1 , nc
            j = indexb(n,m)
            call read1(41,temp1,npanel)
            do kk = 2, nblade
             call read1(41,temp1,npanel)
             do i = 1 , npanel
                  RRHS(I)=RRHS(I)-TEMP1(I)*STRGTH(J,KK)
             enddo
            enddo
         enddo
      enddo
      
      if(ihub .ne. 0) then
         do n = 1, nhbx 
          do m = 1 , mhbt
             j = indexh(n,m)
             call read1(41,temp1,npanel)
             do kk = 2, nblade
                  call read1(41,temp1,npanel)
                  do i = 1 , npanel
                     RRHS(I)=RRHS(I)-TEMP1(I)*STRGTH(J,KK)
                  enddo
             enddo
            enddo
         enddo
      endif 

      if(iduct .ne. 0) then
         do m = 1, mduct 
          do n = 1 , nduct
             j = indexd(n,m)
             call read1(41,temp1,npanel)
             do kk = 2, nblade
                  call read1(41,temp1,npanel)
                  do i = 1 , npanel
                     RRHS(I)=RRHS(I)-TEMP1(I)*STRGTH(J,KK)
                  enddo
             enddo
            enddo
         enddo
      endif 

      if(itun .ne. 0) then
         do n = 1, naxt 
          do m = 1 , mtunel
             j = indextn(n,m)
             call read1(41,temp1,npanel)
             do kk = 2, nblade
                  call read1(41,temp1,npanel)
                  do i = 1 , npanel
                     RRHS(I)=RRHS(I)-TEMP1(I)*STRGTH(J,KK)
                  enddo
             enddo
            enddo
         enddo
      endif 
      
C-----------------------------------------------------------------------
C     Compute wake influence from the other blades
C-----------------------------------------------------------------------
C.....TEMP5: here is DPhi(k)............................................
C.....TEMP1 is the wake infl. coef.'s...................................

      DO 100 KK=2,NBLADE
         IO=80+KK
         REWIND IO
         IREC=NTPOS(KK)
         NREAD=NWMIN*MR
         CALL READ2(46,IREC,TEMP5,NREAD)
         DO 90 L=1,NWMIN
            DO 80 M=MR,1,-1
               
               idx = indexw2(L,m)
               if(l .lt. nwmin) then
                  idx1 = indexw2(L+1,m)
               endif
               
               CALL READ1(IO,TEMP1,NPANEL)
               
               DO 75 I=1,NPANEL
                  if(L .lt. nwmin) then
                     RRHS(I)=RRHS(I)
     %                      -TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
                  else
                     RRHS(I)=RRHS(I)-TEMP1(I)*TEMP5(IDX)
                  endif
 75            CONTINUE
 80         CONTINUE
 90      CONTINUE
 100  CONTINUE
      
      IF(IDUCT .NE. 0) THEN
         DO KK=2,NBLADE
            IO=500+KK
            REWIND IO
            IREC=NTPOS(KK)
            CALL READ2(49,IREC,TEMPD4,NPWAKED)
            DO L=1,NDWK
               DO M=1, MDUCT
                  
                  IDX = INDEXWD(L,M)
                  IF(L .LT. NDWK) THEN
C                     IDX1 = INDEXWD(N+1,M)    ! It's a bug. --S.KIM 2021.
                     IDX1 = INDEXWD(L+1,M)
                  ENDIF
                  
                  CALL READ1(IO,TEMP1,NPANEL)
                  
                  DO I=1,NPANEL
                     IF(L .LT. NDWK) THEN
                        RRHS(I)=RRHS(I)
     %                       -TEMP1(I)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
                     ELSE
                        RRHS(I)=RRHS(I)-TEMP1(I)*TEMPD4(IDX)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

C-----------------------------------------------------------------------
C     Compute the influence of shedding vortices of the key blade
C-----------------------------------------------------------------------
 105  CONTINUE

      IREC=NTPOS(1)-1
      IF(IREC.LE.0)IREC=IREC+TWOPI/DELTAT
      
C.... Read last time step's shed vorticities in the wake (TEMP5).........
      NREAD=NWMIN*MR
      CALL READ2(46,IREC,TEMP5,NREAD)

C....Read the key blade wake inf. func. (TEMP1).........................
      REWIND 81

C....L=1 (first wake panel of the key blade)............................
      DO 106 M=MR,1,-1
         
C......IDX=M: t(k-1) solutions..........................................
         IDX = indexw2(1,m)
         CALL READ1(81,TEMP1,NPANEL)
         DO 102 I=1,NPANEL
            RRHS(I)=RRHS(I)-(HALF*TEMP1(I)-WSUBIF(I,M))*TEMP5(IDX)
 102     CONTINUE
 106  CONTINUE

C.....L>1...............................................................
      DO 130 L=2,NWMIN
         DO 120 M=MR,1,-1

C........IDX: t(k-1);  IDX1: t(k-2).....................................
            idx = indexw2(L-1,m)
            idx1 = indexw2(L,m)
            CALL READ1(81,TEMP1,NPANEL)
            
            DO 110 I=1,NPANEL
               RRHS(I)=RRHS(I)-TEMP1(I)*HALF*(TEMP5(IDX)+TEMP5(IDX1))
 110        CONTINUE
 120     CONTINUE
 130  CONTINUE

C-----------------------------------------------------------------------
C     Add the influence of the rest of the wake
C      (assume with steady strength)
C-----------------------------------------------------------------------

      DO 150 I=1,NPANEL
         DO 140 M=MR,1,-1
            RRHS(I)=RRHS(I)-WSTINF(I,M)*DPHI(M,0)
 140     CONTINUE
 150  CONTINUE


      IF(IDUCT .NE. 0) THEN
         CALL READ2(49,IREC,TEMPD4,NPWAKED)

         REWIND 501

C....L=1 (first wake panel of the key duct)............................

         DO M=1, MDUCT
            IDX = INDEXWD(1,M)
            CALL READ1(501,TEMP1,NPANEL)
            DO I=1,NPANEL
               RRHS(I)=RRHS(I)-(HALF*TEMP1(I)-WSUBIFD(I,M))*TEMPD4(IDX)
            ENDDO
         ENDDO

C.....L>1...............................................................
         DO L=2,NDWK
            DO M=1, MDUCT

               IDX = INDEXWD(L-1,M)
               IDX1 = INDEXWD(L,M)

               CALL READ1(501,TEMP1,NPANEL)
            
               DO I=1,NPANEL
               RRHS(I)=RRHS(I)-TEMP1(I)*HALF*(TEMPD4(IDX)+TEMPD4(IDX1))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      
      
C-----------------------------------------------------------------------
C     Compute source influence from all the blades
C-----------------------------------------------------------------------
 155  CONTINUE

C.....call inflow, and update STRGTH as the strength of sources.........

      IF(ITERGB .EQ. 0) THEN
         CALL INFLOW
      ELSEIF(ITERGB .EQ. 1 .AND. ICON .EQ. 5) THEN
         CALL INFLOWH
      ENDIF

      REWIND 42
      
      DO M = MR , 1, -1
         DO N = 1 , NC
            J = INDEXB(N,M)
            DO KK = 1, NBLADE
               CALL READ1(42,TEMP1,NPANEL)
               DO I = 1 , NPANEL
                  RRHS(I)=RRHS(I)+TEMP1(I)*STRGTH(J,KK)
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      IF(IHUB .NE. 0) THEN
         DO N = 1 , NHBX
            DO M = 1, MHBT
               J = INDEXH(N,M)
               DO KK = 1 , NBLADE
                  CALL READ1(42,TEMP1,NPANEL)
                  DO I = 1 , NPANEL
                     RRHS(I)=RRHS(I)+TEMP1(I)*STRGTH(J,KK)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF(IDUCT .NE. 0) THEN
         DO M = 1 , MDUCT
            DO N = 1, NDUCT
               J = INDEXD(N,M)
               DO KK = 1 , NBLADE
                  CALL READ1(42,TEMP1,NPANEL)
                  DO I = 1 , NPANEL
                     RRHS(I)=RRHS(I)+TEMP1(I)*STRGTH(J,KK)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF(ITUN .NE. 0) THEN
         DO N = 1 , NAXT
            DO M = 1, MTUNEL
               J = INDEXTN(N,M)
               DO KK = 1 , NBLADE
                  CALL READ1(42,TEMP1,NPANEL)
                  DO I = 1 , NPANEL
                     RRHS(I)=RRHS(I)+TEMP1(I)*STRGTH(J,KK)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.        
c      IF(IAN .EQ. 2) THEN
c         KK = 1
c         DO N = 1 , NTHX
c            DO M = 1, MCVT
c               J = INDEXT(N,M)
c               CALL READ1(42,TEMP1,NPANEL)
c               DO I = 1 , NPANEL
c                  RRHS(I)=RRHS(I)+TEMP1(I)*STRGTH(J,KK)
c               ENDDO
c            ENDDO
c         ENDDO
c
c         DO N = 1 , NCVX
c            DO M = 1, MCVT
c               J = INDEXC(N,M)
c                  CALL READ1(42,TEMP1,NPANEL)
c                  DO I = 1 , NPANEL
c                     RRHS(I)=RRHS(I)+TEMP1(I)*STRGTH(J,KK)
c                  ENDDO
c            ENDDO
c         ENDDO
c      ENDIF
C/e S.N.KIM | Aug. 2018.

      IF((ICON.EQ.5).AND.(IHULL.EQ.1)) NBLADE=1

      RETURN
      END 
