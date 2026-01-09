      SUBROUTINE WAKEPLT
************************************************************************
*                                                                      *
*  Date     Revision or Comment                                        *
*  -------- -------------------                                        *
*  CM041798 This subroutine was created (moved from propcav)           *
*                                                                      *
*                                                                      *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C---------This Section is for looking at the wake at each timestep------
C---------of the cavitating solution to make sure that everything-------

      OPEN(800,FILE="wake.plt",STATUS='UNKNOWN')

      WRITE(800,*) 'TITLE="Plot Wake Geometry"'
      WRITE(800,*) 'VARIABLES="x","y","z"'

      IF(NSUB.NE.0) THEN
         IF(NSUB.EQ.1.and.NWSUB1.EQ.1) THEN
           WRITE(800,5305) (MR+1)*(NWSUB1*NSUB+1),(MR)*(NWSUB1*NSUB)
         ELSE
           WRITE(800,5305) (MR+1)*NWSUB1*NSUB,(MR)*(NWSUB1*NSUB-1)
         ENDIF
 5305    FORMAT(1x,'ZONE T="Subpanels" N=',I5,' E=',I5,
     *          ' F=FEPOINT ET=QUADRILATERAL')
         DO 226 M=1,MR+1
           IF(NSUB.EQ.1.and.NWSUB1.EQ.1) THEN
            DO 777 N123=1, NWSUB1*NSUB+1
               L123 = (MR-M+1)*NWSUB1*NSUB+N123
               XWS1 = XWS(N123,M)
               YWS1 = YWS(N123,M)
               ZWS1 = ZWS(N123,M)
               IF(ISP.NE.0) CALL ROTATE2(0,SPANGLE,XWS1,YWS1)
               WRITE(800,*) XWS1,YWS1,ZWS(N123,M)
 777        CONTINUE
           ELSE
            DO 225 N123=1, NWSUB1*NSUB
               L123=(MR-M+1)*NWSUB1*NSUB+N123
               XWS1 = XWS(N123,M)
               YWS1 = YWS(N123,M)
               ZWS1 = ZWS(N123,M)
               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XWS1,YWS1)
               WRITE(800,*) XWS1,YWS1,ZWS(N123,M)
 225        CONTINUE
           ENDIF
 226     CONTINUE

         DO 228 M=1,MR
          IF(NSUB.EQ.1.and.NWSUB1.EQ.1) THEN
            DO 888 N123=1,NWSUB1*NSUB
               WRITE(800,*) (M-1)*NWSUB1*(NSUB+1)+N123,(M-1)*NWSUB1
     *              *(NSUB+1)+N123+1,
     *              NWSUB1*(NSUB+1)*(M)+N123+1,NWSUB1*(NSUB+1)*(M)+N123
 888        CONTINUE

          ELSE
            DO 227 N123=1,NWSUB1*NSUB-1
               WRITE(800,*) (M-1)*NWSUB1*NSUB+N123,(M-1)*NWSUB1
     *              *NSUB+N123+1,
     *              NWSUB1*NSUB*(M)+N123+1,NWSUB1*NSUB*(M)+N123
 227        CONTINUE
          ENDIF
 228     CONTINUE
      END IF

C----------Section added to look at wake before it is passed to geo3dw--
C----------will plot it in wake.plt-----------------------------CM010998
C----------NNTSUM is number of points, NANSUM is number of panels.------
C----------these numbers are required for tecploting--------------------
      DO 11 M=1,MR
         IF(NSW(M).EQ.NSW(M+1))THEN
            NWPAN(M)=NSW(M)-1
         ELSE
            NWPAN(M)=MIN(NSW(M),NSW(M+1))
            NWPAN(M)=NWPAN(M)-1
         ENDIF
 11   CONTINUE


      nntsum=0
      nansum=0
      do 15 m=1,mr
         nansum=nansum+nwpan(m)
         nntsum=nntsum+(nwpan(m)+1)*2
 15   continue

C-----Modified routine to also plot the wake of the other blades.
      DELK=TWOPI/NBLADE
 5001 format('ZONE T="wake',i2,'",N=',I5,' E=',I5,' F=FEPOINT',
     *     ' ET=QUADRILATERAL')
      do k=1,nblade

         write(800,5001) k,nntsum,nansum

         T=-DELK*FLOAT(K-1)

         DO M=1,MR+1
            DO N=1,NSW(M)
               XW1(N,M) = XW(N,M)
               YW1(N,M)=YW(N,M)*COS(T)-ZW(N,M)*SIN(T)
               ZW1(N,M)=YW(N,M)*SIN(T)+ZW(N,M)*COS(T)
               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XW1(N,M),YW1(N,M))
            END DO
         END DO

C---------Add points to file (adding only quadrilaterals to wake). A lot
C---------of points will be redundant for simplicity's sake.--CM010998--

         do 17 m=1,mr
            do 16, n=1,nwpan(m)+1
               write(800,*) xw1(n,m),yw1(n,m),zw1(n,m)
               write(800,*) xw1(n,m+1),yw1(n,m+1),zw1(n,m+1)
 16         continue
 17      continue

C---------Add connectivities to file---------------------------CM010998-

         nncount=-2
         do 19 m=1,mr
            nncount=nncount+2
            if(m.gt.1)then
               nncount=nncount+2*(nwpan(m-1))
            endif
            do 18 n=1,nwpan(m)
               write(800,*) 1+(n-1)*2+nncount, 2+(n-1)*2+nncount,
     *              4+(n-1)*2+nncount,3+(n-1)*2+nncount
 18         continue
 19      continue

      end do

      CLOSE(800)

      RETURN

C---------End of Wake plot Section!-------------------------------------
      END



      SUBROUTINE WAKEPLT_CHECK
************************************************************************
*                                                                      *
*  Date     Revision or Comment                                        *
*  -------- -------------------                                        *
*  CM041798 This subroutine was created (moved from propcav)           *
*                                                                      *
*                                                                      *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C---------This Section is for looking at the wake at each timestep------
C---------of the cavitating solution to make sure that everything-------


!S.N.Kim - See if Sub Panel Numer is well chageed into 1.
!      IF(IVISC.NE.1) THEN
!       WRITE(*,*) 'Subpanel ='
!       WRITE(*,*) NSUB
!       WRITE(*,*) 'Subpanel1 ='
!       WRITE(*,*) NWSUB1
!      END IF
!S.N.Kim

      IF(NSUB.NE.0) THEN
         IF(NSUB.EQ.1.and.NWSUB1.EQ.1) THEN
          WRITE(3367,5305)(MR+1)*(NWSUB1*NSUB+1),(MR)*(NWSUB1*NSUB)
         ELSE
           WRITE(3367,5305) (MR+1)*NWSUB1*NSUB,(MR)*(NWSUB1*NSUB-1)
         ENDIF
 5305    FORMAT(1x,'ZONE T="Subpanels" N=',I5,' E=',I5,
     *          ' F=FEPOINT ET=QUADRILATERAL')
         WRITE(3367,*) 'SOLUTIONTIME= ',icavt
         DO 226 M=1,MR+1
           IF(NSUB.EQ.1.and.NWSUB1.EQ.1) THEN
            DO 777 N123=1, NWSUB1*NSUB+1
               L123 = (MR-M+1)*NWSUB1*NSUB+N123
               XWS1 = XWS(N123,M)
               YWS1 = YWS(N123,M)
               ZWS1 = ZWS(N123,M)
               IF(ISP.NE.0) CALL ROTATE2(0,SPANGLE,XWS1,YWS1)
               WRITE(3367,*) XWS1,YWS1,ZWS(N123,M)
 777        CONTINUE
           ELSE
            DO 225 N123=1, NWSUB1*NSUB
               L123=(MR-M+1)*NWSUB1*NSUB+N123
               XWS1 = XWS(N123,M)
               YWS1 = YWS(N123,M)
               ZWS1 = ZWS(N123,M)
               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XWS1,YWS1)
               WRITE(3367,*) XWS1,YWS1,ZWS(N123,M)
 225        CONTINUE
           ENDIF
 226     CONTINUE

         DO 228 M=1,MR
          IF(NSUB.EQ.1.and.NWSUB1.EQ.1) THEN
            DO 888 N123=1,NWSUB1*NSUB
               WRITE(3367,*)(M-1)*NWSUB1*(NSUB+1)+N123,(M-1)*NWSUB1
     *              *(NSUB+1)+N123+1,
     *              NWSUB1*(NSUB+1)*(M)+N123+1,NWSUB1*(NSUB+1)*(M)+N123
 888        CONTINUE

          ELSE
            DO 227 N123=1,NWSUB1*NSUB-1
               WRITE(3367,*) (M-1)*NWSUB1*NSUB+N123,(M-1)*NWSUB1
     *              *NSUB+N123+1,
     *              NWSUB1*NSUB*(M)+N123+1,NWSUB1*NSUB*(M)+N123
 227        CONTINUE
          ENDIF
 228     CONTINUE
      END IF

C----------Section added to look at wake before it is passed to geo3dw--
C----------will plot it in wake.plt-----------------------------CM010998
C----------NNTSUM is number of points, NANSUM is number of panels.------
C----------these numbers are required for tecploting--------------------
      DO 11 M=1,MR
         IF(NSW(M).EQ.NSW(M+1))THEN
            NWPAN(M)=NSW(M)-1
         ELSE
            NWPAN(M)=MIN(NSW(M),NSW(M+1))
            NWPAN(M)=NWPAN(M)-1
         ENDIF
 11   CONTINUE


      nntsum=0
      nansum=0
      do 15 m=1,mr
         nansum=nansum+nwpan(m)
         nntsum=nntsum+(nwpan(m)+1)*2
 15   continue

C-----Modified routine to also plot the wake of the other blades.
      DELK=TWOPI/NBLADE
 5001 format('ZONE T="wake',i2,'",N=',I5,' E=',I5,' F=FEPOINT',
     *     ' ET=QUADRILATERAL')
      do k=1,nblade

         write(3367,5001) k,nntsum,nansum
         WRITE(3367,*) 'SOLUTIONTIME= ',icavt

         T=-DELK*FLOAT(K-1)

         DO M=1,MR+1
            DO N=1,NSW(M)
               XW1(N,M) = XW(N,M)
               YW1(N,M)=YW(N,M)*COS(T)-ZW(N,M)*SIN(T)
               ZW1(N,M)=YW(N,M)*SIN(T)+ZW(N,M)*COS(T)
               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XW1(N,M),YW1(N,M))
            END DO
         END DO

C---------Add points to file (adding only quadrilaterals to wake). A lot
C---------of points will be redundant for simplicity's sake.--CM010998--

         do 17 m=1,mr
            do 16, n=1,nwpan(m)+1
               write(3367,*) xw1(n,m),yw1(n,m),zw1(n,m)
               write(3367,*) xw1(n,m+1),yw1(n,m+1),zw1(n,m+1)
 16         continue
 17      continue

C---------Add connectivities to file---------------------------CM010998-

         nncount=-2
         do 19 m=1,mr
            nncount=nncount+2
            if(m.gt.1)then
               nncount=nncount+2*(nwpan(m-1))
            endif
            do 18 n=1,nwpan(m)
               write(3367,*) 1+(n-1)*2+nncount, 2+(n-1)*2+nncount,
     *              4+(n-1)*2+nncount,3+(n-1)*2+nncount
 18         continue
 19      continue
      end do

      RETURN

C---------End of Wake plot Section!-------------------------------------
      END

      SUBROUTINE WAKEPLT_UNSTEADY
************************************************************************
*                                                                      *
*  Date     Revision or Comment                                        *
*  -------- -------------------                                        *
*  SNKIM021617 This subroutine is provided to visualize unsteady       *
*              movement of shedding vortex                             *
*                                                                      *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      INTEGER I, J, K, N, II

      ALLOCATABLE :: FCOFMAT(:)

      ALLOCATE(FCOFMAT(3000))

      DO N = 1, NBLADE
C --- Wake Geometry Plotting.
        WRITE(8829,*) 'ZONE DATAPACKING=BLOCK I=', nwpanel+1, ',J=', mrp
        WRITE(8829,*) 'SOLUTIONTIME=', IUNSPLOT
        DO I = 1, MRP
          DO J = 1, NWPANEL + 1
            WRITE(8829,*) XWW(J,I,N)
          ENDDO
        ENDDO
        DO I = 1, MRP
          DO J = 1, NWPANEL + 1
            WRITE(8829,*) YWW(J,I,N)
          ENDDO
        ENDDO
        DO I = 1, MRP
          DO J = 1, NWPANEL + 1
            WRITE(8829,*) ZWW(J,I,N)
          ENDDO
        ENDDO
C --- End of Wake Geometry Plotting.

C --- Wake Strength Plotting.
        IREC = NTPOS(N)
c        IF (N.EQ.1.AND.ICAVT.EQ.1) THEN
c          IREC = IREC - 1
c          IF (IREC.EQ.0) IREC = IREC + (360 / NDLTAT)
c        ENDIF
        IF (IUNSPLOT.EQ.1) THEN 
          TEMP5 = - 0.1 
        ELSE
          CALL READ2(46,IREC,TEMP5,NWPANEL*MR)
        ENDIF
        DO I = 1, MR + 1
          DO J = 1, NWPANEL + 1
            I35 = INDEXW2(J,I)
            WRITE(8829,*) -TEMP5(I35)
          ENDDO
        ENDDO
C --- End of Wake Strength Plotting.

C --- Wake Inf. Coef. Plotting.
        IO = 80 + N
        REWIND IO
        DO I = 1, NWMIN + 1
          DO J = MR + 1, 1, -1
            II = (MR+1-J)*(NWMIN+1)+I
            IF (I.EQ.NWMIN+1.OR.J.EQ.MR+1) THEN 
              TEMP9 = 0.0
              GO TO 1010
            ENDIF
            CALL READ1(IO,TEMP1,NPANEL)
            TEMP9 = 0.0
            DO K = 1, NPANEL
              TEMP9 = TEMP9 + TEMP1(K)
            ENDDO             
 1010 CONTINUE
            FCOFMAT(II) = TEMP9
          ENDDO
        ENDDO
        DO I = 1, MR + 1
          DO J = 1, NWMIN + 1
            II = (MR+1-I)*(NWMIN+1)+J
            WRITE(8829,*) FCOFMAT(II)
          ENDDO
        ENDDO
C --- End of Wake Inf. Coef. Plotting.
      ENDDO

      IUNSPLOT = IUNSPLOT + 1

      RETURN
C---------End of Wake plot Section!-------------------------------------
      END SUBROUTINE WAKEPLT_UNSTEADY

