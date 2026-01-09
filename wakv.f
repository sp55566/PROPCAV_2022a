      SUBROUTINE wakv
**********************************************************************
c Used to calculate influence coefficents on wake for viscous run
c by XM Yu 01/19/2012 *
***********************************************************************
      use m_WKNP
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
!     COMMON /WKNP/ NWIDX(MBZ),NWSEC(MBZ)
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)
CSH--------------------------------------------
C-----------------------------------------------------------------------
C     Prepare parameters
C-----------------------------------------------------------------------
C.....NWMIN is the number of wake panels at each radius that the........
C.....vortices shedding.................................................
   
      nwv=nwmin+nwsub-nsub
      npwakev=nwv*mr

      do i=1,npanel
         do j=1,npwakev
            bsbw(i,j)=0.0
         enddo
      enddo 
      do i=1,npwakev
         do j=1,npwakev
            bsww(i,j)=0.0
         enddo
      enddo 

      DO 50 M=MR,1,-1

C-----------------------------------------------------------------------
C
C        Prepare panel geometries of     (N,M+1) 2 *-----* 3 (N+1,M+1)
C          the wake, clockwise viewing             |     |
C          from the suction side of the            |     |
C          blade                           (N,M) 1 *-----* 4 (N+1,M)
C
C-----------------------------------------------------------------------

         DO 30 N=1,Nwv
           L=(mr-m)*nwv+n
           XGW(L,1,1)=xwvs3d(N,M)
           XGW(L,1,2)=ywvs3d(N,M)
           XGW(L,1,3)=zwvs3d(N,M)
           XGW(L,2,1)=xwvs3d(N,M+1)
           XGW(L,2,2)=ywvs3d(N,M+1)
           XGW(L,2,3)=zwvs3d(N,M+1)
           XGW(L,3,1)=xwvs3d(N+1,M+1)
           XGW(L,3,2)=ywvs3d(N+1,M+1)
           XGW(L,3,3)=zwvs3d(N+1,M+1)
           XGW(L,4,1)=xwvs3d(N+1,M)
           XGW(L,4,2)=ywvs3d(N+1,M)
           XGW(L,4,3)=zwvs3d(N+1,M)
 30      CONTINUE
 50   continue

      CALL GEO3DW(NPWAKEv,XGW,CHRLEWS,IER)

      IF(IER.EQ.0) THEN
         WRITE(*,*) 'UNACCEPTABLE PANEL IN INFWAK'
         STOP
      END IF
! XM YU 06/2012 
! store the wake panel information
C/s S.N.KIM | added if statement, and 'xctpwv's at other wakes (kk>1)
C           | are considered in the subroutine, 'infwaks.f'.
C/e S.N.KIM | Aug. 2018.
      do j=1,npwakev
         if(ian.eq.2) then
           xctpwv(j,1,1)=xctw(j,1)
           xctpwv(j,2,1)=xctw(j,2)
           xctpwv(j,3,1)=xctw(j,3)
         else
           rcpw=sqrt(xctw(j,2)**2+xctw(j,3)**2)
           thpw=atan2(xctw(j,3),xctw(j,2)) 
           xctpwv(j,1,1)=xctw(j,1)
           xctpwv(j,2,1)=xctw(j,2)
           xctpwv(j,3,1)=xctw(j,3)
           if (nblade.gt.1) then
              do kk=2,nblade
                 xctpwv(j,1,kk)=xctw(j,1)
                 xctpwv(j,2,kk)=rcpw*cos(thpw+delk*(kk-1))
                 xctpwv(j,3,kk)=rcpw*sin(thpw+delk*(kk-1))
              enddo
           endif
         endif
      enddo
! XM YU 06/2012 
    
      DO 290 m=mr,1,-1
         DO 270 L1=1,nwv
            L=(mr-m)*nwv+l1

C.....Transfer data to the common block /GEOM/..........................
            DO 110 K=1,4
               XV(K)=XVPW(L,K)
               YV(K)=YVPW(L,K)
               SIDE(K)=SIDW(L,K)
110         CONTINUE
            DO 130 K=1,15
               S(K)=SSW(L,K)
130         CONTINUE
C
C.....XM1,XM2,XM3,XM4 for Morino's formulation..........................
            XM1(1)=xwvs3d(L1,M)
            XM1(2)=ywvs3d(L1,M)
            XM1(3)=zwvs3d(L1,M)
            XM2(1)=xwvs3d(L1,M+1)
            XM2(2)=ywvs3d(L1,M+1)
            XM2(3)=zwvs3d(L1,M+1)
            XM3(1)=xwvs3d(L1+1,M+1)
            XM3(2)=ywvs3d(L1+1,M+1)
            XM3(3)=zwvs3d(L1+1,M+1)
            XM4(1)=xwvs3d(L1+1,M)             
            XM4(2)=ywvs3d(L1+1,M)             
            XM4(3)=zwvs3d(L1+1,M) 

            do 230 kk=1,nblade
C.....Add the influence coefficients from all panels....................
               DO 210 I=1,npanel
                  
C.....Transfer control points to the local coordinate...................
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 170 K=1,3
                     XLOC=XLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,1,K)
                     YLOC=YLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,KK)-XCTW(L,K))*DIRW(L,3,K)
 170              CONTINUE

C.....Computethe induced potentials.....................................
                  IMR=1
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     *                 FSX,FSY,FDX,FDY,FDZ,0,IMR)
                  
                  bsbw(I,L)=bsbw(i,L)+fs
c IMAGE MODEL
                  if ((icon.eq.5).and.(wingimag.ne.0)) then

                     call rudimv_bw(L,i,kk,imr,s_add)
                     bsbw(I,L)=fs+s_add
c                  else
c                     bsbw(I,L)=fs
                  endif
c IMAGE MODEL

                  if (kk.eq.1) then
                      kt=nc*(mr-m)+1
                      if(i.ge.kt.and.i.lt.kt+nc) then
                         ii=i+1-kt
                         bv(ii,nc+l1,m)=fs
                      endif
                  endif
 210           CONTINUE
 230        continue

            kk=1
            do 95 i1=1,nwv
                  i=(mr-m)*nwv+i1
C.................Transfer control points to the local coordinate
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 75 K=1,3
                     XLOC=XLOC+(XCTW(I,K)-XCTW(L,K))*DIRW(L,1,K)
                     YLOC=YLOC+(XCTW(I,K)-XCTW(L,K))*DIRW(L,2,K)
                     ZLOC=ZLOC+(XCTW(I,K)-XCTW(L,K))*DIRW(L,3,K)
75                CONTINUE
C
C.................Compute the induced potentials
ccc .. since only source influence needed it is unnecessary to use
ccc .. hyperboloid panels near te edge.  Answers are almost identical
ccc .. for source influence but are not for dipole influence.
                  IMR=1
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     1                      FSX,FSY,FDX,FDY,FDZ,0,IMR)
                  csig(i1,nc+l1,m)=fs/(4.*pi)
95          continue

          do 102 kk=1,nblade
           do 103 mm=mr,1,-1
            do 105 i1=1,nwv
                  i=(mr-mm)*nwv+i1
C.................Transfer control points to the local coordinate
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 77 K=1,3
                     XLOC=XLOC+(XCTpWv(I,K,kk)-XCTW(L,K))*DIRW(L,1,K)
                     YLOC=YLOC+(XCTpWv(I,K,kk)-XCTW(L,K))*DIRW(L,2,K)
                     ZLOC=ZLOC+(XCTpWv(I,K,kk)-XCTW(L,K))*DIRW(L,3,K)
77                CONTINUE
C
C.................Compute the induced potentials
ccc .. since only source influence needed it is unnecessary to use
ccc .. hyperboloid panels near te edge.  Answers are almost identical
ccc .. for source influence but are not for dipole influence.
                  IMR=1
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEWS(L),FS,FD,
     1                      FSX,FSY,FDX,FDY,FDZ,0,IMR)
                  bsww(I,L)=bsww(I,L)+fs
c IMAGE MODEL
                  if ((icon.eq.5).and.(wingimag.ne.0)) then

                     call rudimv_ww(L,i,kk,imr,s_add)
                     bsww(I,L)=fs+s_add
c                  else
c                     bsww(I,L)=fs
                  endif
c IMAGE MODEL

105          continue
103         continue
102        continue

 270     CONTINUE
 290  CONTINUE

      RETURN
      END
