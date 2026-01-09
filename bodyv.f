      SUBROUTINE bodyv
************************************************************************
c This subroutine is used to calculate the influnece coefficients on body
c for  viscous case
c By XM YU 01/19/2012
c last version 6/2012
***********************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION XM1(3),XM2(3),XM3(3),XM4(3),XMC(3)
!s--YE TIAN  -- 06/20/2013---
      IMR0=1
!e--YE TIAN  -- 06/20/2013---

! XM YU 6/2012
      do i=1,npanel
         do j=1,npanel
            bsbb(i,j)=0.0
         enddo
      enddo

      npwakev=nwv*mr
      do i=1,npwakev
         do j=1,npanel
            bswb(i,j)=0.0
         enddo
      enddo
! XM YU 6/2012

      DO 160 M=MR,1,-1
         DO 150 N=1,NC
            J=INDEXB(N,M)
            DO 100 K=1,4
               XV(K)=XVP(J,K)
               YV(K)=YVP(J,K)
               SIDE(K)=SID(J,K)
 100        CONTINUE
            DO 110 K=1,15
               S(K)=SS(J,K)
 110        CONTINUE

            XM1(1)=XB(N,M)
            XM1(2)=YB(N,M)
            XM1(3)=ZB(N,M)
            XM2(1)=XB(N,M+1)
            XM2(2)=YB(N,M+1)
            XM2(3)=ZB(N,M+1)
            XM3(1)=XB(N+1,M+1)
            XM3(2)=YB(N+1,M+1)
            XM3(3)=ZB(N+1,M+1)
            XM4(1)=XB(N+1,M)
            XM4(2)=YB(N+1,M)
            XM4(3)=ZB(N+1,M)   
            
             do 140 kk=1,nblade
               DO 130 I=1,npanel
                  XLOC=ZERO
                  YLOC=ZERO
                  ZLOC=ZERO
                  DO 120 K=1,3
                     XLOC=XLOC+(XCTP(I,K,kk)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(XCTP(I,K,kk)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(XCTP(I,K,kk)-XCT(J,K))*DIR(J,3,K)
 120              CONTINUE
                  
c                  IMR=IMR0
                  IMR=1
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),
     *                 FS,FD,FSX,FSY,FDX,FDY,FDZ,0,IMR)
                  
                  bsbb(i,j)=bsbb(i,j)+fs
c IMAGE MODEL it is assumed only one blade when image model is applied.
                  if (icon.eq.5.and.wingimag.ne.0) then

                      call rudimv_bb(j,i,kk,imr,s_add)
                      bsbb(i,j)=fs+s_add
c                  else
c                      bsbb(i,j)=fs
                  endif
c IMAGE MODEL
                  if (kk.eq.1) then
                      kt=nc*(mr-m)+1
                      if(i.ge.kt.and.i.lt.kt+nc) then
                        if(j.ge.kt.and.j.lt.kt+nc) then
                           ii=i+1-kt
                           jj=j+1-kt 
                           bv(ii,jj,m)=fs
                        endif
                      endif
                  endif

 130           CONTINUE
 140           continue 
      
ccc......compute viscous influence (zone 2)
            do 340 i=1,nwv
               ii=(mr-m)*nwv+i
                kk=1
cC.................Transfer control points to local coordinate
                  XLOC=0.
                  YLOC=0.
                  ZLOC=0.
                  DO 320 K=1,3
                     XLOC=XLOC+(xctw(ii,K)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(xctw(ii,K)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(xctw(ii,K)-XCT(J,K))*DIR(J,3,K)

320               CONTINUE
C.................Compute the induced potentials due to the blades
                  IMR=IMR0
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                      FDX,FDY,FDZ,0,IMR)
C
C.................Near field use Morino's formulation
                  IF(IMR.EQ.2) THEN
                     DO 332 IXYZ=1,3
                        XMC(IXYZ)=xctw(ii,IXYZ)
332                  CONTINUE
                     CALL HYPOT(XM1,XM2,XM3,XM4,XMC,FD,FS)
                  END IF
C.................Write down all the inf. func's (except self inf.)
C.................  which are > 6.27
                     IF(ABS(FD).GT.6.27) THEN
                        FD=0.0
                     END IF
ccc .......... save needed influence functions
                csig(i,n,m)=fs/(4.*pi)
                awv3(i,n,m)=fd/(4.*pi)
340        continue

         do 372 kk=1,nblade
           do 371 mm=mr,1,-1
            do 370 i=1,nwv
               ii=(mr-mm)*nwv+i
cC.................Transfer control points to local coordinate
                  XLOC=0.
                  YLOC=0.
                  ZLOC=0.
                  DO 390 K=1,3
                     XLOC=XLOC+(xctpwv(ii,K,kk)-XCT(J,K))*DIR(J,1,K)
                     YLOC=YLOC+(xctpwv(ii,K,kk)-XCT(J,K))*DIR(J,2,K)
                     ZLOC=ZLOC+(xctpwv(ii,K,kk)-XCT(J,K))*DIR(J,3,K)

390               CONTINUE
C.................Compute the induced potentials due to the blades
                  IMR=1
                  CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),FS,FD,FSX,FSY,
     *                      FDX,FDY,FDZ,0,IMR)
                     bswb(ii,j)=bswb(ii,j)+fs
c IMAGE MODEL
                  if (icon.eq.5.and.wingimag.ne.0) then

                     call rudimv_wb(j,ii,kk,imr,s_add)
                     bswb(ii,j)=fs+s_add
c                  else
c                     bswb(ii,j)=fs
                  endif
c IMAGE MODEL
 370        continue
 371       continue
 372      continue

 150     CONTINUE
 160  CONTINUE

      RETURN
      END
