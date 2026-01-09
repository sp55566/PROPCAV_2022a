C     THE FOLLOWING SUBROUTINES ARE COPIED FROM CAV2DBL BY HONG SUN


      SUBROUTINE SETUP1(N1,nw,nt,B,C,XC,YC,XP,YP,XN,YN,DZ,nz)

C------------------------- MIT - PAN2D ---------------------------------
C     CALCULATES INFLUENCE COEFFICIENTS & SETS-UP SIMULTANEOUS EQUATIONS
C-----------------------------------------------------------------------   
      DIMENSION XC(nz),YC(nz),XP(nz),YP(nz),XN(nz),YN(nz),DZ(nz) 
      DIMENSION B(NZ,NZ),C(NZ,NZ)

      CPLOC=0.5

      DO 30 N=1,N1                                                       
      DO 20 M=1,NW 
         ICOD=0                                                            
         IF(M.EQ.N) ICOD=1                                                 
         CALL POTCON2(XP(N),YP(N), XC(M),YC(M),XC(M+1),YC(M+1),XN(M),
     *               YN(M),DZ(M),B(N,M),C(N,M),ICOD,CPLOC)
20       CONTINUE
30    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE SETUP2(nw,n2,nt,b,c,xc,yc,xp,yp,xn,yn,dz,nz)
c--------------------------------------------------------------------
      dimension xc(nz),yc(nz),xp(nz),yp(nz),xn(nz),yn(nz),dz(nz)
      dimension b(nz,nz),c(nz,nz)

      cploc = 0.5
      do 1 i=1,nw
            ni=nt+i
            do 2 j=1,n2
                  icod = 0
                  if(nt+i.eq.j)icod=1
                  call potcon2(xp(ni),yp(ni),xc(j),yc(j),xc(j+1)
     *,yc(j+1),xn(j),yn(j),dz(j),b(i,j),c(i,j),icod,cploc)                    
 2            continue
 1      continue
      return
      end
C
C      
C
      SUBROUTINE POTCON2(XP,YP, XL,YL,XR,YR,XN,YN,D,B,C,ICOD,CPL)
C-----------------------------------------------------------------------
C        INFLUENCE COEFFICIENT AT (XP,YP)
C                 BY A CONSTANT STRENGTH SEGMENT (XL,YL)-(XR,YR)
C                 OF LENGTH D WITH NORMAL (XN,YN)
C         B     ; INDUCED POTENTIAL DUE TO SOURCE
C         C     ; INDUCED POTENTIAL DUE TO NORMAL DIPOLE
C        (UD,VD); INDUCED VELOCITY BY NORMAL DIPOLE  ! never used
C        (US,VS); INDUCED VELOCITY BY SOURCE    ! never used
c         if iflag = false compute source and dipole influence functions
c         if iflag = true compute cphi (a velocity influenc function
c                    defined as integral along panel of d2/dn**2 (log r)
c                    the coefficient comes back to setup in c (ie. c = cphi)
C----------------------------------------------------------------------

      DATA PI/3.1415962653589793/

      IF(ICOD.EQ.1) GO TO 20   ! n=m 
      DX=XP-XL
      DY=YP-YL
      R1S=DX*DX+DY*DY
      R1=SQRT(R1S)
      TX=DX*YN-DY*XN
      TY=DX*XN+DY*YN
      DX2=XP-XR
      DY2=YP-YR
      R2S=DX2*DX2+DY2*DY2
      R2=SQRT( R2S)
      PH1=ATAN2(TX,TY)
      PH2=ATAN2(TX-D,TY)
      DPH=PH2-PH1                 
      IF( ABS(DPH).GT.PI) DPH=DPH+SIGN(2*PI,PH1)   ! dph is included angle
      B=-TX*LOG(R2/R1)+D*(LOG(R2)-1.)-TY*DPH
      C=DPH

      GO TO 90
 20   XD=CPL*D
      DMX=D-XD
      B=XD*LOG(XD)+DMX*LOG(DMX)-D
      C=PI
      go to 90 

 90   RETURN
      END
