      SUBROUTINE PROPLT
C************************ MIT PSF10.02 PROPLT **************************
C
C                PROPELLER STEADY FLOW ANALYSIS PROGRAM
C                        PANEL METHOD SOLUTION
C
C
C       Copyright (c) Massachusetts Institute of Technology 1989
C
C                     Release Date: 31 January 1989
C
C***********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION XB1(NBPZ,MBPZ),YB1(NBPZ,MBPZ),ZB1(NBPZ,MBPZ),
     %          XH1(NHPZ,MHPZ),YH1(NHPZ,MHPZ),ZH1(NHPZ,MHPZ)
      dimension dum1(3),dum2(3),dum3(3)

      CALL CHRLEN(FN,LENCH)

 7001 Format(1x,'TITLE="Cavitating Propeller"')
 7002 Format(1x,'VARIABLES="x","y","z"')
 7003 Format(1x,'ZONE T="Propeller", N=',I5,', E=',I5,
     *     ' F=FEPOINT, ET=QUADRILATERAL')
      DELK=TWOPI/NBLADE

      WRITE(700,7003) (MR+2)*(NC+1)*NBLADE,(MR+1)*NC*NBLADE

C-----------------------------------------------------------------------
C       The following line is added so the vectors can be properly shown
C       on the key blade.                                       JY093098
C-----------------------------------------------------------------------
Ct      WRITE(62,7003) MR*(NC+1),(MR-1)*NC
C-----------------------------------------------------------------------

      DO 220 K=1, NBLADE

C...........Plot the blades.............................................

C-------The direction of rotation is wrong (JY011000)-------------------
         T=-DELK*FLOAT(K-1)

         DO 140 M=1, MR+1
            DO 130 N=1,NC+1
               XX0 = XB(N,M)
               YY0 = YB(N,M)
               ZZ0 = ZB(N,M)
               IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
               XB1(N,M) = XX0
               YB1(N,M)=YY0*COS(T)-ZB(N,M)*SIN(T)
               ZB1(N,M)=YY0*SIN(T)+ZB(N,M)*COS(T)
               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XB1(N,M),YB1(N,M))
 130        CONTINUE
 140     CONTINUE

C............Generate CONMRX(3,NC+1,MR) and then VECMRX(3,NC+1,MR)......
C............Please note that CONMRX & VECMRX contain >>entire<< blade..
C............So at N=1, we are at the T.E. on pressure side of blade....
C............cm.........................................................

         IF(K.EQ.1) THEN
            DO 145 M=1, MR
               DO 143 N=1, NC+1
                  CONMRX(1,N,M)=
     *                 (XB1(N,M)+XB1(N,M+1))/2
                  CONMRX(2,N,M)=
     *                 (YB1(N,M)+YB1(N,M+1))/2
                  CONMRX(3,N,M)=
     *                 (ZB1(N,M)+ZB1(N,M+1))/2
 143           CONTINUE
 145        CONTINUE

            NTOTAL=0
         ENDIF

C------------First, write out the hub end coordinates of the blade------
C------------CM 040297--------------------------------------------------
         M=1
         DO 150 N=1,NC+1

            WRITE(700,7004) XB1(N,M),YB1(N,M),ZB1(N,M)

 7004       FORMAT(1x,F12.9,2x,F12.9,2x,F12.9)
            IF(K.EQ.1)THEN
               NTOTAL=NTOTAL+1
            ENDIF
 150     CONTINUE

C------------Next, write out the mid-point coordinates of the blade-----
C------------CM 040297--------------------------------------------------

         DO 155 M=1, MR
            DO 153 N=1, NC+1
               WRITE(700,7004) (XB1(N,M)+XB1(N,M+1))/2,
     *              (YB1(N,M)+YB1(N,M+1))/2,
     *              (ZB1(N,M)+ZB1(N,M+1))/2

               IF(K.EQ.1)THEN
                  NTOTAL=NTOTAL+1
               ENDIF
 153        CONTINUE
 155     CONTINUE

C------------Finally, write out the tip end coords CM 040297------------
         M=MR+1
         DO 158 N=1, NC+1
            WRITE(700,7004) XB1(N,M), YB1(N,M), ZB1(N,M)
            IF(K.EQ.1)THEN
               NTOTAL=NTOTAL+1
            ENDIF
 158     CONTINUE


C.............Generate the normal vector matrix.........................
C -- Changed by HSLEE(980324)
C -- Normal vector was calculated at the wrong coordinate
C -- Normal vector should be defined at the mid point of panel side
C -- Include dummy variable DUM to calculate the panel mid point

         IF (K.EQ.1) THEN

            DO 174 M=1, MR
               DO 173 N=1, NC+1

                if(n .eq. 1) then
                 dum1(1) = half*(xb1(2,m+1)+xb1(2,m)
     %                    -xb1(1,m+1)-xb1(1,m) )
                 dum1(2) = half*(yb1(2,m+1)+yb1(2,m)
     %                    -yb1(1,m+1)-yb1(1,m) )
                 dum1(3) = half*(zb1(2,m+1)+zb1(2,m)
     %                    -zb1(1,m+1)-zb1(1,m) )
                elseif(n .eq. nc+1) then
                     dum1(1) = half*(xb1(nc+1,m) +xb1(nc+1,m+1)
     %                    -xb1(nc,m) - xb1(nc,m+1))
                     dum1(2) = half*(yb1(nc+1,m) +yb1(nc+1,m+1)
     %                    -yb1(nc,m) - yb1(nc,m+1))
                     dum1(3) = half*(zb1(nc+1,m) +zb1(nc+1,m+1)
     %                    -zb1(nc,m) - zb1(nc,m+1))
              else
                 nn1 = n+1
                 nn2 = n-1
                     dum1(1) = half*(xb1(nn1,m+1)+xb1(nn1,m)
     %                    -xb1(nn2,m+1)-xb1(nn2,m) )
                     dum1(2) = half*(yb1(nn1,m+1)+yb1(nn1,m)
     %                    -yb1(nn2,m+1)-yb1(nn2,m) )
                     dum1(3) = half*(zb1(nn1,m+1)+zb1(nn1,m)
     %                    -zb1(nn2,m+1)-zb1(nn2,m) )
                endif


              dum2(1) = xb1(n,m+1) - xb1(n,m)
              dum2(2) = yb1(n,m+1) - yb1(n,m)
              dum2(3) = zb1(n,m+1) - zb1(n,m)

              call expro(dum1,dum2,dum3)

              do kk = 1 , 3
                     vecmrx(kk,n,m) = dum3(kk)
              enddo

 173           Continue
 174        Continue


C.........Now lets make our matrix, VECMRX a unit vector matrix!........


            Do 178 M=1, MR
               Do 177 N=1, (NC+1)
                  VECLEN=SQRT((VECMRX(1,N,M)**2)+
     *                 (VECMRX(2,N,M)**2)+
     *                 (VECMRX(3,N,M)**2))
                  IF(VECLEN.LE.0) THEN
                     VECLEN=1.0
                  ENDIF
                  VECMRX(1,N,M)=
     *                 (VECMRX(1,N,M)/VECLEN)
                  VECMRX(2,N,M)=
     *                 (VECMRX(2,N,M)/VECLEN)
                  VECMRX(3,N,M)=
     *                 (VECMRX(3,N,M)/VECLEN)

Ct                  WRITE(62,7020) CONMRX(1,N,M),
Ct     *                 CONMRX(2,N,M),
Ct     *                 CONMRX(3,N,M),
Ct     *                 VECMRX(1,N,M),
Ct     *                 VECMRX(2,N,M),
Ct     *                 VECMRX(3,N,M)
Ct 7020             Format(1x,F12.9,2x,F12.9,2x,F12.9,2x,
Ct     *                 F12.9,2x,F12.9,2x,F12.9)

 177           Continue
 178        Continue

         END IF

 220  Continue

C-----------------------------------------------------------------------
C       The following line is added so the vectors can be properly shown
C       on the key blade.                                       JY093098
C-----------------------------------------------------------------------
Ct      DO 222 M=1, MR-1
Ct         DO 221 N=1, NC
Ct            WRITE(62, 7005) ((M-1)*(NC+1)+N),((M-1)*(NC+1)+N+1),
Ct     *           (M*(NC+1)+N+1),(M*(NC+1)+N)
Ct 221     CONTINUE
Ct 222  CONTINUE
C-----------------------------------------------------------------------

      DO 228 K=1, NBLADE
         DO 226 M=1, MR+1
            DO 224 N=1, NC
               WRITE(700, 7005) ((M-1)*(NC+1)+N)+(K-1)*NTOTAL,
     *              ((M-1)*(NC+1)+N+1)+(K-1)*NTOTAL,
     *              (M*(NC+1)+N+1)+(K-1)*NTOTAL,
     *              (M*(NC+1)+N)+(K-1)*NTOTAL
 7005          FORMAT(1x,I5,2x,I5,2x,I5,2x,I5)
 224        CONTINUE
 226     CONTINUE
 228  CONTINUE

      open(610,file='blade.dat',status='unknown')
      write(610,*) nc+1,mr+1
      do m = 1, mr+1
         do n = 1, nc+1
            write(610,*) xb(n,m),yb(n,m),-zb(n,m)
         enddo
      enddo
      close(610)
C
C......Plot the hub.....................................................
C
      IF(IHUB.NE.0) THEN
         OPEN(610,FILE='hub.plt',STATUS='UNKNOWN')
         WRITE(610,*) ' TITLE = "Plot Hub Geometry" '
         WRITE(610,7002)
         Write(610,7030) ((NHBX+1)*(MHBT+1))*NBLADE,
     *        ((NHBX)*(MHBT))*NBLADE
 7030    Format(1x,'ZONE T="Main part of Hub", N=',I6,', E=',I6,
     *        ', F=FEPOINT, ET=QUADRILATERAL')
 7040    FORMAT(1x, F12.9,5x,F12.9,5x,F12.9)

         NTOTAL=0
         DO 232 K=1,NBLADE

            T=-DELK*FLOAT(K-1)

            DO 234 M=1,MHBT+1
               DO 236 N=1,NHBX+1
                  XX0 = XH(N,M)
                  YY0 = YH(N,M)
                  ZZ0 = ZH(N,M)
                  IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
                  XH1(N,M) = XX0
                  YH1(N,M)=YY0*COS(T)-ZH(N,M)*SIN(T)
                  ZH1(N,M)=YY0*SIN(T)+ZH(N,M)*COS(T)
                IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XH1(N,M),YH1(N,M))
 236           CONTINUE
 234        CONTINUE

            DO 240 N=1,NHBX+1
               DO 230 M=1,MHBT+1
                  WRITE(610,7040) XH1(N,M),YH1(N,M),ZH1(N,M)
                  IF(K.EQ.1) NTOTAL=NTOTAL+1
 230           CONTINUE
 240        CONTINUE
 232     CONTINUE

         DO 238 K=1,NBLADE
            DO 245 N=1, NHBX
               DO 243 M=1,MHBT
                  Write(610,7005) ((N-1)*(MHBT+1)+M)+(K-1)*NTOTAL,
     *                 ((N-1)*(MHBT+1)+M+1)+(K-1)*NTOTAL,
     *                 (N*(MHBT+1)+M+1)+(K-1)*NTOTAL,
     *                 (N*(MHBT+1)+M)+(K-1)*NTOTAL
 243           CONTINUE
 245        CONTINUE
 238     CONTINUE
         CLOSE(610)
      END IF


      IF(IDUCT .NE. 0) THEN
         open(998,file='ductsct.plt',status='unknown')

         write(998,*) 'variables =x,y,z'
         WRITE (998,*)
     %        'ZONE T="FACE", I=', nducthp, ', J = ',
     %        mductp , ', K = 1 '

         DO M = 1 , MDUCTP
            do n = 1, nducthp
               write(998,*) xd(n,m),yd(n,m),zd(n,m)
            enddo
         enddo

         WRITE (998,*)
     %        'ZONE T="BACK", I=', nducthp, ', J = ',
     %        mductp , ', K = 1 '
         DO M = 1 , MDUCTP
            do n = 1, nducthp
               write(998,*) xd(nducth+n,m),yd(nducth+n,m),zd(nducth+n,m)
            enddo
         enddo

         WRITE (998,*)
     %        'ZONE T="HUb", I=', nhbx+1, ', J = ',
     %        mhbt+1 , ', K = 1 '
         DO M=1,MHBT+1
            DO N=1,NHBX+1
               write(998,*) xh(n,m),yh(n,m),zh(n,m)
            enddo
         enddo

         close(998)
      ENDIF

      write(*,*) ' End of writing geometry plot '

      RETURN
C))))))))))))))))))))) End of subroutine PROPLT ((((((((((((((((((((((((
      END












