  	SUBROUTINE RUDIM (iv1,iv2,iv3,iv4,xv5,xv6,xv7,xv8,xv12,xv13)
csingh-------Subroutine added on 08/26/2008 (Sowmitra Singh)------------
csingh-------iv1 = J (panel index) -------------------------------------
csingh-------iv2 = I (panel at which the influence is considered)-------
csingh-------iv3 = KK (Blade index)-------------------------------------
csingh-------iv4 = IMR (determines if Rpan needs to be called or Hypot)-
csingh-------xv5 = XM1 (stores coordinates of panel N,M)----------------
csingh-------xv6 = XM2 (Stores coordinates of panel N,M+1)--------------
csingh-------xv7 = XM3 (Stores coordinates of panel N+1,M+1)------------
csingh-------xv8 = XM4 (Stores coordinates of panel N+1,M)--------------
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC' 
        DIMENSION xv5(3),xv6(3),xv7(3),xv8(3),XMC(3)

        xv12 = ZERO
        xv13 = ZERO
        xv9 = 1
        xv14 = YB(1,1)
csingh------xv14 (lower wall), xv9 (upper wall) ------------------------
        xv18 = xv9 - xv14 
        iv15 = wingimag
csingh-----sets of images that are considered on each wall-------------- 
        Do 100 iv16=1,iv15
csingh--------Considering top wall ----- -------------------------------
            XLOC=ZERO
            YLOC=ZERO 
            ZLOC=ZERO
            Do 10 K=1,3
               IF (K.EQ.2) THEN
                  xv17 = xv9 - XCTP(iv2,K,iv3)
                  IF (mod(iv16,2).GT.0) THEN 
                     xv10 = xv9 + (iv16-1)*xv18 + xv17 
                  ELSE
                     xv10 = xv9 + iv16*xv18 - xv17 
                  ENDIF
csingh-------xv10 = y coordinate of the image --------------------------

                  XLOC = XLOC + (xv10-XCT(iv1,K))*DIR(iv1,1,K)
                  YLOC = YLOC + (xv10-XCT(iv1,K))*DIR(iv1,2,K)
                  ZLOC = ZLOC + (xv10-XCT(iv1,K))*DIR(iv1,3,K)
               ELSE
                  XLOC = XLOC + (XCTP(iv2,K,iv3)-XCT(iv1,K))*
     *                   DIR(iv1,1,K)
                  YLOC = YLOC + (XCTP(iv2,K,iv3)-XCT(iv1,K))*
     *                   DIR(iv1,2,K)
                  ZLOC = ZLOC + (XCTP(iv2,K,iv3)-XCT(iv1,K))*
     *                   DIR(iv1,3,K)
               ENDIF
10          CONTINUE

            CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                 FSY,FDX,FDY,FDZ,0,iv4)

            IF (iv4.EQ.2) THEN
csingh-----HYPOT being called ------------------------------------------
               DO 20 iv11=1,3
                  IF (iv11.EQ.2) THEN
                     xv17 = xv9 - XCTP(iv2,iv11,iv3)
                     IF (mod(iv16,2).GT.0) THEN 
                        XMC(iv11) = xv9 + (iv16-1)*xv18 + xv17 
                     ELSE
                        XMC(iv11) = xv9 + iv16*xv18 - xv17 
                     ENDIF
csingh--------y co-ordinate of the image--------------------------
                  ELSE
                     XMC(iv11)=XCTP(iv2,iv11,iv3)
                  ENDIF
20             CONTINUE 
               CALL HYPOT(xv5,xv6,xv7,xv8,XMC,FD,FS)
            ENDIF

            xv12 =  xv12 + FD
            xv13 =  xv13 + FS
csingh--------Considering bottom wall ----- ----------------------------
            XLOC=ZERO
            YLOC=ZERO 
            ZLOC=ZERO
            Do 30 iv20=1,3
               IF (iv20.EQ.2) THEN
                  xv19 = XCTP(iv2,iv20,iv3) - xv14 
                  IF (mod(iv16,2).GT.0) THEN 
                     xv10 = xv14 - (iv16-1)*xv18 - xv19 
                  ELSE
                     xv10 = xv14 - iv16*xv18 + xv19 
                  ENDIF
csingh-------xv10 = y coordinate of the image --------------------------

                  XLOC = XLOC + (xv10-XCT(iv1,iv20))*DIR(iv1,1,iv20)
                  YLOC = YLOC + (xv10-XCT(iv1,iv20))*DIR(iv1,2,iv20)
                  ZLOC = ZLOC + (xv10-XCT(iv1,iv20))*DIR(iv1,3,iv20)
               ELSE
                  XLOC = XLOC + (XCTP(iv2,iv20,iv3)-XCT(iv1,iv20))*
     *                   DIR(iv1,1,iv20)
                  YLOC = YLOC + (XCTP(iv2,iv20,iv3)-XCT(iv1,iv20))*
     *                   DIR(iv1,2,iv20)
                  ZLOC = ZLOC + (XCTP(iv2,iv20,iv3)-XCT(iv1,iv20))*
     *                   DIR(iv1,3,iv20)
               ENDIF
30         CONTINUE

           CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                FSY,FDX,FDY,FDZ,0,iv4)

           IF (iv4.EQ.2) THEN
csingh-----HYPOT being called ------------------------------------------
              DO 40 iv21=1,3
                 IF (iv21.EQ.2) THEN
                    xv19 = XCTP(iv2,iv21,iv3) - xv14 
                    IF (mod(iv16,2).GT.0) THEN 
                       XMC(iv21) = xv14 - (iv16-1)*xv18 - xv19
                    ELSE
                       XMC(iv21) = xv14 - iv16*xv18 + xv19
                    ENDIF
csingh--------y co-ordinate of the image-------------------------------
                 ELSE
                    XMC(iv21)=XCTP(iv2,iv21,iv3)
                 ENDIF
40            CONTINUE 

              CALL HYPOT(xv5,xv6,xv7,xv8,XMC,FD,FS)
           ENDIF

csingh------------------------------------------------------------------

           xv12 =  xv12 + FD
           xv13 =  xv13 + FS

100     CONTINUE

        RETURN
        END
