       SUBROUTINE rudimv_bb (iv1,iv2,iv3,iv4,xv13)
c-------Subroutine added on 04/2012 (XM YU)------------
csingh-------iv1 = J (panel index) -------------------------------------
csingh-------iv2 = I (panel at which the influence is considered)-------
csingh-------iv3 = KK (Blade index)-------------------------------------
csingh-------iv4 = IMR (determines if Rpan needs to be called or Hypot)-
c            xv13 = source influence from images
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC' 

        iv3=1
        iv4=1
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

           xv13 =  xv13 + FS

100     CONTINUE

        RETURN
        END

       SUBROUTINE rudimv_wb (iv1,iv2,iv3,iv4,xv13)
c-------Subroutine added on 04/2012 (XM YU)------------
csingh-------iv1 = J (panel index) -------------------------------------
csingh-------iv2 = I (panel at which the influence is considered)-------
csingh-------iv3 = KK (Blade index)-------------------------------------
csingh-------iv4 = IMR (determines if Rpan needs to be called or Hypot)-
c            xv13 = source influence from images
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC' 
        INCLUDE 'PUFCAVC.INC' 

        iv3=1
        iv4=1
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
                  xv17 = xv9 - xctw(iv2,K)
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
                  XLOC = XLOC + (xctw(iv2,K)-XCT(iv1,K))*
     *                   DIR(iv1,1,K)
                  YLOC = YLOC + (xctw(iv2,K)-XCT(iv1,K))*
     *                   DIR(iv1,2,K)
                  ZLOC = ZLOC + (xctw(iv2,K)-XCT(iv1,K))*
     *                   DIR(iv1,3,K)
               ENDIF
10          CONTINUE

            CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                 FSY,FDX,FDY,FDZ,0,iv4)

            xv13 =  xv13 + FS
csingh--------Considering bottom wall ----- ----------------------------
            XLOC=ZERO
            YLOC=ZERO 
            ZLOC=ZERO
            Do 30 iv20=1,3
               IF (iv20.EQ.2) THEN
                  xv19 = xctw(iv2,iv20) - xv14 
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
                  XLOC = XLOC + (xctw(iv2,iv20)-XCT(iv1,iv20))*
     *                   DIR(iv1,1,iv20)
                  YLOC = YLOC + (xctw(iv2,iv20)-XCT(iv1,iv20))*
     *                   DIR(iv1,2,iv20)
                  ZLOC = ZLOC + (xctw(iv2,iv20)-XCT(iv1,iv20))*
     *                   DIR(iv1,3,iv20)
               ENDIF
30         CONTINUE

           CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                FSY,FDX,FDY,FDZ,0,iv4)

           xv13 =  xv13 + FS

100     CONTINUE

        RETURN
        END

       SUBROUTINE rudimv_bw (iv1,iv2,iv3,iv4,xv13)
c-------Subroutine added on 04/2012 (XM YU)------------
csingh-------iv1 = J (panel index) -------------------------------------
csingh-------iv2 = I (panel at which the influence is considered)-------
csingh-------iv3 = KK (Blade index)-------------------------------------
csingh-------iv4 = IMR (determines if Rpan needs to be called or Hypot)-
c            xv13 = source influence from images
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC' 
        INCLUDE 'PUFCAVC.INC' 

        iv3=1
        iv4=1
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
                  xv17 = xv9 - xctp(iv2,K,iv3)
                  IF (mod(iv16,2).GT.0) THEN 
                     xv10 = xv9 + (iv16-1)*xv18 + xv17 
                  ELSE
                     xv10 = xv9 + iv16*xv18 - xv17 
                  ENDIF
csingh-------xv10 = y coordinate of the image --------------------------

                  XLOC = XLOC + (xv10-xctw(iv1,K))*DIRw(iv1,1,K)
                  YLOC = YLOC + (xv10-XCTw(iv1,K))*DIRw(iv1,2,K)
                  ZLOC = ZLOC + (xv10-XCTw(iv1,K))*DIRw(iv1,3,K)
               ELSE
                  XLOC = XLOC + (xctp(iv2,K,iv3)-XCTw(iv1,K))*
     *                   DIRw(iv1,1,K)
                  YLOC = YLOC + (xctp(iv2,K,iv3)-XCTw(iv1,K))*
     *                   DIRw(iv1,2,K)
                  ZLOC = ZLOC + (xctp(iv2,K,iv3)-XCTw(iv1,K))*
     *                   DIRw(iv1,3,K)
               ENDIF
10          CONTINUE

            CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                 FSY,FDX,FDY,FDZ,0,iv4)

            xv13 =  xv13 + FS
csingh--------Considering bottom wall ----- ----------------------------
            XLOC=ZERO
            YLOC=ZERO 
            ZLOC=ZERO
            Do 30 iv20=1,3
               IF (iv20.EQ.2) THEN
                  xv19 = xctp(iv2,iv20,iv3) - xv14 
                  IF (mod(iv16,2).GT.0) THEN 
                     xv10 = xv14 - (iv16-1)*xv18 - xv19 
                  ELSE
                     xv10 = xv14 - iv16*xv18 + xv19 
                  ENDIF
csingh-------xv10 = y coordinate of the image --------------------------

                  XLOC = XLOC + (xv10-XCTw(iv1,iv20))*DIRw(iv1,1,iv20)
                  YLOC = YLOC + (xv10-XCTw(iv1,iv20))*DIRw(iv1,2,iv20)
                  ZLOC = ZLOC + (xv10-XCTw(iv1,iv20))*DIRw(iv1,3,iv20)
               ELSE
                  XLOC = XLOC + (xctp(iv2,iv20,iv3)-XCTw(iv1,iv20))*
     *                   DIRw(iv1,1,iv20)
                  YLOC = YLOC + (xctp(iv2,iv20,iv3)-XCTw(iv1,iv20))*
     *                   DIRw(iv1,2,iv20)
                  ZLOC = ZLOC + (xctp(iv2,iv20,iv3)-XCTw(iv1,iv20))*
     *                   DIRw(iv1,3,iv20)
               ENDIF
30         CONTINUE

           CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                FSY,FDX,FDY,FDZ,0,iv4)

           xv13 =  xv13 + FS

100     CONTINUE

        RETURN
        END

       SUBROUTINE rudimv_ww (iv1,iv2,iv3,iv4,xv13)
c-------Subroutine added on 04/2012 (XM YU)------------
csingh-------iv1 = J (panel index) -------------------------------------
csingh-------iv2 = I (panel at which the influence is considered)-------
csingh-------iv3 = KK (Blade index)-------------------------------------
csingh-------iv4 = IMR (determines if Rpan needs to be called or Hypot)-
c            xv13 = source influence from images
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC' 
        INCLUDE 'PUFCAVC.INC' 

        iv3=1
        iv4=1
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
                  xv17 = xv9 - xctw(iv2,K)
                  IF (mod(iv16,2).GT.0) THEN 
                     xv10 = xv9 + (iv16-1)*xv18 + xv17 
                  ELSE
                     xv10 = xv9 + iv16*xv18 - xv17 
                  ENDIF
csingh-------xv10 = y coordinate of the image --------------------------

                  XLOC = XLOC + (xv10-xctw(iv1,K))*DIRw(iv1,1,K)
                  YLOC = YLOC + (xv10-XCTw(iv1,K))*DIRw(iv1,2,K)
                  ZLOC = ZLOC + (xv10-XCTw(iv1,K))*DIRw(iv1,3,K)
               ELSE
                  XLOC = XLOC + (xctw(iv2,K)-XCTw(iv1,K))*
     *                   DIRw(iv1,1,K)
                  YLOC = YLOC + (xctw(iv2,K)-XCTw(iv1,K))*
     *                   DIRw(iv1,2,K)
                  ZLOC = ZLOC + (xctw(iv2,K)-XCTw(iv1,K))*
     *                   DIRw(iv1,3,K)
               ENDIF
10          CONTINUE

            CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                 FSY,FDX,FDY,FDZ,0,iv4)

            xv13 =  xv13 + FS
csingh--------Considering bottom wall ----- ----------------------------
            XLOC=ZERO
            YLOC=ZERO 
            ZLOC=ZERO
            Do 30 iv20=1,3
               IF (iv20.EQ.2) THEN
                  xv19 = xctw(iv2,iv20) - xv14 
                  IF (mod(iv16,2).GT.0) THEN 
                     xv10 = xv14 - (iv16-1)*xv18 - xv19 
                  ELSE
                     xv10 = xv14 - iv16*xv18 + xv19 
                  ENDIF
csingh-------xv10 = y coordinate of the image --------------------------

                  XLOC = XLOC + (xv10-XCTw(iv1,iv20))*DIRw(iv1,1,iv20)
                  YLOC = YLOC + (xv10-XCTw(iv1,iv20))*DIRw(iv1,2,iv20)
                  ZLOC = ZLOC + (xv10-XCTw(iv1,iv20))*DIRw(iv1,3,iv20)
               ELSE
                  XLOC = XLOC + (xctw(iv2,iv20)-XCTw(iv1,iv20))*
     *                   DIRw(iv1,1,iv20)
                  YLOC = YLOC + (xctw(iv2,iv20)-XCTw(iv1,iv20))*
     *                   DIRw(iv1,2,iv20)
                  ZLOC = ZLOC + (xctw(iv2,iv20)-XCTw(iv1,iv20))*
     *                   DIRw(iv1,3,iv20)
               ENDIF
30         CONTINUE

           CALL RPAN (XLOC,YLOC,ZLOC,CHRLEPS(iv1),FS,FD,FSX,
     *                FSY,FDX,FDY,FDZ,0,iv4)

           xv13 =  xv13 + FS

100     CONTINUE

        RETURN
        END
