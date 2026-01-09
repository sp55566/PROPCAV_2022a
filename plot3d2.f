C
C      Subroutine to plot tip vortex cavity convergence check
C      10/13/99  HSLEE
C---------------------------------------------------------------
      subroutine plot3d2
C --------------------------------------------------------------
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC'
        INCLUDE 'PUFCAVC.INC'

      if(nrev .eq. 1) then
         iunt = 910
      elseif(nrev .eq. 2) then
         iunt = 911
C/s S.N.KIM | increases the number of files plotting unsteady wake geo.
C           | at each revolution.
      elseif(nrev .eq. 3) then
         iunt = 913
      elseif(nrev .eq. 4) then
         iunt = 914
      elseif(nrev .eq. 5) then
         iunt = 915
      elseif(nrev .eq. ntrev) then
         iunt = 916
C/e S.N.KIM | Aug. 2018.
      endif

      ddelta = -tstep

      ittt = itstep 

Ct      nbt1 = nblade
      nbt1 = 1

c      ntt = nblade*mr*nc + nbt1*mr*nwpanel+
c     %            (nthx+ncvx)*mcvt
      ntt = nblade*mr*nc + nbt1*mr*nwpanel
 
c      nnode = nblade*mrp*ncp + (ncvxp+nthxp)*mcvtp 
c     %            + nbt1*(nwpanel+1)*mrp
      nnode = nblade*mrp*ncp + nbt1*(nwpanel+1)*mrp

      if(ihub .ne. 0) then
         ntt = ntt + nblade*nhbx*mhbt
         nnode = nnode + nblade*(nhbx+1)*(mhbt+1)
      endif

      write(iunt,100) ittt,NNODE, NTT

      i1 = 0
        nstb = i1
      do m = 1, mrp
        do n = 1 , ncp
           i1 = i1 + 1
           yyd = yb(n,m)*cos(ddelta)-zb(n,m)*sin(ddelta)
           zzd = yb(n,m)*sin(ddelta)+zb(n,m)*cos(ddelta)
           write(iunt,300) xb(n,m),yyd,zzd
        enddo
      enddo

       do ik = 1 , nblade-1

        ddeltk = -delk * ik + ddelta

        do m = 1 , mrp
          do n = 1 , ncp
             i1 = i1 + 1
             yyd = yb(n,m)*cos(ddeltk)-zb(n,m)*sin(ddeltk)
             zzd = yb(n,m)*sin(ddeltk)+zb(n,m)*cos(ddeltk)
             write(iunt,300) xb(n,m),yyd,zzd
          enddo
        enddo
      enddo

      if(ihub .ne. 0) then
         nsthb = i1
         do m = 1, mhbt+1
            do n = 1 , nhbx+1
             i1 = i1 + 1
             yyd = yh(n,m)*cos(ddelta)-zh(n,m)*sin(ddelta)
             zzd = yh(n,m)*sin(ddelta)+zh(n,m)*cos(ddelta)
             write(iunt,300) xh(n,m),yyd,zzd
            enddo
         enddo

         do ik = 1 , nblade-1
            
            ddeltk = -delk * ik + ddelta
            
            do m = 1 , mhbt+1
             do n = 1 , nhbx+1
                i1 = i1 + 1
                yyd = yh(n,m)*cos(ddeltk)-zh(n,m)*sin(ddeltk)
                zzd = yh(n,m)*sin(ddeltk)+zh(n,m)*cos(ddeltk)
                write(iunt,300) xh(n,m),yyd,zzd
             enddo
            enddo
         enddo
      endif


        nstw = i1
      do m = 1 , mrp
        do n = 1 , nwpanel+1 
           i1 = i1 + 1
           yyd = yw(n,m)*cos(ddelta)-zw(n,m)*sin(ddelta)
           zzd = yw(n,m)*sin(ddelta)+zw(n,m)*cos(ddelta)
           write(iunt,300) xw(n,m),yyd,zzd
        enddo
      enddo

c        do ik = 2 , nblade
c          do m = 1 , mrp
c            do n = 1 , nwpanel+1 
c             i1 = i1 + 1
c             yyd = yww(n,m,ik)*cos(ddelta)-zww(n,m,ik)*sin(ddelta)
c             zzd = yww(n,m,ik)*sin(ddelta)+zww(n,m,ik)*cos(ddelta)
c              write(iunt,300) xww(n,m,ik),yyd,zzd
c            enddo
c          enddo
c        enddo

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
c          nsth = i1
c        do m = 1 , mcvtp
c          do n = 1 , nthxp 
c           i1 = i1 + 1
c           yyd = ych(n,m)*cos(ddelta)-zch(n,m)*sin(ddelta)
c           zzd = ych(n,m)*sin(ddelta)+zch(n,m)*cos(ddelta)
c            write(iunt,300) xch(n,m),yyd,zzd
c          enddo
c        enddo
c
c          nstc = i1
c        do m = 1 , mcvtp
c          do n = 1 , ncvxp 
c           i1 = i1 + 1
c           yyd = yvc(n,m)*cos(ddelta)-zvc(n,m)*sin(ddelta)
c           zzd = yvc(n,m)*sin(ddelta)+zvc(n,m)*cos(ddelta)
c            write(iunt,300) xvc(n,m),yyd,zzd
c          enddo
c        enddo
C/e S.N.KIM | Aug. 2018.

C ---- Connectivity

       DO IK = 1 , NBLADE
        DO M = 1 , MR
          DO N = 1 , NC
             J1 = CONECT(N,M,1,IK,NSTB)
             J2 = CONECT(N+1,M,1,IK,NSTB)
             J3 = CONECT(N+1,M+1,1,IK,NSTB)
             J4 = CONECT(N,M+1,1,IK,NSTB)
             WRITE(iunt,200) J1,J2,J3,J4
          ENDDO
        ENDDO
       ENDDO

      IF(IHUB .NE. 0) THEN
         DO IK = 1 , NBLADE
            DO M = 1 , MHBT
             DO N = 1 , NHBX
                J1 = CONECT(N,M,2,IK,NSTHB)
                J2 = CONECT(N+1,M,2,IK,NSTHB)
                J3 = CONECT(N+1,M+1,2,IK,NSTHB)
                J4 = CONECT(N,M+1,2,IK,NSTHB)
                WRITE(iunt,200) J1,J2,J3,J4
             ENDDO
            ENDDO
         ENDDO
      ENDIF

c       DO IK = 1 , NBLADE
      IK = 1
        DO M = 1 , MR
          DO N = 1 , NWPANEL 
             J1 = CONECT(N,M,3,IK,NSTW)
             J2 = CONECT(N+1,M,3,IK,NSTW)
             J3 = CONECT(N+1,M+1,3,IK,NSTW)
             J4 = CONECT(N,M+1,3,IK,NSTW)
             WRITE(iunt,200) J1,J2,J3,J4
          ENDDO
        ENDDO
c       ENDDO

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
c        DO M = 1 , MCVT
c          DO N = 1 , NTHX 
c             J1 = CONECT(N,M,4,0,NSTH)
c             J2 = CONECT(N+1,M,4,0,NSTH)
c             J3 = CONECT(N+1,M+1,4,0,NSTH)
c             J4 = CONECT(N,M+1,4,0,NSTH)
c             WRITE(iunt,200) J1,J2,J3,J4
c          ENDDO
c        ENDDO
c
c        DO M = 1 , MCVT
c          DO N = 1 , NCVX 
c             J1 = CONECT(N,M,5,0,NSTC)
c             J2 = CONECT(N+1,M,5,0,NSTC)
c             J3 = CONECT(N+1,M+1,5,0,NSTC)
c             J4 = CONECT(N,M+1,5,0,NSTC)
c             WRITE(iunt,200) J1,J2,J3,J4
c          ENDDO
c        ENDDO
C/e S.N.KIM | Aug. 2018.

100     FORMAT('ZONE T="ANG=',I4,'", F=FEPOINT,',
     %         'ET=QUADRILATERAL,N=',I5,',E=',I5)
200      FORMAT(4(1x,I6))
300     FORMAT(3(1X,F12.6))

        RETURN
        END



C --------------------------------------
        FUNCTION CONECT(N,M,IP,KK0,NST)
C ---------------------------------------
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC'
        INCLUDE 'PUFCAVC.INC'
C  KK0 = kk0-th blade or wake
C
C  IP = 1 : Blade
C       2 : HUB
C       3 : Wake
C       4 : tip Bulb
C       5 : Tip Vortex
C ---------------------------------

        IF(IP .EQ. 1) Then
          CONECT = NST + (KK0-1)*(MRP*NCP) 
     %                  + (M - 1)*NCP+N
        ELSEIF(IP .EQ. 2) Then
          CONECT = NST + (KK0-1)*((MHBT+1)*(NHBX+1)) 
     %                  + (M - 1)*(NHBX+1)+N
        ELSEIF(IP .EQ. 3) THEN
          CONECT = NST + (KK0-1)*(NWPANEL+1)*MRP
     %                    +(M-1)*(NWPANEL+1)+N 
        ELSEIF(IP .EQ. 4) THEN
          CONECT = NST + (M-1)*NTHXP + N 
        ELSEIF(IP .EQ. 5) THEN
          CONECT = NST + (M-1)*NCVXP + N 
        ENDIF

        RETURN
        END

