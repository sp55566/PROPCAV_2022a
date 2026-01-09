C
C      Subroutine to plot tip vortex cavity convergence check
C      10/13/99  HSLEE
C---------------------------------------------------------------
      SUBROUTINE PLOT3D
C --------------------------------------------------------------
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC'
        INCLUDE 'PUFCAVC.INC'
      
      ddelta = -tstep

      ittt = itstep 

      write(912,100)
      write(912,210) ittt,ncp, mrp

      do m = 1 , mrp
        do n = 1 , ncp
           yyd = yb(n,m)*cos(ddelta)-zb(n,m)*sin(ddelta)
           zzd = yb(n,m)*sin(ddelta)+zb(n,m)*cos(ddelta)
           write(912,300) xb(n,m),yyd,zzd
        enddo
      enddo

      write(912,220) ittt,nwpanel+1, mrp

      do m = 1 , mrp
        do n = 1 , nwpanel+1 
           yyd = yw(n,m)*cos(ddelta)-zw(n,m)*sin(ddelta)
           zzd = yw(n,m)*sin(ddelta)+zw(n,m)*cos(ddelta)
           write(912,300) xw(n,m),yyd,zzd
        enddo
      enddo

        write(912,230) ittt,nthxp, mcvtp 

        do m = 1 , mcvtp
          do n = 1 , nthxp 
           yyd = ych(n,m)*cos(ddelta)-zch(n,m)*sin(ddelta)
           zzd = ych(n,m)*sin(ddelta)+zch(n,m)*cos(ddelta)
            write(912,300) xch(n,m),yyd,zzd
          enddo
        enddo

        write(912,240) ittt,ncvxp, mcvtp 
        do m = 1 , mcvtp
          do n = 1 , ncvxp 
           yyd = yvc(n,m)*cos(ddelta)-zvc(n,m)*sin(ddelta)
           zzd = yvc(n,m)*sin(ddelta)+zvc(n,m)*cos(ddelta)
            write(912,300) xvc(n,m),yyd,zzd
          enddo
        enddo

100     FORMAT(' VARIABLES =X,Y,Z')
210     FORMAT('ZONE T="P,AG=',i4,'",
     %            I=',I4,',J = ',I4,',F=POINT')
220     FORMAT('ZONE T="W, AG=',i4,'",
     %            I=',I4,',J = ',I4,',F=POINT')
230     FORMAT('ZONE T="B, AG=',i4,'",
     %            I =',I4,',J = ',I4,',F=POINT')
240     FORMAT('ZONE T="C, AG=',i4,'",
     %                 I =',I4,',J = ',I4,',F=POINT')

300     FORMAT(3(1X,F12.6))

        RETURN
        END

