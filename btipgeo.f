c----------------------------------------------------------------------
        SUBROUTINE BTIPGEO 
c----------------------------------------------------------------------
c
C   Original version : April 28 1999
c   latest revision  : April 28 1999
c
c   purpose          : Generate tip cavity shape at the blade tip 
C                      Assume shape to be circle 
c
c   usage            : call BTIPGEO 
c
c   arguments :
c     Input   
c     Output
c
c   required  routines
c
c----------------------------------------------------------------------
c
        INCLUDE 'PUFCAV.INC'
        dimension slop(nhmxp,3), slop2(nhmxp,3),slop3(nhmxp,3)
        dimension unn(3),un1(3),un2(3),un3(3),slopd(nhmxp,3)
        dimension xchc(nhmxp),ychc(nhmxp),zchc(nhmxp)

      do n = 1 , nhp
        xch(n,1) = xb(nh+n,mrp) 
        ych(n,1) = yb(nh+n,mrp)
        zch(n,1) = zb(nh+n,mrp)
      enddo

      do m = 1 , mcvtp
          xch(1,m) = xb(nhp,mrp)      
          ych(1,m) = yb(nhp,mrp)      
          zch(1,m) = zb(nhp,mrp)      
      enddo

      xchc(1) = xch(1,1)
      ychc(1) = ych(1,1)
      zchc(1) = zch(1,1)

      dtt = twopi / real(mcvt)

        xlen = 0.0
        do n = 1 , nh
          dx1 = xb(nhp+n,mrp) - xb(nhp+n-1,mrp)
          dy1 = yb(nhp+n,mrp) - yb(nhp+n-1,mrp)
          dz1 = zb(nhp+n,mrp) - zb(nhp+n-1,mrp)
          xlen = xlen + sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
        enddo

        do n = 1 , nh
          un1(1) = xb(nhp+n,mrp) - xb(nhp+n-1,mrp)
          un1(2) = yb(nhp+n,mrp) - yb(nhp+n-1,mrp)
          un1(3) = zb(nhp+n,mrp) - zb(nhp+n-1,mrp)
          duu = sqrt(un1(1)**2 + un1(2)**2 + un1(3)**2)
          slop(n,1) = un1(1) / duu
          slop(n,2) = un1(2) / duu
          slop(n,3) = un1(3) / duu

C --- Check unit vector along the wake mr-(mr+1) line at n = i

          un2(1) = xb(nhp+n,mrp) - 
     %               half * (xb(nhp+n,mr)+xb(nhp-n,mr))
          un2(2) = yb(nhp+n,mrp) - 
     %               half * (yb(nhp+n,mr)+yb(nhp-n,mr))
          un2(3) = zb(nhp+n,mrp) - 
     %               half * (zb(nhp+n,mr)+zb(nhp-n,mr))

          duu = sqrt(un2(1)**2 + un2(2)**2 + un2(3)**2)
          slop2(n,1) = un2(1) / duu
          slop2(n,2) = un2(2) / duu
          slop2(n,3) = un2(3) / duu

C --- Find unit vector normali to slop & slop2 vector

          call expro(un1,un2,unn)

          duu = sqrt(unn(1)**2 + unn(2)**2 + unn(3)**2)

          unn(1) = unn(1) / duu
          unn(2) = unn(2) / duu
          unn(3) = unn(3) / duu

          slopd(n,1) = unn(1)
          slopd(n,2) = unn(2)
          slopd(n,3) = unn(3)

          call expro(unn,un1,un3)

          duu = sqrt(un3(1)**2 + un3(2)**2 + un3(3)**2)

          slop3(n,1) = un3(1)/duu
          slop3(n,2) = un3(2)/duu
          slop3(n,3) = un3(3)/duu
        enddo

        dxlen = 0.0
        do n = 1 , nh
          dx1 = xb(nhp+n,mrp) - xb(nhp+n-1,mrp)
          dy1 = yb(nhp+n,mrp) - yb(nhp+n-1,mrp)
          dz1 = zb(nhp+n,mrp) - zb(nhp+n-1,mrp)
          dxlen = dxlen + sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
        fact = dxlen / xlen
        rc2 = rc1 * fact
          xchc(n+1) = rc2 * slop3(n,1) + xb(nhp+n,mrp)
          ychc(n+1) = rc2 * slop3(n,2) + yb(nhp+n,mrp)
          zchc(n+1) = rc2 * slop3(n,3) + zb(nhp+n,mrp)
          do m = 1 , mcvt+1
             dtt2 = dtt * real(m-1)
             ss2 = rc2 * (1.-cos(dtt2))
             ss3 = rc2 * sin(dtt2)

             dxx = slop3(n,1)*ss2 + slopd(n,1)*ss3
             dyy = slop3(n,2)*ss2 + slopd(n,2)*ss3
             dzz = slop3(n,3)*ss2 + slopd(n,3)*ss3

             xch(n+1,m) = xb(nhp+n,mrp) + dxx
             ych(n+1,m) = yb(nhp+n,mrp) + dyy
             zch(n+1,m) = zb(nhp+n,mrp) + dzz
          enddo
        enddo

      do m = 1 , mcvt
        xch(nhp,m) = xvc(1,m)
        ych(nhp,m) = yvc(1,m)
        zch(nhp,m) = zvc(1,m)
        enddo
  
      return
      end
 
