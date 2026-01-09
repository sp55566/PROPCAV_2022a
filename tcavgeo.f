C ----------------------------------------
        SUBROUTINE TCAVGEO
C ----------------------------------------
        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC'

        dimension slop(ncavmp,3), slop2(ncavmp,3),slop3(ncavmp,3)
        dimension unn(3),un1(3),un2(3),un3(3),slopd(ncavmp,3)
        dimension cleng(ncavmp),fact(ncavmp)

C --- Check unit vector along the wake i-(i+1) line at m = mrp

c      if(icavt .eq. 0) then
        rc1 = radini 
c      elseif(icavt .ge. 1 .and. icavt .le. 4) then
c        sigmal = sigma / advco**2
c          rc1 = abs(gamma) / twopi / sqrt(sigmal)
c        endif

        dtt = twopi / real(mcvt)

        write(*,*) ' gamma : ', gamma
        WRITE(*,*) ' Tip radius :', rc1

        cleng(1) = 0.0
        do n = 2 , ncvxp
           dx1 = xw(n,mrp) - xw(n-1,mrp)
           dy1 = yw(n,mrp) - yw(n-1,mrp)
           dz1 = zw(n,mrp) - zw(n-1,mrp)
           dxyz = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
           cleng(n) = cleng(n-1) + dxyz
        enddo

C       Hong corrected the initialization of FACT 09/13/07
c      fact = 1.0
        DO N = 1, NCVXP
           FACT(N) = 1.0 
        ENDDO 

        do m = 1 , mcvtp
          xvc(1,m) = xb(1,mrp)
          yvc(1,m) = yb(1,mrp)
          zvc(1,m) = zb(1,mrp)
        enddo

        do n = 1 , ncvxp
          if(n .eq. 1) then
            un1(1) = xw(2,mrp) - xw(1,mrp)
            un1(2) = yw(2,mrp) - yw(1,mrp)
            un1(3) = zw(2,mrp) - zw(1,mrp)
          else
            un1(1) = xw(n,mrp) - xw(n-1,mrp)
            un1(2) = yw(n,mrp) - yw(n-1,mrp)
            un1(3) = zw(n,mrp) - zw(n-1,mrp)
          endif
          duu = sqrt(un1(1)**2 + un1(2)**2 + un1(3)**2)
          slop(n,1) = un1(1) / duu
          slop(n,2) = un1(2) / duu
          slop(n,3) = un1(3) / duu

C --- Check unit vector along the wake mr-(mr+1) line at n = i

          un2(1) = xw(n,mrp) - xw(n,mr)
          un2(2) = yw(n,mrp) - yw(n,mr)
          un2(3) = zw(n,mrp) - zw(n,mr)

          duu = sqrt(un2(1)**2 + un2(2)**2 + un2(3)**2)
          slop2(n,1) = un2(1) / duu
          slop2(n,2) = un2(2) / duu
          slop2(n,3) = un2(3) / duu

C --- Find unit vector normal to slop & slop2 vector

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

        do n = 1 , ncvxp
          xvcc(n) = rc1 *fact(n) * slopd(n,1) + xw(n,mrp)
          yvcc(n) = rc1 *fact(n) * slopd(n,2) + yw(n,mrp)
          zvcc(n) = rc1 *fact(n) * slopd(n,3) + zw(n,mrp)
        enddo

C -- Local coordinate s1,s2 & s3 coordinate along the wake line
C -- unit vector of s1 relative to xyz --> slop(n,k)
C --                s2             xyz --> slop3(n,k)
C --                s3             xyz --> slopd(n,k)

C -- Set circle with radius r at the s2 & s3 coordinate

        do n = 1 , ncvxp
          do m = 1 , mcvt
             dtt2 = dtt * real(m-1)
             ss2 = rc1 * fact(n) * (1.-cos(dtt2))
             ss3 = rc1 * fact(n) * sin(dtt2)

             dxx = slop3(n,1)*ss2 + slopd(n,1)*ss3
             dyy = slop3(n,2)*ss2 + slopd(n,2)*ss3
             dzz = slop3(n,3)*ss2 + slopd(n,3)*ss3

             xvc(n,m) = xw(n,mrp) + dxx
             yvc(n,m) = yw(n,mrp) + dyy
             zvc(n,m) = zw(n,mrp) + dzz
          enddo
        enddo

        do n = 1 , ncvxp
          xvc(n,mcvtp) = xvc(n,1)
          yvc(n,mcvtp) = yvc(n,1)
          zvc(n,mcvtp) = zvc(n,1)
        enddo

      call btipgeo
      
        RETURN
        END


