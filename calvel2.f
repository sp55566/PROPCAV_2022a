************************************************************************
      SUBROUTINE CALVEL2
************************************************************************
*     CALVEL2: Velocity and coordinate trajectory calculation          *
*                  at the Tip Vortex cavity Center                     *
*                                                                      *
*  !!! WARNING : This subroutine has to be compiled with Fortran 90    *
*  !!!           (You can't compile with Fortran 77)                   *
*                                                                      *
*  Date of last Revision       Revision                                *
*  ---------------------       --------                                *
*  10-12-99   HSLee   Calculate total velocity on the tip vortex center*
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

        DO 250 N = 1 , NCVX
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
          DO 251 M = 1 , MCVT 
           L = INDEXC(N,M)
           UXI = DPDUC(N,M)
           UETA = (DPDVC(N,M) - DPDUC(N,M)*SINPHI(L))/COSPHI(L)
           uxctot(n,m) = (vxic(n,m)+uxi) * dir(l,1,1)+
     &                  (vetac(n,m)+ueta)*dir(l,2,1)
           uyctot(n,m) = (vxic(n,m)+uxi) * dir(l,1,2)+
     &                  (vetac(n,m)+ueta)*dir(l,2,2)
           uzctot(n,m) = (vxic(n,m)+uxi) * dir(l,1,3)+
     &                  (vetac(n,m)+ueta)*dir(l,2,3)
         sum1 = sum1 + uxctot(n,m) 
         sum2 = sum2 + uyctot(n,m)
         sum3 = sum3 + uzctot(n,m)
251       continue      
       vxmean(n) = sum1/real(mcvt)
       vymean(n) = sum2/real(mcvt)
       vzmean(n) = sum3/real(mcvt)
250      CONTINUE

        do n = 1 , ncvx

          cmx(n)= 0.0
          cmy(n)= 0.0
          cmz(n)= 0.0

          do m = 1 , mcvt
           nm = indexc(n,m)
             cmx(n) = cmx(n) + xct(nm,1) 
             cmy(n) = cmy(n) + xct(nm,2)
             cmz(n) = cmz(n) + xct(nm,3)
          enddo

        cmx(n) = cmx(n) / real(mcvt)
        cmy(n) = cmy(n) / real(mcvt)
        cmz(n) = cmz(n) / real(mcvt)

        enddo

        RETURN
        END
  
