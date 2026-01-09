LN:1914
C----- usual logarithmic differences
       XLOG = DLOG(X2/X1)
!s--- YE TIAN ----
       if((U2.lt.0.0).or.(U1.le.0.0)) write(*,*) U2,U1, 'xbl.f:1911'
       U2=ABS(U2)
       U1=ABS(U1)
!e--- YE TIAN ----
       ULOG = DLOG(U2/U1) !--------> A
       TLOG = DLOG(T2/T1)
       HLOG = DLOG(HS2/HS1)
C       XLOG = 2.D0*(X2-X1)/(X2+X1)
C       ULOG = 2.D0*(U2-U1)/(U2+U1)
C       TLOG = 2.D0*(T2-T1)/(T2+T1)
C       HLOG = 2.D0*(HS2-HS1)/(HS2+HS1)
       DDLOG = 1.D0
      ENDIF

LN:2317
!s-- YE TIAN 6/10/2013-----
      if (RT.lt.0.0) then
        write(*,*) 'RT=',RT,'xbl.f:2319'
        RT = abs(RT)
      end if
!e-- YE TIAN 6/10/2013-----
      GR = DLOG10(RT)  !------------> B
      GR_RT = 1.D0 / (2.3026D0*RT)


LN:2552
C---- Turbulent HS correlation   ( from Swafford profile family )
!s--YE TIAN ----
      RT = abs(RT)
!e--YE TIAN ----
      GRT = DLOG(RT)

LN:1452
C---- equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
      HKB = HK2 - 1.D0
      USB = 1.D0 - US2
!s--YE TIAN -- 06/10/2013----
      if (CTCON*HS2*HKB**3 / (USB*H2*HK2**2).lt.0.0) then
        write(*,*) CTCON, HS2, HKB, USB, H2, HK2, 'xbl.f:1456'
        HS2 = abs(HS2)
      end if
!e--YE TIAN -- 06/10/2013----
      CQ2     =
     &    DSQRT( CTCON*HS2*HKB**3 / (USB*H2*HK2**2) )
