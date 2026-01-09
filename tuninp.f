C ======================================================================
      SUBROUTINE TUNINP
C
C     Tunnel geometry with circular cross section
C     
C     ITUNGEO=1  Straight Panels on tunnel walls 
C                (not recommended if TUNRAD < 1.5)
C     ITUNGEO=2  Panels using blade pitch angle on tunnel walls 
C                (recommended if TUNRAD is close to 1.0)
C ======================================================================
      INCLUDE 'PUFCAV.INC'
      
c     read tunnel data
c
      OPEN(9,FILE='tunnel.dat',STATUS='OLD')

      
      READ(9,*) ITUNGEO, TUNRAD
      WRITE(*,*) ' CIRCULAR TUNNEL WALL  RADIUS / RPROP: ', TUNRAD

      READ(9,*) XFOR, XAFT
      WRITE(*,*) ' XFOR, XAFT FOR THE TUNNEL:     ', Xfor, Xaft

      READ(9,*) Xmidf, xmida
      WRITE(*,*) ' XMIDF, XMIDA (Mid Part)  :     ', Xmidf,xmida

      READ(9,*) NAXf, NAXm, NAXa, MTUNEL, NSIDE

      WRITE(*,*) ' No. of wall grid (NAXf,_m,_a): ', Naxf,Naxm,Naxa
      WRITE(*,*) ' No. of circumferential grid :  ', MTUNEL
      WRITE(*,*) ' No. of radial grid on LID   :  ', NSIDE

      CLOSE(9)

      RETURN
      END
