      FUNCTION IDXWAK(LIN,MIN)
      use m_WKNP
      INCLUDE 'PUFCAV.INC'
!     COMMON /WKNP/ NWIDX(MBZ),NWSEC(MBZ)
      IDXWAK=NWIDX(MIN)+LIN
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of function IDXFWAK>>>>>>>>>>>>>>>>>>>>>>>>>>
      END
