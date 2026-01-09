      SUBROUTINE CAVVOL
************************************************************************
*     Compute the CAVity VOLume                                        *
*     05-14-92 NF                                                      *
*   Date of revision        Revision                                   *
*   ----------------        --------                                   *
*   05-18-92 NF     -found two bugs: deltav was not zeroed every strip *
*                    and A1 was written instead of A2 inside DO 70...  *
*   11-15-98 JY     -modified subroutine to obtain a more accurate     *
*                    estimate of cavity volume.                        *
*   03-05-99 JY     -Modified subroutine to allow cavity to grow on    *
*                    both the back and face of the foil.               *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
*   09-17-99 JY     -Modified the format of *.vol so that the user can *
*                    plot the contribution of face, back, and super-   *
*                    cavitation.                                       *
*   07-12-00 JY     -Modified subroutine to obtain a more accurate     *
*                    estimate of cavity volume.                        *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      CVOL=ZERO
C-----------------------------------------------------------------------
C     Print cavity length to the screen
C-----------------------------------------------------------------------
       ISR=1
       IF(IFACE.EQ.2) ISR=2
       DO 10 II=1,ISR
          IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN        
             IDR=1
             K=1
             ISF=0
             CVOLB=ZERO
          ELSE
             IDR=2
             K=-1
             ISF=1
             CVOLF=ZERO
          END IF
          
          IF(IDR.EQ.1) THEN
             WRITE(*,'(a4,1x,10(F6.3,1X))') 'Back',(CAVL(M,IDR),M=1,MR)
          ELSE IF(IDR.EQ.2) THEN
             WRITE(*,'(a4,1x,10(F6.3,1X))') 'Face',(CAVL(M,IDR),M=1,MR)
          END IF

C-----------------------------------------------------------------------
C     Compute the cavity volume at the current timestep
C-----------------------------------------------------------------------
          DO 20 M=1,MR
             DELA1=ZERO

             M0M=M0(M,IDR)
             M0M1=M0M-1

             DO 30 N=1,JCV(M,IDR)+NSPP(M,IDR)
                N1=M0M1+K*N+ISF
                L=INDEXB(N1,M)
                IF(N.LE.JCV(M,IDR)) THEN
                   DSN=SS(L,1)
                   DHM=HALF*(HT(N+1,M,IDR)+HT(N,M,IDR))
                ELSE 
                   DSN=SS(L,1)*FLP(M,IDR)
                   DHM=HALF*HT(N,M,IDR)
                END IF
                DELA1=DELA1+DHM*DSN
 30          CONTINUE
             CVOL=CVOL+DELA1
             IF(IDR.EQ.1) CVOLB=CVOLB+DELA1
             IF(IDR.EQ.2) CVOLF=CVOLF+DELA1
 20       CONTINUE

 10    CONTINUE

       CVOLW=ZERO
       DO 40 M=1,MR
          DELA1=ZERO
          IF(NNWC(M).GT.0) THEN
             DO 50 N=1,NNWC(M)+NSPS(M,NWDIR(M))
                N1=JCV(M,NWDIR(M))+N
                L=NTRA*(MR-M)+N
                IF(N.LE.NNWC(M)) THEN
                   DSN=SSW(L,1)
                   DHM=HALF*(HT(N1+1,M,NWDIR(M))+HT(N1,M,NWDIR(M)))
                ELSE 
                   DSN=SSW(L,1)*FLS(M,NWDIR(M))
                   DHM=HALF*HT(N1,M,NWDIR(M))
                END IF
                DELA1=DELA1+DHM*DSN
 50          CONTINUE
          END IF
          CVOL=CVOL+DELA1
          CVOLW=CVOLW+DELA1
 40    CONTINUE

       THETA=-FLOAT(NTSTEP-1)*DELTAT
       ANGOUT=THETA*(-180)/PI

       IF(ISTEADY.EQ.0) THEN
          WRITE(58,2000) NTSTEP,CVOL,CVOLB,CVOLF,CVOLW
       ELSE
          WRITE(58,1000) ANGOUT,CVOL,CVOLB,CVOLF,CVOLW
       END IF
 1000  FORMAT(1X,F8.2,2X,E13.6,2X,E13.6,2X,E13.6,2X,E13.6)
 2000  FORMAT(1X,I3,2X,E13.6,2X,E13.6,2X,E13.6,2X,E13.6)
!Allen Du 01/15/2018 output the cavity volume
       if(opt_vol_flag==1.and.opt_vol<cvol) opt_vol=cvol
      
       RETURN
C<<<<<<<<<<<<<<<<<<<<end of subroutine CAVVOL>>>>>>>>>>>>>>>>>>>>>>>>>>>
       END
