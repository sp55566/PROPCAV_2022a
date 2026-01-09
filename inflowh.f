      SUBROUTINE INFLOWH
************************************************************************
*     INFLOW: INFLOW velocities                                        *
*      --- Compute the source strength of each blade at current time   *
*            step                                                      *
*      --- Compute the total inflow velocities at the key blade        *
*     ---obtained from gbflow at the control points --Shreenaath *
*                                                                      *
* Date        Revision                                                 *
* ----        --------                                                 *

************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

C-----------------------------------------------------------------------
C         If running hydrofoil case, with total inflow velocities obtained from gbflow
C-----------------------------------------------------------------------
          if((icon.eq.5).OR.(ICON.EQ.6).OR.(ICON.EQ.8))then
             open(unit=20,file='S98.DAT',status='unknown')

             DO J=1,NPANEL
                READ(20,*)VOX1(J),VOY1(J),VOZ1(J)
             ENDDO  
          endif
             close(20)

             DO 100 KK = 1 , NBLADE
                             DTBLA=DELK*FLOAT(KK-1)
            DO 40 J=1,NPANEL
               VOZ1(J)=VOZ1(J)+tan(ALPHA*rad)
           VELY=VEL(J,2)*COS(DTBLA)-VEL(J,3)*SIN(DTBLA) 
           VELZ=VEL(J,2)*SIN(DTBLA)+VEL(J,3)*COS(DTBLA) 
           BUG=VOX1(J)*VEL(J,1)+VOY1(J)*VELY+VOZ1(J)*VELZ

          if(kk.eq.1) then
            kk0=1
            DPDN(J)=-BUG
          else
            kk0=nblade+2-kk
          end if
          STRGTH(J,KK0)=-BUG

   40   CONTINUE

        IF(KK.EQ.1) THEN
C-----------------------------------------------------------------------
C           On-coming velocity in local coordinate system
C-----------------------------------------------------------------------
          DO 70 N=1,NC
            DO 60 M=1,MR
               L=INDEXB(N,M)
               VXIB(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
     *                     +VOZ1(L)*DIR(L,1,3)
               VETAB(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
     *                      +VOZ1(L)*DIR(L,2,3)
               VINFSB(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
   60       CONTINUE
   70     CONTINUE

          IF(IHUB.NE.0) THEN
             DO 90 N=1,NHBX
                DO 80 M=1,MHBT
                   L=INDEXH(N,M)
                   VXIH(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
     *                  +VOZ1(L)*DIR(L,1,3)
                   VETAH(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
     *                  +VOZ1(L)*DIR(L,2,3)
                   VINFSH(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
 80             CONTINUE
 90          CONTINUE
          END IF
       END IF
 100  CONTINUE

C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
cC -- Begin Tip HSLEE(10/12/99)
c
c      If(IAN .eq. 2) then
c
c             DO N=1,NTHX
c               DO M=1, MCVT
c                 L=INDEXT(N,M)
c                 VXITH(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
c     *                     +VOZ1(L)*DIR(L,1,3)
c                 VETATH(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
c     *                      +VOZ1(L)*DIR(L,2,3)
c                 VINFSTH(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
c               ENDDO
c             ENDDO
c
c             DO N=1,NCVX
c               DO M=1, MCVT
c                 L=INDEXC(N,M)
c                 VXIC(N,M)=VOX1(L)*DIR(L,1,1)+VOY1(L)*DIR(L,1,2)
c     *                     +VOZ1(L)*DIR(L,1,3)
c                 VETAC(N,M)=VOX1(L)*DIR(L,2,1)+VOY1(L)*DIR(L,2,2)
c     *                      +VOZ1(L)*DIR(L,2,3)
c                 VINFSC(N,M)=VOX1(L)**2+VOY1(L)**2+VOZ1(L)**2
c               ENDDO
c             ENDDO
c
c      Endif
C/e S.N.KIM | Aug. 2018.

cC -- End Tip (10/12/99)

      RETURN
C>>>>>>>>>>>>>>>>>>>>>>End of subroutine INFLOW>>>>>>>>>>>>>>>>>>>>>>>>>
      END
