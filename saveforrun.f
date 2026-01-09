C =================================
      SUBROUTINE SAVEFORRUN
C =================================

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      OPEN(755,FILE='SAVE1.DAT',STATUS='UNKNOWN')
      OPEN(756,FILE='SAVE2.DAT',STATUS='UNKNOWN',
     %     FORM='UNFORMATTED')

      DO K = 1 , 2
         DO M = 1 , MR
            DO N = 1, NTPREV
               WRITE(755,*) CAVLSA(M,N,K),NLEP(M,N,K),
     %              NOCAV(M,N,K)
            ENDDO
         ENDDO
      ENDDO

      DO N = 1 , NPANZ
         CALL WRITE1(756,POTM(1,N),NSTEP)
      ENDDO

      CLOSE(755)
      CLOSE(756)

      RETURN
      END


C =================================
      SUBROUTINE READFORRUN
C =================================

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

      OPEN(755,FILE='SAVE1.DAT',STATUS='UNKNOWN')
      OPEN(756,FILE='SAVE2.DAT',STATUS='UNKNOWN',
     %     FORM='UNFORMATTED')

      DO K = 1 , 2
         DO M = 1 , MR
            DO N = 1, NTPREV
               READ(755,*) CAVLSA(M,N,K),NLEP(M,N,K),
     %              NOCAV(M,N,K)
            ENDDO
         ENDDO
      ENDDO

      DO N = 1 , NPANZ
         CALL READ1(756,POTM(1,N),NSTEP)
      ENDDO

      CLOSE(755)
      CLOSE(756)

      RETURN
      END
