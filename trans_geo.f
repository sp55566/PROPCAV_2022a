C ===============================================
      SUBROUTINE TRANS_GEO1(IP)
C ===============================================
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVC.INC'

      DO N = 1 , NCP
         DO M = 1 , MRP
            CALL ROTATE2(IP,SPANGLE,XB(N,M),YB(N,M))
         ENDDO
      ENDDO

      IF(IHUB .NE. 0) THEN
         DO N = 1, NHBX+1
            DO M = 1 , MHBT+1
               CALL ROTATE2(IP,SPANGLE,XH(N,M),YH(N,M))
            ENDDO
         ENDDO
         IF(IHUB .EQ. 1 .OR. IHUB .EQ. 3 .OR. 
     %      IHUB .EQ. 5 .OR. IHUB .EQ. 7) THEN
            DO M = 1, 2
               DO N = 1 , MHBT+1
                  CALL ROTATE2(IP,SPANGLE,XHDK(N,M),YHDK(N,M))
               ENDDO
            ENDDO   
         ENDIF

         IF(IHUB .EQ. 3 .OR. IHUB .EQ. 4 .OR. IHUB .EQ. 5) THEN
            DO M = 1, 2
               DO N = 1 , MHBT+1
                  CALL ROTATE2(IP,SPANGLE,XHDKDT(N,M),YHDKDT(N,M))
               ENDDO
            ENDDO   
         ENDIF

      ENDIF

      IF(IDUCT .NE. 0) THEN
         DO N = 1, NDUCT+1
            DO M = 1 , MDUCT+1
               CALL ROTATE2(IP,SPANGLE,XD(N,M),YD(N,M))
            ENDDO
         ENDDO
      ENDIF
      
      RETURN
      END

C ===============================================
      SUBROUTINE TRANS_GEO2(IP)
C ===============================================
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVC.INC'

      DO M = 1, MRP
         DO N = 1 , NSW(M)
            CALL ROTATE2(IP,SPANGLE,XW(N,M),YW(N,M))
         ENDDO
      ENDDO            
      
      DO M = 1 , 2
         DO N = 1 , MHBT+1
            CALL ROTATE2(IP,SPANGLE,XWDK(N,M),YWDK(N,M))
         ENDDO
      ENDDO            
      
      DO M = 1 , MRP
         DO N = 1 , NSUB
            DO L = 1 , NWSUB1
               NIDX=(N-1)*NWSUB1+L+1
               CALL ROTATE2(IP,SPANGLE,XWS(NIDX,M),YWS(NIDX,M))
            ENDDO
         ENDDO      
      ENDDO
      
      RETURN
      END
