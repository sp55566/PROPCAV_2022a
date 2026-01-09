      SUBROUTINE DUVLC(mm)

C**********************************************************************
C     This subroutine is used to update the inviscid velocity component
c     when coupled with 2D boundary layer solver XFOIL (fully wetted)
C     *
C     By XM YU                          Dec, 2011                  *
C**********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      INCLUDE 'PUFBL.INC' 

      DIMENSION dtemp(nzvis)

        NTW(mm)=nc+nwsub-nsub+nwmin

        NWv = NTW(mm)-NC

        do im=1,nc
           do in=1,nc
            av(in,im)=avl(in,im,mm)
           enddo 
        enddo

        do i=1,nwv
           do j=1,nc
              awv(i,j)=awv3(i,j,mm)
           enddo
        enddo

        CALL INVRSE(AV,AVINV,NC,NZVIS) !AVINV

      do i=1,ntw(mm)-1 
         duvl(i,mm)=0.0
      enddo
     
      ca = 0.0
      cb = 0.0
      dh1= 0.0
      DO 1000 M = 1, MR
        ntw(m)=nc+nwsub-nsub+nwmin
       
        beta = 0.0
        csigma = 0.0
        if ((m.ne.mm).and.(lvflag(m).ne.1)) then
           do ni=1,nc
             i=indexb(ni,mm)
             do nj=1,nc
                j=indexb(nj,m)
                beta(ni,nj)=bsbb(i,j)
             enddo
           enddo

           do ni=1,nc
             i=indexb(ni,mm)
             do nj=1,nwv
c                j=idxwak(nj,m)
                j=(mr-m)*nwv+nj
                beta(ni,nc+nj)=bsbw(i,j)
             enddo
           enddo

           do ni=1,nwv
c              i=idxwak(ni,mm)
              i=(mr-mm)*nwv+ni
              do nj=1,nc
                 j=indexb(nj,m)
                 csigma(ni,nj)=bswb(i,j)
              enddo
           enddo
                      
           do ni=1,nwv
c              i=idxwak(ni,mm)
              i=(mr-mm)*nwv+ni
              do nj=1,nwv
c                 j=idxwak(nj,m)
                 j=(mr-m)*nwv+nj
                 csigma(ni,nc+nj)=bsww(i,j)
              enddo
           enddo
        endif   

        if ((m.eq.mm).and.(lvflag(m).ne.1)) then
           do ni=1,nc
             i=indexb(ni,mm)
             do nj=1,nc
                j=indexb(nj,m)
                beta(ni,nj)=bsbb(i,j)-bv(ni,nj,mm)
             enddo
           enddo

           do ni=1,nc
             i=indexb(ni,mm)
             do nj=1,nwv
c                j=idxwak(nj,m)
                j=(mr-m)*nwv+nj
                beta(ni,nc+nj)=bsbw(i,j)-bv(ni,nc+nj,mm)
             enddo
           enddo

           do ni=1,nwv
c              i=idxwak(ni,mm)
              i=(mr-mm)*nwv+ni
              do nj=1,nc
                 j=indexb(nj,m)
                 csigma(ni,nj)=bswb(i,j)-csig(ni,nj,mm)*4*PI
              enddo
           enddo
                      
           do ni=1,nwv
c              i=idxwak(ni,mm)
              i=(mr-mm)*nwv+ni
              do nj=1,nwv
c                 j=idxwak(nj,m)
                 j=(mr-m)*nwv+nj
                 csigma(ni,nc+nj)=bsww(i,j)-csig(ni,nc+nj,mm)*4*PI
              enddo
           enddo
        endif   
c XM YU 06/2012
        nstag=nc/2+1
c XM YU 06/2012

        CALL MAT2DOT(AVINV,BETA,CA,NC,NZVIS,NC,NZVIS,NTW(M),NZVIS)
        
        DO 140 I=2,NC
          DO 130 J=1, NTW(M)
            DH(I,J)=(CA(I,J)-CA(I-1,J))/(SP(I,mm)-SP(I-1,mm))
            IF(I.LT.NSTAG)THEN
              DH(I,J)=-DH(I,J) 
            ENDIF
 130      CONTINUE
 140    CONTINUE


        CALL MAT2DOT(AWV,AVINV,CB,NWv,NWZ,NC,NBZ,NC,NBZ)

        CALL MAT2DOT(CB,BETA,DH1,NWv,NWZ,NC,NBZ,NTW(M),NZVIS)

        DO 160 I=1,NWv
          DO 150 J=1,NTW(M)
c            DH1(I,J)=(DH1(I,J)+CSIGMA(I,J))/2./PI
            DH1(I,J)=(DH1(I,J)+CSIGMA(I,J))/4.0/PI
 150      CONTINUE
 160    CONTINUE
        
        DO 180 I=2,NWv
          DO 170 J=1,NTW(M)
            DH(NC+I,J)=(DH1(I,J)-DH1(I-1,J))/(SP(NC+I,mm)-SP(NC+I-1,mm))
 170      CONTINUE
 180    CONTINUE

        DO 190 J=1,NTW(M)
          DH(1,J)   = DH(2,J)+(DH(2,J)-DH(3,J))*(SC(1,mm)-SC(2,mm))
     &      / (SC(2,mm)-SC(3,mm))
          DH(NC+1,J)= DH(NC,J)+(DH(NC,J)-DH(NC-1,J))
     &      * (SC(NC+1,mm)-SC(NC,mm))/(SC(NC,mm)-SC(NC-1,mm))
          IF(J.EQ.1)   DH(1,J)=-DH(1,J)
          IF(J.EQ.NC)  DH(NC+1,J)=-DH(NC+1,J)
          IF(J.LT.NSTAG) DH(NC+1,J)=-DH(1,J)
          IF(J.GE.NSTAG) DH(1,J)=-DH(NC+1,J)
 190    CONTINUE
       
        DO 210 I=1,NTW(M)-1
          DO 200 J=2,NTW(M)-1
            DD(I,J) = DH(I,J-1)/(SC(J,mm)-SC(J-1,mm)) 
     &        - DH(I,J)/(SC(J+1,mm)-SC(J,mm))
            IF(I.LT.NSTAG.AND.J.LE.NC+1) THEN
              DD(I,J) = -DD(I,J)
            ENDIF
            IF(J.GT.NC+1) DD(1,J)=-DD(1,J) 
 200      CONTINUE
 210    CONTINUE

        DO 220 I=1,NTW(M)-1
          IF(I.GT.NC+1)THEN
            DD(I,1)=DD(I,NC+1)
          ELSE
            DD(I,1)=DD(NC+2-I,NC+1)
          ENDIF
 220    CONTINUE
C     CALCULATE DVL 
        DO 250 I = 1,NTW(M)-1
          DO 240 J = 1,NTW(M)-1
            IF(I.LE.NC+1.AND.J.LE.NC+1)THEN
              DVL(I,J)=DD(NC+2-I,NC+2-J)
            ELSE IF(I.LE.NC+1.AND.J.GT.NC+1)THEN
              DVL(I,J)=DD(NC+2-I,J)
            ELSE IF(I.GT.NC+1.AND.J.LE.NC+1)THEN
              DVL(I,J)=DD(I,NC+2-J)
            ELSE IF(I.GT.NC+1.AND.J.GT.NC+1)THEN
              DVL(I,J)=DD(I,J)
            ENDIF
 240      CONTINUE
 250    CONTINUE

C     PUT TO VISCAL SIGN CONVENTION 
        DO 270 I=1,NTW(M)
          DO 260 J=1,NTW(M)
            IF(I.NE.NC+1) DVL(I,J)=-DVL(I,J)
            IF(I.GT.NC+1.AND.J.LE.NC) DVL(I,J)=-DVL(I,J)
            IF(J.GT.NC+1.AND.I.LE.NSTAG) DVL(I,J)=-DVL(I,J)
            IF(J.GT.NC+1.AND.I.EQ.NC+1) DVL(I,J)=-DVL(I,J)
 260      CONTINUE
 270    CONTINUE

        DO 280 J=1,NTW(M)
          DVL(NC+2,J)= DVL(NC+1,J)
 280    CONTINUE
        DVL(1,NC+1)=-DVL(1,NC+1)
        DVL(NC+1,NC+1)=-DVL(NC+1,NC+1)
        DVL(NC+2,NC+1)=-DVL(NC+2,NC+1)

        DO 290 J=NSTAG+1,NC
          DVL(1,J)=-DVL(1,J)
          DVL(NC+1,J)=-DVL(NC+1,J)
          DVL(NC+2,J)=-DVL(NC+2,J)
 290    CONTINUE

        do i=1,ntw(m)-1
           do j=1,ntw(m)-1
               dvl(i,j)=-dvl(i,j)
           enddo
        enddo
       
        if (lvflag(m).eq.1) then
          do j = 1, ntw(m)-1
            umass(j,m) = 0.0
          end do
        end if
        sum=0.0 
        do i=1,ntw(m)-1
           if (i.le.nc/2) then
              coei=1
           else
              coei=-1
           endif
           do j=1,ntw(m)-1
              if (j.le.nc/2) then
                 coej=1
              else
                 coej=-1
              endif    

!s--- YE TIAN --- 06/20/2013---
!             if(isnan(dvl(i,j))) then
              ! YE TIAN 09/01/2013 --- gfortran doesn't support isnan
              ! x.ne.x is an alternative for this function
              if(dvl(i,j).ne.dvl(i,j)) then  
                write(*,*) 'duvl.f:285:dvl',i,j
                stop
              end if
!             if(isnan(umass(j,m))) then
              if(umass(j,m).ne.umass(j,m)) then
                write(*,*) 'duvl.f:289:umass',i,j,m,mm
                write(*,*) 'lvflag(m)',lvflag(m)
                write(*,*) umass(j,m)
                write(*,*) dvl(i,j)
                write(*,*) dh(1,1), dh(2,1), dh(3,1)
                stop
              end if
!e--- YE TIAN --- 06/20/2013---
              sum=sum-coei*coej*dvl(i,j)*umass(j,m)
           enddo
           dtemp(i)=sum
           sum=0.0
        enddo
        do i =1,ntw(m)-1
           duvl(i,mm)=duvl(i,mm)+dtemp(i)
        enddo
       
 1000 CONTINUE                  !!     LOOP 1000 ENDED

      RETURN
      END
