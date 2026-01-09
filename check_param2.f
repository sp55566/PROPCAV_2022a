C ==================================
       SUBROUTINE CHECK_PARAM2
C ==================================

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       dimension ierh(10),ierw(10)

       ierh = 0
       ierw = 0

       icount = 3

       if(ihub .gt. 0) then

          icount = icount + 1

          CALL CLEAR(IERH,10)

          write(99,*)
          write(99,*) ' --------------------------------------'
          write(99,411) icount
 411      format('   (',i1,') Checking....  HUB Part          ')
          write(99,*) ' --------------------------------------'
          write(99,*)
          
          if(MHBT .gt. mhbz) ierh(1) = 1
          write(99,110) mhbz , mhbt,ierh(1)

          if(nblade*mhbz .lt. 20) then
             write(99,111) nblade , mhbz 
          endif

          if(NHBU .gt. nhuz) ierh(2) = 1
          write(99,115) nhuz , nhbu,ierh(2)

          if(NHBDT .gt. nhdz) ierh(3) = 1
          write(99,120) nhdz , nhbdt,ierh(3)
          if(ierh(3) .eq. 1) write(99,121) 

          if(rhub .ne. xr(1)) then
             write(99,125) xr(1), rhub
          end if

          x1 =  xhbd + xhbt + xb(1,1)
          if(x1 .ge. xw(nsw(1),1)) then
             ierh(4) = 1
             write(99,130) xw(nsw(1),1), x1 , ierh(4)
          endif

C          if(RULT .gt. 1.0) then
C             ierh(5) = 1
C             write(99,135) RULT 
C          endif

C          write(*,*) ' rhult =', rhult, rhub
C
C          if(RHULT .gt. RHUB) then
C             ierh(6) = 1
C             write(99,140) RHULT 
C          endif



       endif
      
 110   format(' Max mhbt   = ',i3,'    Input mhbt   = ',i3,
     %      '    Err = ',i1) 
 111   format(' Nblade     = ',i3,'    MHBT         = ',i3 /
     %      '          <==== WARNING ===> '/
     %      ' Too little number of panels between blades') 
 115   format(' Max nhbu   = ',i3,'    Input nhbu   = ',i3,
     %      '    Err = ',i1) 
 120   format(' Max nhbdt  = ',i3,'    Input nhbdt  = ',i3,
     %      '    Err = ',i1) 
 121   format(' Warning! -- Check downstream hub length ----- ')
 125   format(' INPUT XR(1)= ',f6.4,'  INPUT RHUB   = ',f6.4,
     %      '          <==== WARNING ===> '/
     *      ' Extrapolate to determine RZ(1)')
 130   format(' XW_end     = ',f6.4,'  XH at HUB END= ',f8.4,
     %      '  Err = ',i1 /
     %      '            ==== ERROR ====  '/
     %      '      Change XHBD and XHBT Value  '/
     %      '                    or                    '/
     %      '      Change XULT and XUWDK ') 
 135   format(' RULT = ',f6.3,' ***** RULT =< 1.0 ','   Err = 1') 
 140   format(' RHULT = ',f6.3,' ***** RHULT =< RHUB ','     Err = 1') 


       CALL CLEAR(IERW,10)
       
C/s S.N.KIM | Tip vortex model is omitted in PROPCAV released in 2018.
c       IF(IAN.EQ.2) THEN
c          icount = icount + 1
c
c          write(99,*)
c          write(99,*) ' --------------------------------------'
c          write(99,412) icount
c 412      format('   (',i1,') Checking.... Tip Vortex Part    ')
c          write(99,*) ' --------------------------------------'
c          write(99,*)
c
c          if(mcvt .gt. mcavm) ierw(1) = 1
c          write(99,250) mcavm, mcvt , ierw(1)
c          
c 250      format(' Max MCVT   = ',i3,'    Input MCVT = ',i3,
c     %         '    Err = ',i1) 
c       END IF
C/e S.N.KIM | Aug. 2018.

       icount = icount + 1
       write(99,*)
       write(99,*) ' --------------------------------------'
       write(99,413) icount
 413   format('   (',i1,') Checking.... Time Step Size     ')
       write(99,*) ' --------------------------------------'
       write(99,*)

       NREC = 360 / ndltat
       IF(NREC.GT.NSTEP) ierw(2)=1
       WRITE(99,310) NSTEP,NREC,IERW(2)
       NDUM = 360 / NDLTAT
       IF(NDUM*NDLTAT .NE. 360) NDUM = NDUM + 1 

       IF(IERW(2).EQ.1) WRITE(99,315) NREC, NDUM

 310   FORMAT(' Max NSTEP     = ',i3,'    Calculated NSTEP= ',i3,
     %      '    Err = ',i1) 
 315   FORMAT('            ===== ERROR ======  '/
     *        '      Change NSTEP in PARAM.INC to ',i3,/
     *        '                   OR               '/
     *        '      Change DELTAT = ',I3 )

       icount = icount + 1
       write(99,*)
       write(99,*) ' --------------------------------------'
       write(99,414) icount
 414   format('   (',i1,') Checking.... Cavitating Analysis')
       write(99,*) ' --------------------------------------'
       write(99,*)

       IF(NCTIME.GT.0) THEN
          IF(NTZ.LE.NPANZ) IERW(3)=1
          WRITE(99,320) 
 320      FORMAT(' NCTIME > 0 ==> Please set NSCWZ=MBZ*NWZ ')
       END IF

       isum = 0
       do i = 1 , 10
          isum = isum + ierh(i)
          isum = isum + ierw(i)
       enddo
       if(isum .gt. 0) then
          write(*,*)
          write(*,*)
          write(*,*) '********** ERROR in INPUT DATA File ! ***********'
          write(*,*) '*                                               *'
          write(*,*) '*       Check <ERR.LOG> File to check ERROR     *'
          write(*,*) '*                                               *'
          write(*,*) '*************************************************'
          write(*,*)
          stop
       endif
       CLOSE(99)
       
       RETURN
       END
