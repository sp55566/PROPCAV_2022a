C ------------------------------------------
       subroutine wtprs
C ------------------------------------------
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'

       ittt =itstep
       do m = 1 , mcvt
          write(21,230) ittt,m
          sum = 0.0
          do nn = 1 , ncvx
             nm = indexc(nn,m)
             sum = sum + delu(nm)
          enddo
          
          dum1 = 0.0
          do nn = 1 , ncvx
             if(nn .eq. 1) then
                nm = indexc(nn,m)
                dum1 = dum1 + half * delu(nm)
             else
                nm = indexc(nn,m)
                nm1 = indexc(nn-1,m)
                dum1 = dum1 + half *(delu(nm) + delu(nm1))
             endif
             dum = dum1 / sum
             
             write(21,240) dum,-cpc(nn,m)*advco**2,
     %            vtc(nn,m),pot(nm)
          enddo
       enddo
       
 230   FORMAT('ZONE T = "T=',i3,' M= ',i3,'"')
 240   FORMAT(4(1X,F16.6))
       
       return
       END


