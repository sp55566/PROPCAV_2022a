C --------------------------------
       FUNCTION INDEXW2(N,M)
C --------------------------------
C       INDEXW2: INDEX of the Wake panels of "GWAKETIP"
C
        INCLUDE 'PUFCAV.INC'

        if((ian.ne.6).and.(ian.ne.2)) then
           indexw2=mr*(n-1)+m
        else
           indexw2 = (mr-m) * nwpanel + n
        end if

        RETURN
        END
