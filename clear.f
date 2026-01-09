C ================================
      subroutine clear(a,n)
C ================================
      dimension a(n)

      do i = 1 , n
        a(i) = 0.0
      enddo

      return
      end

C ================================
      subroutine iclear(ia,n)
C ================================
      dimension ia(n)

      do i = 1 , n
        ia(i) = 0
      enddo

      return
      end


      
