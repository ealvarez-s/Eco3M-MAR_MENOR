      program notebook_rivers_add_line
      implicit none
      integer nriver,loop,loop1,loop2,loop3
      character*100 ligne1,ligne2
      loop2=13
      loop3=4

      ligne2='0     ! Boundary condition scheme at river grid input cell          RIVER_OBCTYPE(:)'
!'
      open(unit=3,file='notebook_rivers')
      open(unit=4,file='notebook_rivers_new')

! 0     ! Number of river inputs

       read(3,*)nriver
      write(4,*)nriver,' ! Number of river inputs'

      do loop1=1,nriver

        do loop=1,loop2
          read(3,'(a)')ligne1
         write(4,'(a)')ligne1
        enddo
        write(4,'(a)')ligne2
        do loop=1,loop3
          read(3,'(a)')ligne1
         write(4,'(a)')ligne1
        enddo

      enddo ! loop1


      close(3)
      close(4)


      end
