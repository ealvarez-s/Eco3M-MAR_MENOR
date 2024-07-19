      program toto
      implicit none
      integer :: i=1
      character*50 texte

      open(unit=3,file='notebook_rivers')

       i=1
       read(3,'(i4,a50)')i,texte
       write(6,*)i,'   ',texte

       i=1
       read(3,*)
       read(3,*)

       read(3,'(i4,a50)')i,texte
       write(6,*)i,'   ',texte

       read(3,'(i4,a50)')i,texte
       write(6,*)i,'   ',texte

      close(3)

      write(6,*)i,'   ',texte

      end
