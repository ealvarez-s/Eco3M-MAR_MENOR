      SUBROUTINE LightParameters
!______________________________________________________________________
! 3d ecosystem model
!
! LAST REVISION: 23 MARCH 2008
!
! Implementation: Caroline Ulses
!                 NIOO-CEME
!______________________________________________________________________
 
 
!______________________________________________________________________
!
! 1. Reads the light parameters in a notebook
!   1.a Zooplankton
!   1.b Phytoplankton
!   1.c Remineralisation
! 3. Writes them in an output file
!
!______________________________________________________________________
 
!______________________________________________________________________
!
! Modifications:
! 21/02/07: include symphonie2007.h dans la version 2004_v7
! 23/03/08: Passage à SYMPHONIE2008 (Module principal)
!
!______________________________________________________________________
 
 
      USE ModuleComBlock
      USE MODULE_PRINCIPAL
      use module_parallele
      IMPLICIT NONE
 
      INTEGER LDIR
 
!======================================================================
! 1. Reading of the light parameter notebook: notebook_light
! BEGINNING:
!======================================================================
 
 
      if(par%rank==0)WRITE(6,'(A,A)')'Ready to read notebook_light'            &
                        ,NOMFICHIER(23)
 
      OPEN(UNIT=3,FILE=NOMFICHIER(23)) ! open the notebook
 
      DO K=1,11
      READ(3,*)
      ENDDO
 
      READ(3,*)AlbedoBio
      READ(3,*)pHeat
      READ(3,*)pLONG
 
      READ(3,*)
 
      READ(3,*)kro
      READ(3,*)kbo
      READ(3,*)krp
      READ(3,*)kbp
      READ(3,*)lr
      READ(3,*)lb
#ifdef mes
      READ(3,*)kbn
#endif
 
!======================================================================
! 1. Reading of the light  parameter notebook: notebook_light
! END:
!======================================================================
 
 
!======================================================================
! 2. Writing of the light parameters in an output file
! BEGINNING :
!======================================================================
      if(par%rank==0)print*,'Avant nomfichier LightParameters'
      TEXTE90=nomfichier(21)
      OPEN(UNIT=3,FILE=TEXTE90)
      READ(3,*)
      READ(3,*)
      READ(3,'(A)')DIRGRAPH                                            !13/04/06
      CLOSE(3)
 
      DO 20 K=1,90
        IF(DIRGRAPH(K:K).EQ.' ') THEN
        LDIR=K-1
        GOTO 21
        ENDIF
   20 CONTINUE
   21 CONTINUE
 
      TEXTE90=DIRGRAPH(1:LDIR)//'/light_parameters'
      if(par%rank==0)WRITE(6,*)'Writing of the parameters in the directory'         &
                ,TEXTE90
 
      if(par%rank==0)then
      OPEN(unit=2,FILE=TEXTE90)
 
      WRITE(2,*)
 
      WRITE(2,*)'AlbedoBio     =       ',AlbedoBio,'    -'
      WRITE(2,*)'pHeat         =       ',pHeat ,'    -'
      WRITE(2,*)'pLONG         =       ',pLONG ,'    -'
      WRITE(2,*)'kro           =       ',kro   ,'    -'
      WRITE(2,*)'kbo           =       ',kbo   ,'    -'
      WRITE(2,*)'krp           =       ',krp   ,'    -'
      WRITE(2,*)'kbp           =       ',kbp   ,'    /mmolN/m'
      WRITE(2,*)'lr            =       ',lr    ,'    -'
      WRITE(2,*)'lb            =       ',lb    ,'    -'
#ifdef mes
      WRITE(2,*)'kbn           =       ',kbn   ,'    -'
#endif
      CLOSE(2)
      endif
 
!======================================================================
! 3. Writing of the light parameters in an output file
! END:
!======================================================================
 
      END SUBROUTINE LightParameters
