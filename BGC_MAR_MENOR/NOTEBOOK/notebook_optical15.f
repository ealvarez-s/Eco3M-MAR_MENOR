&notebook_optical

!.............................................
! LIGHT ATTENUATION:
! z-profil of light attenuation is of the form:
! A1*exp(kpar1*z)+A2*exp(kpar2*z)

light_rat1=0.58       ! A1
light_att1=0.35       ! par1 (m)

light_rat2=0.42       ! A2=1-A1 
!light_att2=23.        ! par2 (m) 
 light_att2=15.        ! par2 (m) 
!texte250='../../../GLOBMED/BATHYMASK/symphonie-symtools-comodo.nc' ! lon-lat-dependent from file
! light_att2_val1=10. 
! light_att2_h1=0.    ! light_att2=light_att2_val1 if h<light_att2_h1
! light_att2_val2=20. 
! light_att2_h2=200.  ! light_att2=light_att2_val2 if h>light_att2_h2 , linear in between


!.............................................
! ALBEDO
!texte30='apel1987'
!texte30='br1982'
!texte30='br1986'
texte30='constant'
albedo_val=0.08   ! if texte30=='constant', give the value of the constant
!albedo_val=0.066 ! Maraldi et al, Ocean Science, 2013

/


 Explications:
 L'attenuation de la lumiere en fonction 
 de la profondeur comptée à partir de la sse 
 est parametree de la maniere suivante:
 
 A1*exp(z/L1)+A2*exp(z/L2)

 Avec 0<A1<1  0<A2<1 et A1+A2=1
 A1 est LIGHT_RAT1
 A2 est LIGHT_RAT2

 L1 est LIGHT_ATT1 exprime en metres
 L2 est LIGHT_ATT2 exprime en metres

Reference bibliographique associée:
Paulson C. A., and J. J. Simpson, 1977: Irradiance !measurements in the upper ocean. J. Phys. Oceanogr., 7, !952-956.
Maraldi et al, Ocean Science 2013

c...................................................
Informations generales:
 Déduit de la figure 9.44 de "Principles of Ocean Physics"
 par J.R. Apel montrant la fraction de l'irradiance de surface
 en fonction de la profondeur. Figure elle même extraite
 d'une étude de Jerlov dans Marine Optics en 1976.
      LIGHT_ATT= 7.  metres CAS III
      LIGHT_ATT= 11. metres CAS II
      LIGHT_ATT= 22. metres CAS I

Les eaux côtières turbides (sediment, matiere biologique...)
sont caracterisées par des valeurs encore plus faibles,  
eventuellement inferieures au metre.
c...................................................


Valeurs possibles:
Configuration intercomparaison projet opa-symphonique
(c.a.d. valeurs dans opa):
0.35            ! LIGHT_ATT1
0.58            ! LIGHT_RAT1
23.             ! LIGHT_ATT2
0.42            ! LIGHT_RAT2

Version originale du modele:
22.             ! LIGHT_ATT1
1.              ! LIGHT_RAT1
22.             ! LIGHT_ATT2
0.              ! LIGHT_RAT2
