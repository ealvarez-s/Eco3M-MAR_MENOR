#! /bin/bash

#conda activate belich-legos-r

INDIR3=/home/alvarez/EVA/Eco3M-MAR_MENOR/BGC_MAR_MENOR/GRAPHIQUES/RUN_test_48
OUTDIR=/home/alvarez/EVA/Eco3M-MAR_MENOR/BGC_MAR_MENOR/GRAPHIQUES/RUN_test_48
CODEDIR=/home/alvarez/EVA/Eco3M-MAR_MENOR/scripts_postproc

Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} z_eu ${OUTDIR}

Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} PAR ${OUTDIR}
Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} pH ${OUTDIR}

Rscript ${CODEDIR}/paste_forcing_variables_3D.R ${INDIR3} tem ${OUTDIR}
Rscript ${CODEDIR}/paste_forcing_variables_3D.R ${INDIR3} sal ${OUTDIR}

#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} zoonanoc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} zoomicroc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} zoomesoc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} synec ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} synen ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} synep ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} synechl ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nanoc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nanon ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nanop ${OUTDIR}   
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nanochl ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} diac ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} dian ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} diap ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} diachl ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} diasi ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} bactc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopn ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopp ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopchl ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopsi ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} lmopc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} lmopn ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} lmopp ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} lmopsi ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} modc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} modn ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} modp ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nitrate ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} ammonium ${OUTDIR}  
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} phosphate ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} silice ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} oxygen ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} odu ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} dic ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} talk ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} density ${OUTDIR}
