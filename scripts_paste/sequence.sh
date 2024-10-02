#! /bin/bash

#conda activate to_run_R

INDIR3=/tmpdir/alvarez/MAR_MENOR/BGC_MAR_MENOR/GRAPHIQUES/RUN_reference_10days
OUTDIR=/tmpdir/alvarez/MAR_MENOR/BGC_MAR_MENOR/GRAPHIQUES/RUN_reference_10days

Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_forcing_variables_2D.R ${INDIR3} ssr ${OUTDIR}

#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} zoonanoc ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} zoomicroc ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} zoomesoc ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} synec ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} synen ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} synep ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} synechl ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} nanoc ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} nanon ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} nanop ${OUTDIR}   
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} nanochl ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} diac ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} dian ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} diap ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} diachl ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} diasi ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} bactc ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} smopc ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} smopn ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} smopp ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} smopchl ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} smopsi ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} lmopc ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} lmopn ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} lmopp ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} lmopsi ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} modc ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} modn ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} modp ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} nitrate ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} ammonium ${OUTDIR}  
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} phosphate ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} silice ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} oxygen ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} odu ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} dic ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} talk ${OUTDIR}
#Rscript /tmpdir/alvarez/MAR_MENOR/scripts_paste/paste_state_variables_3D.R ${INDIR3} density ${OUTDIR}
