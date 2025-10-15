#! /bin/bash

#conda activate belich-legos-r

INDIR3=/home/alvarez/EVA/Eco3M-MAR_MENOR/OCEANCOLOUR_MED_BGC/
OUTDIR=/home/alvarez/EVA/Eco3M-MAR_MENOR/OCEANCOLOUR_MED_BGC/
CODEDIR=/home/alvarez/EVA/Eco3M-MAR_MENOR/OCEANCOLOUR_MED_BGC/

#Rscript ${CODEDIR}/get_Chla_multi.R ${INDIR3} A 37.7 -0.8
Rscript ${CODEDIR}/get_Chla_multi.R ${INDIR3} B 37.7 -0.8
#Rscript ${CODEDIR}/get_Chla_multi.R ${INDIR3} C 37.7 -0.8

#Rscript ${CODEDIR}/get_Chla_multi.R ${INDIR3} A1-6 37.7 -0.8
#Rscript ${CODEDIR}/get_Chla_multi.R ${INDIR3} B1-12 37.7 -0.8

Rscript ${CODEDIR}/get_SPM_sentinel.R ${INDIR3} B 37.7 -0.8
