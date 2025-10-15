#! /bin/bash

#conda activate belich-legos-r

OUTDIR=/home/alvarez/EVA/Eco3M-MAR_MENOR/BGC_MAR_MENOR/LIST/
# CODEDIR=/home/alvarez/EVA/Eco3M-MAR_MENOR/meteo/

Rscript list_meteo_files.R

cp liste_ecmwf_dp2m ${OUTDIR}/liste_ecmwf_dp2m
cp liste_ecmwf_dp2m ${OUTDIR}/liste_ecmwf_ir
cp liste_ecmwf_dp2m ${OUTDIR}/liste_ecmwf_p0m
cp liste_ecmwf_dp2m ${OUTDIR}/liste_ecmwf_rain
cp liste_ecmwf_dp2m ${OUTDIR}/liste_ecmwf_ssr
cp liste_ecmwf_dp2m ${OUTDIR}/liste_ecmwf_t2m
cp liste_ecmwf_dp2m ${OUTDIR}/liste_ecmwf_u10m
cp liste_ecmwf_dp2m ${OUTDIR}/liste_ecmwf_v10m
