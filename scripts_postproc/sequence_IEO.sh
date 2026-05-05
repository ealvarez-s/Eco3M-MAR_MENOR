#! /bin/bash

#conda activate seamless-bb-r

INDIR3=/home/alvarez/EVA/Eco3M-MAR_MENOR/BGC_MAR_MENOR/GRAPHIQUES/RUN_test_T35_T37_7y_SWAT_EXPpH_50_15
OUTDIR=/home/alvarez/EVA/Eco3M-MAR_MENOR/BGC_MAR_MENOR/GRAPHIQUES/RUN_test_T35_T37_7y_SWAT_EXPpH_50_15
CODEDIR=/home/alvarez/EVA/Eco3M-MAR_MENOR/scripts_postproc/low-res

#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} z_eu ${OUTDIR}

# Benthic fluxes
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} Pefflux2d ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} Siefflux2d ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} NO3efflux2d ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} NH4efflux2d ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} DICefflux2d ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} 02influx2d ${OUTDIR}
# Benthic deposition
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} CDepo ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} NDepo ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} PDepo ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} SiDepo ${OUTDIR}
# Benthic mineralization
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} CMinB ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} NMinB ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} PMinB ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} SiMinB ${OUTDIR}
# Benthic nitrif/denitrif
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} NitrificationB ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} DenitrificationB ${OUTDIR}
# Benthic pools
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} CBFDetB ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} CBSDetB ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} NBDetB ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} PBDetB ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} SiBDetB ${OUTDIR}
# Air-sea fluxes
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} O2Flux ${OUTDIR}
#Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} CO2_Flux ${OUTDIR}
# Diagnostics
#Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} PAR ${OUTDIR}      # instantaneous !!change
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} PAR ${OUTDIR}     # instantaneous
Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} pH ${OUTDIR}       # averaged
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} pH ${OUTDIR}      # instantaneous
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} density ${OUTDIR} # instantaneous
# CSYS
Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} pHT ${OUTDIR}     # averaged
Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} iCO2 ${OUTDIR}    # averaged
Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} iHCO3 ${OUTDIR}   # averaged
Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} iCO3 ${OUTDIR}    # averaged
Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} pCO2aq ${OUTDIR}  # averaged
Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} omegaCa ${OUTDIR} # averaged
Rscript ${CODEDIR}/paste_supp_variables_2D.R ${INDIR3} pCO2 ${OUTDIR}    # instantaneous
# Water column 3D fluxes
#Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} gpp3d ${OUTDIR}   # averaged
#Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} resp3d ${OUTDIR}   # averaged
#Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} exuc3d ${OUTDIR}   # averaged
#Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} nitrif3d ${OUTDIR}  # averaged
#Rscript ${CODEDIR}/paste_supp_variables_3D.R ${INDIR3} uptnit3d ${OUTDIR}  # averaged
# Phys
Rscript ${CODEDIR}/paste_forcing_variables_3D.R ${INDIR3} tem ${OUTDIR}
Rscript ${CODEDIR}/paste_forcing_variables_3D.R ${INDIR3} sal ${OUTDIR}
# BGC state
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} zoonanoc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} zoomicroc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} zoomesoc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} synec ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} synen ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} synep ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} synechl ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nanoc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nanon ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nanop ${OUTDIR}   
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nanochl ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} diac ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} dian ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} diap ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} diachl ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} diasi ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} bactc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopn ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopp ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopchl ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} smopsi ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} lmopc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} lmopn ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} lmopp ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} lmopsi ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} modc ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} modn ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} modp ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} nitrate ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} ammonium ${OUTDIR}  
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} phosphate ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} silice ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} oxygen ${OUTDIR}
#Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} odu ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} dic ${OUTDIR}
Rscript ${CODEDIR}/paste_state_variables_3D.R ${INDIR3} talk ${OUTDIR}
