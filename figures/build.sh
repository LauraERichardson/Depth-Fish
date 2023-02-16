#!/bin/bash

set -ex

dir /input/

export CLOUD_RUN="TRUE"
#Rscript -e "R.version; installed.packages(fields = 'Version')[,'Version']"
Rscript pop_effect_figs.R # doesn't work with new int terms in s()
#Rscript depth_dens_plot.R
#Rscript slope_dens_plot.R
#Rscript proba_of_increase_plot.R
#Rscript proba_of_increase_slope.R
#Rscript island_ecoregion_depth_plot.R
#Rscript ternery_plot.R
#Rscript Table1.R
#Rscript TablesS7_S8.R
#Rscript summary_tables.R
#Rscript dens_prior_post_plot.R


if test -n "$(find . -maxdepth 1 -name '*.csv' -print -quit)"
then
    cp *.csv ${OUTPUTDIR}/.
fi


if test -n "$(find . -maxdepth 1 -name '*.txt' -print -quit)"
then
    cp *.txt ${OUTPUTDIR}/.
fi

if test -n "$(find . -maxdepth 1 -name '*.pdf' -print -quit)"
then
    cp *.pdf ${OUTPUTDIR}/.
fi

if test -n "$(find . -maxdepth 1 -name '*.html' -print -quit)"
then
    cp *.html ${OUTPUTDIR}/.
fi

if test -n "$(find . -maxdepth 1 -name '*.png' -print -quit)"
then
    cp *.png ${OUTPUTDIR}/.
fi
