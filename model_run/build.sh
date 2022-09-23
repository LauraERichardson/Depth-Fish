#!/bin/bash

set -ex

export CLOUD_RUN="TRUE"
Rscript ../data_prep/prepare_data.R && Rscript run_models.R