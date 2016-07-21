#!/bin/bash
../stan/turnover_horseshoe variational \
  id=1 \
  init=0 \
  data file=../data/data_dump/trait_info.data.R \
  output file=../data/mcmc_out/turnover_hs_advi.csv &
wait
../stan/turnover_model variational \
  id=1 \
  init=0 \
  data file=../data/data_dump/trait_w_gaps.data.R \
  output file=../data/mcmc_out/turnover_md_advi.csv &
wait
