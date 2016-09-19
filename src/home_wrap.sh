#!/bin/bash
../stan/turnover_basic variational \
  id=1 \
  init=0 \
  data file=../data/data_dump/trait_info.data.R \
  output file=../data/mcmc_out/turnover_basic_advi.csv &
wait
#../stan/turnover_full variational \
#  id=1 \
#  init=0 \
#  data file=../data/data_dump/trait_w_gaps.data.R \
#  output file=../data/mcmc_out/turnover_full_advi.csv &
wait
