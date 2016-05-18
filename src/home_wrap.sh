#!/bin/bash
FILES=../data/data_dump/trait_info*
for f in $FILES;
do
  ../stan/turnover_horseshoe variational \
    id=1 \
    random seed=420 \
    data file=$f \
    output file=../data/mcmc_out/turnover_hs_advi.csv &
  wait
done
FILES=../data/data_dump/trait_w_gaps*
for f in $FILES;
do
  ../stan/turnover_model variational \
    id=1 \
    random seed=420 \
    data file=$f \
    output file=../data/mcmc_out/turnover_md_advi.csv &
  wait
done
